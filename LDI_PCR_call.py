#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Created at Mon Mar 14 15:05:11 2016 by Kimmo Palin <kpalin@merit.ltdk.helsinki.fi>
"""

from collections import namedtuple
SAtype = namedtuple("Alignment",
                    ["rname", "pos", "strand", "CIGAR", "mapQ", "NM"])


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description=
        "Filter a bam file for Long-Distance-Inverse PCR reads. These reads should have something mapping on given region and exactly two other fragments (Mapping close by on the same chromosome. Not implemented.)")

    parser.add_argument(
        "input",
        nargs=1,
        help="Input bam file for filtering reads [default:%(default)s]", )
    parser.add_argument(
        "-o",
        "--output",
        default="LDI_PCR.bam",
        help="Output bam for LDI-PCR reads [default:%(default)s]", )

    parser.add_argument(
        "-f",
        "--filtered",
        default=None,
        help="Output bam for discarded reads [default:%(default)s]", )

    parser.add_argument(
        "-r",
        "--region",
        required=True,
        help=
        "Region of the LDI-PCR primers required for mapping [default:%(default)s]")

    parser.add_argument(
        "-R",
        "--reasons",
        default=None,
        help=
        "For reads mapping to the target region, write read names, filtering decisions and reasons to this file [default:%(default)s]")

    parser.add_argument(
        "-V",
        "--verbose",
        default=False,
        const=True,
        nargs="?",
        help=
        "Be more verbose with output [and log to a file] [default:%(default)s]")

    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s:%(funcName)s:%(levelname)s:%(message)s')
        if args.verbose is not True:
            log_file_handler = logging.FileHandler(args.verbose)
            log_file_handler.setFormatter(logging.getLogger().handlers[
                0].formatter)
            logging.getLogger().addHandler(log_file_handler)

    logging.info(str(args))
    return args


def parse_supplementary_alignments(SAtag):
    r = []
    for SA in SAtag.rstrip(";").split(";"):
        p = SA.split(",")
        p[1] = long(p[1])
        p[4] = int(p[4])
        p[5] = int(p[5])
        r.append(SAtype(*p))
    return r


class LDIPCR(object):
    """Manage LDI PCR bams
    """

    def __init__(self, inputbam, region, log_reasons=False):
        """

        Arguments:
        - `inputbam`: Bam file of the input reads
        - `region`:   Region of the LDI-PCR primers
        """
        self._inputbam = inputbam
        self._region = region

        self._log_reasons = log_reasons
        self._reasons = dict()

        import pysam
        self._inbam = pysam.Samfile(self._inputbam)

        chrom, pos = region.split(":")
        self._chrom = chrom

        self._begin, self._end = map(int, pos.split("-"))

    def is_valid_LDI_read(self, aln):
        """Return True if, and only if, the input alignment is valid LDI-PCR read, that is it maps to given target region and and has exactly two other alignment fragments mapping to one chromosome.

        Arguments:
        - `read`:
        """

        if aln.has_tag("SA"):
            SAtag = aln.get_tag("SA")
        else:
            self.reason(aln.qname, "No supplementary alignments")
            return False  # No supplementary alignments

        aln_positions = parse_supplementary_alignments(SAtag)

        this_position = SAtype(
            self._inbam.getrname(aln.reference_id), aln.pos, "-"
            if aln.is_reverse else '+', aln.cigarstring, aln.mapq,
            aln.get_tag("NM"))

        aln_positions.append(this_position)

        target_positions = [x
                            for x in aln_positions
                            if x.rname == self._chrom and x.pos <= self._end
                            and x.pos >= self._begin]

        if len(target_positions) == 0:
            self.reason(aln.qname, "No mapping to target location")
            return False  # No mapping to target location

        if len(target_positions) > 2:
            self.reason(aln.qname,
                        "More than two fragments mapping to target location")
            return False

        alt_positions = [x
                         for x in aln_positions
                         if not (x.rname == self._chrom and x.pos <= self._end
                                 and x.pos >= self._begin)]

        if len(alt_positions) == 0 and len(target_positions) == 2:
            self.reason(aln.qname, "Native locus read with a break.")
            return False

        if len(alt_positions) != 2:
            self.reason(
                aln.qname,
                "No exactly two alternative alignments outside target region. Had %d."
                % (len(alt_positions)))

            # Not exactly two alternative alignments outside target region.
            return False

        if alt_positions[0].rname != alt_positions[1].rname:
            self.reason(aln.qname,
                        "Fragments mapping to different chromosomes.")
            # Mapping to different chromosomes
            return False

        self.reason(aln.qname, "Valid LDI-PCR read", fail=False)

        return True

    def iter_LDI_reads(self, ):
        """Return iterator (generator) over valid LDI-PCR alignments
        """
        alnIt = self._inbam.fetch(region=self._region)

        for aln in alnIt:
            if self.is_valid_LDI_read(aln):
                yield aln

    def flush_reasons(self, reasons):
        """Log filtering decisions to given file
        
        Arguments:
        - `reasons`:
        """
        if reasons is not None and self._log_reasons:
            reasons_file = open(reasons, "w")
            reasons_file.write("read\tfiltered\treason\n")
            reasons_file.writelines(self._reasons.itervalues())
            reasons_file.close()

    def reason(self, read_name, reason, fail=True):
        """Log the reason for failure if asked to
        
        Arguments:
        - `read_name`:
        - `reason`:
        - `fail`:
        """
        if self._log_reasons:
            s = "{}\t{}\t{}\n".format(read_name, fail, reason)
            if read_name in self._reasons:
                assert self._reasons[read_name] == s
            self._reasons[read_name] = s

    def write_filtered(self, outputbam, filteredbam=None):
        """Write filtered reads to outputbam

        Arguments:
        - `outputbam`:
        """
        import pysam
        import logging as log

        valid_reads = set()

        log.info("Finding valid LDI-PCR reads")
        for aln in self.iter_LDI_reads():
            valid_reads.add(aln.qname)

        log.info("Found %d LDI-PCR reads" % (len(valid_reads)))

        outbam = pysam.AlignmentFile(outputbam, "wb", template=self._inbam)
        if filteredbam is not None:
            filtered_bam = pysam.AlignmentFile(
                filteredbam, "wb", template=self._inbam)
            log.info("Outputting filtered reads to %s" % (filteredbam))
        else:
            filtered_bam = None

        out_count, filter_count = 0, 0
        for aln in self._inbam.fetch():
            if aln.qname in valid_reads:
                outbam.write(aln)
                out_count += 1
            else:
                if filtered_bam is not None:
                    filtered_bam.write(aln)
                filter_count += 1

        log.info("Wrote %d alignments to %s and filtered %d alignments" %
                 (out_count, outputbam, filter_count))

        outbam.close()
        if filtered_bam is not None:
            filtered_bam.close()


if __name__ == '__main__':
    args = main()
    cmd = LDIPCR(args.input[0],
                 args.region,
                 log_reasons=args.reasons is not None)

    cmd.write_filtered(args.output, args.filtered)
    cmd.flush_reasons(args.reasons)
