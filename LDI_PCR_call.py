#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Created at Mon Mar 14 15:05:11 2016 by Kimmo Palin <kpalin@merit.ltdk.helsinki.fi>

Implementing Long-Range Inverse PCR variant calling method for long read sequencing (e.g. Oxford Nanopore)

The overall strategy is:

10. Consider supplementary alignments with mapping quality of at least 10
20. Discard all reads that do not have an alignment in the user given primer region
30. Chain supplementary alignments consequtive in both read and genomic coordinates with approximately equal distances
40. Cluster the reads with compatible/similar breakpoints and call the variants (i.e. breakpoints

"""

from collections import namedtuple
SAtype = namedtuple("Alignment",
                    ["rname", "pos", "strand", "CIGAR", "mapQ", "NM", "epos"])


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description=
        """Implementing Long-Range Inverse PCR variant calling method for long read sequencing (e.g. Oxford Nanopore)

The overall strategy is:

10. Consider supplementary alignments with mapping quality of at least 10
20. Discard all reads that do not have an alignment in the user given primer region
30. Chain supplementary alignments consequtive in both read and genomic coordinates with approximately equal distances
40. Cluster the reads with compatible/similar breakpoints and call the variants (i.e. breakpoints)
""")

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
            level=logging.DEBUG,
            format='%(asctime)s:%(funcName)s:%(levelname)s:%(message)s')
        if args.verbose is not True:
            log_file_handler = logging.FileHandler(args.verbose)
            log_file_handler.setFormatter(logging.getLogger().handlers[
                0].formatter)
            logging.getLogger().addHandler(log_file_handler)

    logging.info(str(args))
    return args


def parse_supplementary_alignments(SAtag):
    "Format BAM SA tag to SAtype tupple"
    r = []
    import pysam
    for SA in SAtag.rstrip(";").split(";"):
        p = SA.split(",")
        p[1] = long(p[1])
        p[4] = int(p[4])
        p[5] = int(p[5])

        a = pysam.AlignedSegment()
        a.flag = 0
        a.flag = 16 if SA[2] == "-" else 0

        a.cigarstring = p[3]
        a.reference_start = p[1]
        a.seq = "N" * a.infer_query_length(always=True)
        a.mapq = p[4]

        p.append(a.reference_end)
        SA = SAtype(*p)

        r.append(SA)
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

        self._min_mapQ = 10.0
        self._max_locus_size = 1e6

        self._break_support = {}

    def is_valid_LDI_read(self, aln):
        """Return True if, and only if, the input alignment is valid LDI-PCR read, 


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
            aln.get_tag("NM"), aln.reference_end)

        aln_positions.append(this_position)

        # 10. Consider supplementary alignments with mapping quality of at least 10 (self._min_mapQ)

        aln_positions = [x for x in aln_positions if x.mapQ >= self._min_mapQ]

        if len(aln_positions) <= 1:
            self.reason(aln.qname, "Low quality alignments")
            return False

        # 20. Discard all reads that do not have an alignment in the user given primer region
        # TODO: Be more exact on the regions. Now looking only beginning of the alignment
        target_positions = [x
                            for x in aln_positions
                            if x.rname == self._chrom and x.pos <= self._end
                            and x.pos >= self._begin]

        if len(target_positions) == 0:
            self.reason(aln.qname, "No mapping to target location")
            return False  # No mapping to target location

        alt_positions = [x
                         for x in aln_positions
                         if not (x.rname == self._chrom and x.pos <= self._end
                                 and x.pos >= self._begin)]

        if len(alt_positions) == 0 and len(target_positions) >= 2:
            self.reason(aln.qname, "Native locus read with breaks.")
            return False

        if len(set(x.rname for x in alt_positions)) > 1:
            self.reason(aln.qname,
                        "Fragments mapping to different chromosomes.")
            return False

        def _max_distance(pos):
            if len(pos) < 2:
                return 0
            pos = [x.pos for x in pos]
            pos.sort()
            distances = [np - p for p, np in zip(pos, pos[1:])]
            return max(distances)

        if _max_distance(alt_positions) > self._max_locus_size:
            self.reason(aln.qname, "Fragments mapping far from each other.")
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

    def on_target_region(self, aln):
        """
        
        Arguments:
        - `aln`:
        """

        on_target = True
        on_target &= self._inbam.getrname(aln.reference_id) == self._chrom
        overlap_bp = aln.get_overlap(self._begin, self._end)
        on_target &= (overlap_bp > 0)

        return on_target

    def add_break_info(self, aln):
        """Accumulate breakpoint info from this alignment
        
        Arguments:
        - `aln`:
        """
        import logging as log

        import pysam
        SAtag = aln.get_tag("SA")

        aln_positions = parse_supplementary_alignments(SAtag)

        this_position = SAtype(
            self._inbam.getrname(aln.reference_id), aln.pos, "-"
            if aln.is_reverse else '+', aln.cigarstring, aln.mapq,
            aln.get_tag("NM"), aln.reference_end)

        aln_positions.append(this_position)

        # 10. Consider supplementary alignments with mapping quality of at least 10 (self._min_mapQ)

        aln_positions = [x for x in aln_positions if x.mapQ >= self._min_mapQ]

        aln_segments = []

        for SA in aln_positions:
            a = pysam.AlignedSegment()
            a.reference_id = self._inbam.gettid(SA.rname)
            a.flag = 0
            a.flag = 16 if SA.strand == "-" else 0
            a.seq = aln.seq
            a.cigarstring = SA.CIGAR
            a.reference_start = SA.pos
            a.mapq = SA.mapQ

            if a.mapq >= self._min_mapQ:
                aln_segments.append(a)

        alns_by_query_pos = []
        seqLen = aln.infer_query_length()
        for alns in aln_segments:
            if alns.is_reverse:
                sStart = seqLen - alns.qend
                sEnd = seqLen - alns.qstart
                x = sEnd
                dx = sStart - sEnd
            else:
                sStart, sEnd = alns.qstart, alns.qend
                x = sStart
                dx = sEnd - sStart
            alns_by_query_pos.append((x, alns))

        alns_by_query_pos.sort()

        # Previous alignment segment
        prev_alns = None
        aln_breaks = []
        for _, alns in alns_by_query_pos:
            source3p, source5p, source_query_end, source_query_start = [None
                                                                        ] * 4
            insert_chr, insert3p, insert5p, insert_query_end, insert_query_start = [
                None
            ] * 5
            query_gap = None
            if self.on_target_region(alns):
                if prev_alns is None:
                    # First segment on target
                    pass
                elif self.on_target_region(prev_alns):
                    # Remain on target region
                    pass
                else:
                    # Break from alternative to target region !!REPORT""
                    #print "From  alt:", prev_alns
                    #print "To target:", alns
                    if prev_alns.is_reverse:
                        insert5p = prev_alns.reference_start
                        insert_query_end = seqLen - prev_alns.query_alignment_start
                    else:
                        insert5p = prev_alns.reference_end
                        insert_query_end = prev_alns.query_alignment_end

                    insert_chr = self._inbam.getrname(prev_alns.reference_id)

                    if alns.is_reverse:
                        source3p = alns.reference_end
                        source_query_start = seqLen - alns.query_alignment_end
                    else:
                        source3p = alns.reference_start
                        source_query_start = alns.query_alignment_start

                    insert_strand = '+' if alns.is_reverse == prev_alns.is_reverse else '-'

                    query_gap = source_query_start - insert_query_end
                    d = (aln.query_name, insert_chr, insert5p, insert_strand,
                         self._chrom, source3p, query_gap)
                    aln_breaks.append(d)

                    log.debug("{}  {}:{:d}] {} [{}:{:d} (gap {:d})".format(
                        aln.query_name, insert_chr, insert5p, insert_strand,
                        self._chrom, source3p, query_gap))
                    pass
            else:
                # We're on alternative region
                if prev_alns is None:
                    # First segment is alternative (Missing the target region mapping)
                    pass
                elif prev_alns.reference_id == alns.reference_id:
                    # Remain on the same alternative region, probably jumping to different strand due circularisation
                    pass
                else:
                    assert self.on_target_region(prev_alns)
                    # Break from target to alternative region !! REPORT!!
                    #print "From  target:", prev_alns
                    #print "To       alt:", alns
                    if prev_alns.is_reverse:
                        source3p = prev_alns.reference_start
                        source_query_end = seqLen - prev_alns.query_alignment_start
                    else:
                        source3p = prev_alns.reference_end
                        source_query_end = prev_alns.query_alignment_end

                    insert_chr = self._inbam.getrname(alns.reference_id)

                    if alns.is_reverse:
                        insert5p = alns.reference_end
                        insert_query_start = seqLen - alns.query_alignment_end
                    else:
                        insert5p = alns.reference_start
                        insert_query_start = alns.query_alignment_start

                    insert_strand = '+' if alns.is_reverse == prev_alns.is_reverse else '-'

                    query_gap = source_query_end - insert_query_start

                    log.debug("{}  {}:{:d}] {} [{}:{:d} (gap {:d})".format(
                        aln.query_name, self._chrom, source3p, insert_strand,
                        insert_chr, insert5p, query_gap))
                    d = (aln.query_name, self._chrom, source3p, insert_strand,
                         insert_chr, insert5p, query_gap)
                    aln_breaks.append(d)
                    pass
            prev_alns = alns

        pass
        self._break_support[aln.query_name] = aln_breaks

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
                if not aln.is_supplementary:
                    self.add_break_info(aln)
                #outbam.write(aln)
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
    import itertools as it
    open("breaks.tsv", "w").writelines(
        "\t".join(map(str, x)) + "\n"
        for x in it.chain(*cmd._break_support.itervalues()))
