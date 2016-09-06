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
READ_BREAK_type = namedtuple("ReadBreakPoint",
                             ["qname", "chr3p", "pos3p", "insert_strand",
                              "chr5p", "pos5p", "query_gap"])


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

    def __init__(self,
                 inputbam,
                 region,
                 log_reasons=False,
                 min_mapQ=20.0,
                 max_locus_size=1e5):
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

        self._min_mapQ = min_mapQ
        self._max_locus_size = max_locus_size

        self._break_support = {}
        self._twin_primed_reads = set()

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

    def _find_read_breakpoints(self, alns_by_query_pos, seqLen, aln):
        # Previous alignment segment
        import logging as log
        prev_alns = None
        aln_breaks = []

        is_twin_primed = False
        for _, alns in alns_by_query_pos:
            if prev_alns is not None and self.on_target_region(
                    prev_alns) and self.on_target_region(alns):
                if prev_alns.is_reverse != alns.is_reverse:
                    self._twin_primed_reads.add(aln.query_name)
                    is_twin_primed = True

            prev_alns = alns

        prev_alns = None
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
                    # 3' end of the site on the chromosome where L1 was inserted
                    if prev_alns.is_reverse:
                        insert3p = prev_alns.reference_start
                        insert_query_end = seqLen - prev_alns.query_alignment_start
                    else:
                        insert3p = prev_alns.reference_end
                        insert_query_end = prev_alns.query_alignment_end

                    insert_chr = self._inbam.getrname(prev_alns.reference_id)

                    # Position of the dangling 5' atom that is attached to the insertion site 3' atom.
                    if alns.is_reverse:
                        source5p = alns.reference_end
                        source_query_start = seqLen - alns.query_alignment_end
                    else:
                        source5p = alns.reference_start
                        source_query_start = alns.query_alignment_start

                    insert_strand = '+' if alns.is_reverse == prev_alns.is_reverse else '-'

                    query_gap = source_query_start - insert_query_end
                    d = READ_BREAK_type(aln.query_name, insert_chr, insert3p,
                                        insert_strand, self._chrom, source5p,
                                        query_gap)
                    aln_breaks.append(d)

                    #log.debug("{}  {}:{:d}] {} [{}:{:d} (gap {:d})".format(
                    #    aln.query_name, insert_chr, insert5p, insert_strand,
                    #    self._chrom, source3p, query_gap))
                    pass
            else:
                # We're on alternative region
                if prev_alns is None:
                    # First segment is alternative (Missing the target region mapping)
                    pass
                elif prev_alns.reference_id == alns.reference_id:
                    # Remain on the same alternative region, probably jumping to different strand due circularisation
                    # TODO: Could try to infer Target Site duplication here!!
                    # TODO: Could try to infer restriction cut sites here!!
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

                    #log.debug("{}  {}:{:d}] {} [{}:{:d} (gap {:d})".format(
                    #    aln.query_name, self._chrom, source3p, insert_strand,
                    #    insert_chr, insert5p, query_gap))
                    d = READ_BREAK_type(aln.query_name, self._chrom, source3p,
                                        insert_strand, insert_chr, insert5p,
                                        query_gap)
                    aln_breaks.append(d)
                    pass
            prev_alns = alns
        return aln_breaks

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

        aln_breaks = self._find_read_breakpoints(alns_by_query_pos, seqLen,
                                                 aln)
        if len(aln_breaks) > 2:
            log.debug("Loads of breakpoints on {} {}".format(aln.query_name,
                                                             str(aln_breaks)))
        else:
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
                if not aln.is_supplementary:  # TODO: Get rid of this
                    self.add_break_info(aln)
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

    def call_insertions(self, gap_width=3000):
        """
        """
        import logging as log
        import itertools as it
        from collections import defaultdict

        bts_from_tgt = defaultdict(list)
        bts_to_tgt = defaultdict(list)
        bts_all = defaultdict(list)
        for brk in it.chain(*self._break_support.itervalues()):
            if brk.chr3p == self._chrom and brk.pos3p > self._begin and brk.pos3p < self._end:
                # 3' end of the breakpoint is in the target hence the insertion comes from target
                #bts_to_tgt[(brk.chr5p, brk.insert_strand)].append(brk.pos5p)
                bts_to_tgt[brk.chr5p].append((brk.pos5p, brk.qname))
                bts_all[brk.chr5p].append((brk.pos5p, brk.qname))
            else:
                bts_from_tgt[brk.chr3p].append((brk.pos3p, brk.qname))
                bts_all[brk.chr3p].append((brk.pos3p, brk.qname))

        # Find clusters of breakpoints

        def _cluster_breakpoints(breakpoints):
            from collections import defaultdict, Counter
            ret = defaultdict(list)
            for chrom, breaks in breakpoints.iteritems():
                breaks.sort()
                cluster_start, _ = breaks[0]
                point = cluster_start
                breaks_in_cluster = 0
                break_count = Counter()
                break_reads = set()
                for next_break_pos, read in breaks:
                    if next_break_pos <= point + gap_width:  # More breakpoints in this cluster
                        point = next_break_pos
                        breaks_in_cluster += 1
                        break_count[point] += 1
                        break_reads.add(read)
                    else:
                        log.debug("Cluster {}:{}-{}  len={}#  {}".format(
                            chrom, cluster_start, point, point - cluster_start,
                            breaks_in_cluster))
                        ret[chrom].append(
                            (cluster_start, point, breaks_in_cluster,
                             break_count.most_common(1)[0][0], break_reads))
                        breaks_in_cluster = 1
                        cluster_start = next_break_pos
                        point = next_break_pos

                        break_count = Counter()
                        break_count[point] += 1
                        break_reads = set([read])

                if len(ret[chrom]) == 0 or ret[chrom][-1][0] != cluster_start:
                    # Last breakpoint cluster
                    ret[chrom].append((cluster_start, point, breaks_in_cluster,
                                       break_count.most_common(1)[0][
                                           0], break_reads))

                    log.debug("Cluster {}:{}-{}  len={}#  {}".format(
                        chrom, cluster_start, point, point - cluster_start,
                        breaks_in_cluster))

            return ret

        log.debug("Breakpoints from target")
        from_tgt_clust = _cluster_breakpoints(bts_from_tgt)
        log.debug("Breakpoints To target")
        to_tgt_clust = _cluster_breakpoints(bts_to_tgt)

        log.debug("All breakpoints")
        all_clusts = _cluster_breakpoints(bts_all)

        rclust = {}
        for chrom, clusters in all_clusts.iteritems():
            rclusts = []
            #twin_primed = True
            for clust in clusters:
                if chrom == self._chrom and clust[3] < self._end and clust[
                        3] >= self._begin:
                    continue

                twin_primed = clust[4].isdisjoint(
                    self._twin_primed_reads) == False

                rclusts.append(clust[:4] + (twin_primed, clust[4]))
            rclust[chrom] = rclusts

        all_clusts = rclust
        return dict(from_tgt_clust), dict(to_tgt_clust), dict(all_clusts)
        pass


if __name__ == '__main__':
    args = main()
    cmd = LDIPCR(args.input[0],
                 args.region,
                 log_reasons=args.reasons is not None)

    cmd.write_filtered(args.output, args.filtered)
    cmd.flush_reasons(args.reasons)
    f, t, all_clusters = cmd.call_insertions()

    o = open("insertion_sites.bed", "w")
    STRAND = "."
    for i, (CHR, CLUSTERS) in enumerate(all_clusters.iteritems()):
        o.writelines("{}\t{}\t{}\t{}\t{}\t{}\n".format(
            CHR, MOST_LIKELY, MOST_LIKELY + 1, STRAND, "TwinPrimed" if
            TWINPRIMED else "UnknownPriming", SCORE)
                     for B, E, SCORE, MOST_LIKELY, TWINPRIMED, _ in CLUSTERS)
    o.close()

    import itertools as it
    open("breaks.tsv", "w").writelines(
        "\t".join(map(str, x) + [str(x.qname in cmd._twin_primed_reads)]) +
        "\n" for x in it.chain(*cmd._break_support.itervalues()))
