# Software for calling insertions from long read sequenced LDI-PCR reads

Implements method for calling genomic insertions from Long-Range Inverse PCR sequences as described in "“Detection and Analysis of Somatic LINE-1 Insertions in Colorectal Cancer by Long-distance Inverse-PCR and Nanopore Sequencing” by Pradhan et al. (Submitted)"


	usage: LDI_PCR_call.py [-h] [-o OUTPUT] [-f FILTERED] -r REGION [-R REASONS]
                       [-V [VERBOSE]]
                       input

	Implementing Long-Range Inverse PCR variant calling method for long read
	sequencing (e.g. Oxford Nanopore) The overall strategy is: 10. Consider
	supplementary alignments with mapping quality of at least 20 20. Discard all
	reads that do not have an alignment in the user given primer region 30. Chain
	supplementary alignments consequtive in both read and genomic coordinates with
	approximately equal distances 40. Cluster the reads with compatible/similar
	breakpoints and call the variants (i.e. breakpoints)

	positional arguments:
	input                 Input bam file for filtering reads [default:None]

	optional arguments:
	-h, --help            show this help message and exit
	-o OUTPUT, --output OUTPUT
                          Output bam for LDI-PCR reads [default:LDI_PCR.bam]
	-f FILTERED, --filtered FILTERED
                        Output bam for discarded reads [default:None]
	-r REGION, --region REGION
                        Region of the LDI-PCR primers required for mapping
                        [default:None]
	-R REASONS, --reasons REASONS
                        For reads mapping to the target region, write read
                        names, filtering decisions and reasons to this file
                        [default:None]
	-V [VERBOSE], --verbose [VERBOSE]
                        Be more verbose with output [and log to a file]
                        [default:False]
