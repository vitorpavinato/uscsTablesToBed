"""
This program return BED file(s) from UCSC annotation table.
It returns the exons and introns from the UCSC tables. It also returns the
exon-introns and intron-exons boundaries from the UCSC tables.

The input table should be downloaded from UCSC Table Browser as a tsv file
with the folling fields selected when you define the output format in the
'Retrieve and Display Data' section:

- name
- chrom
- strand
- exonStarts
- exonEnds
- score.
"""

import sys
import argparse
from utils import *


def parse_chromosome(chrom):
    """
    Parse chromsomes type
    """
    if chrom.startswith("chr"):
        return chrom
    else:
        return f"chr{chrom}"


def parseargs():
    """
    Function defines command-line parsing arguments.
    """
    parser = argparse.ArgumentParser("python ucscTableToBed.py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", help="Input UCSC Table file", dest="inputfile", required=True, type=str)
    parser.add_argument("-o", help="Output BED file", dest="outputfile", default=None, type=str)
    parser.add_argument("-j", help="Output BED file containing exon/intron junctions", dest="outputjfile", default=None, type=str)
    parser.add_argument("-c", help="List of chroms to process", dest="chrom_list", nargs="+", required=True, type=parse_chromosome)
    parser.add_argument("-e", help="Exon side in bp", dest="exon_side", default=3, type=int)
    parser.add_argument("-s", help="Intron side in bp", dest="intron_side", default=10, type=int)
    return parser


def main(argv) -> None:
    """
    Main function.
    """

    parser = parseargs()
    if argv[-1] == "":
        argv = argv[0:-1]
    args = parser.parse_args(argv)

    # Define the input and output files
    inputfile = args.inputfile
    outputfile = args.outputfile
    outputjfile = args.outputjfile
    chrom_list = args.chrom_list
    exon_side = args.exon_side
    intron_side = args.intron_side

    # Check the input and output files
    inputfile = check_input_file(inputfile)
    outputfile, outputjfile = check_output_file_name(inputfile, outputfile, outputjfile)

    # Process the wig file to csv
    result = process_table_to_bed(
        inputfile=inputfile,
        chrom_list=chrom_list,
        outputfile=outputfile,
        outputjfile=outputjfile,
        exon_side=exon_side,
        intron_side=intron_side,
    )

    print(result)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        main(["-h"])
    else:
        main(sys.argv[1:])
