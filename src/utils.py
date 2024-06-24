"""
Auxiliary functions for the ucscTableToBed.py
"""

import os
from pathlib import Path
from typing import Tuple


# Check if the input file
def check_input_file(inputfile: str) -> str:
    """
    Simple function to check if the input file was 
    provided by the user and if it exists;
    """

    # Check if input files are provided
    if inputfile == "" or inputfile is None:
        raise ValueError("file must be provided")

    # Check if the file exists
    if not os.path.exists(inputfile):
        raise ValueError("file does not exist")

    return inputfile


# Check if the output file
def check_output_file_name(
    inputfile: str,
    outputfile: str = None,
    outputjfile: str = None
) -> Tuple[str, str]:
    """
    This function checks if the output file name is provided.
    If not, it will use the input file name with the .bed extension
    """

    if outputfile == "" or outputfile is None:
        # Get the path and filename
        path, filename = os.path.split(inputfile)

        # Find the last occurrence of "/" in the path
        last_slash_index = path.rfind("/")

        # Extract the part before the last "/"
        base_path = path[:last_slash_index]

        # Append "tables" to the base path: the output
        # file will be in the "tables" folder
        outputfile_path = base_path + "/bedfiles"

        # Check if the outputfile path exists; create it if it doesn't
        if not Path(outputfile_path).exists():
            Path(outputfile_path).mkdir(parents=True)

        # Define the basename for the outputs from the filename
        basename, _ = filename.strip().split(".txt")

        # Define the output file with a path
        outputfile = outputfile_path + "/" + basename + ".bed"

        # Define the output junctions file with a path
        outputjfile = outputfile_path + "/" + basename + "junctions.bed"

    else:
        # Get the path and filename
        outputfile_path, filename = os.path.split(outputfile)

        # Check if the outputfile path exists; create it if it doesn't
        if not Path(outputfile_path).exists():
            Path(outputfile_path).mkdir(parents=True)

    return outputfile, outputjfile


def process_table_to_bed(
    inputfile: str, chrom_list: list[str],
    outputfile: str, outputjfile: str,
    exon_side: int, intron_side: int,
) -> None:
    """
    Function implementation for bed_from_ucsc_tables.
    """

    # Parameters
    exside = exon_side
    intside = intron_side

    # Lists for outputs
    coding_list = []
    junctions_list = []

    # for line in open(inputfile,'r', encoding='utf-8'):
    with open(inputfile, 'r', encoding='utf-8') as infile:
        for linen, line in enumerate(infile):
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue

            # Skip empty lines
            if len(line) <= 1:
                continue
            
            # temp line list
            temp = line.strip().split()

            # Line fields
            name = temp[0]
            chrom = temp[1]
            strand = temp[2]
            direction = 'f' if strand == '+' else 'r'
            score = temp[5]

            # Process only chroms in the list:
            # This make the inclusion of any chromosome explicit.
            if any(x == chrom for x in chrom_list):
                try:
                    exonstarts = list(map(int, temp[3].split(sep=',')[:-1]))
                except:
                    print(temp)
                    pass
                try:
                    exonends = list(map(int, temp[4].split(sep=',')[:-1]))
                except:
                    print(temp)
                    pass
                
                # Get number of exons and introns in the transcript (line):
                numexons = len(exonstarts) - 1
                numints = numexons - 1

                # Define number of exons-introns and introns-exons junctions:
                numeij = numexons - 1
                numiej = numints

                for i, exonstart in enumerate(exonstarts):
                    exonlength = abs(exonstart - exonends[i])
                    if exonlength > 0:
                        exnum = i if strand == '+' else numexons - i
                        exon_string = f"{chrom}\t{exonstart}\t{exonends[i]}\t{name}_exon_{exnum}_{score}_{chrom}_{exonstart+1}_{direction}\t{score}\t{strand}"
                        coding_list.append(exon_string)
                    if i < len(exonstarts)-1:
                        intronlength = abs(exonends[i] - exonstarts[i+1])
                        if intronlength > 0:
                            intnum = i if strand == '+' else numints - i
                            intronstart = exonends[i]
                            intronend = intronstart + intronlength
                            intron_string = f"{chrom}\t{intronstart}\t{intronend}\t{name}_intron_{intnum}_{score}_{chrom}_{intronstart+1}_{direction}\t{score}\t{strand}"
                            coding_list.append(intron_string)
                            
                            # Defines exon-intron and intron-exon junctions numbers
                            eijnum = i if strand == '+' else numeij - i
                            iejnum = i if strand == '+' else numiej - i

                            # Defines exon-intron and intron-exon junctions numbers coordinates
                            eijstart = exonends[i]-exside
                            eijend = exonends[i]+intside
                            iejstart = intronend-intside
                            iejend = intronend+exside
                            
                            # Defines exon-intron and intron-exon junctions strings
                            eij_string = f"{chrom}\t{eijstart}\t{eijend}\t{name}_eij_{eijnum}_{score}_{chrom}_{eijstart+1}_{direction}\t{score}\t{strand}"
                            junctions_list.append(eij_string)
                            iej_string = f"{chrom}\t{iejstart}\t{iejend}\t{name}_iej_{iejnum}_{score}_{chrom}_{iejstart+1}_{direction}\t{score}\t{strand}"
                            junctions_list.append(iej_string)

    # Write the output files
    with open(outputfile, 'w', encoding='utf-8') as outfile:
        for line in coding_list:
            outfile.write(f"{line}\n")

    with open(outputjfile, 'w', encoding='utf-8') as outjfile:
        for line in junctions_list:
            outjfile.write(f"{line}\n")

    return f"bed files saved in: {outputfile} and {outputjfile}"
