#!/usr/bin/env python

import argparse
import subprocess
import os
import sys
from loguru import logger
from Bio import SeqIO


def check_input(input_fasta, input_bed, output_bed):
    fasta_sequences = SeqIO.parse(open(input_fasta), "fasta")

    # define a dictionary of chromosome keys and values
    description_to_id = {}
    for i in fasta_sequences:
        description_to_id[i.description] = i.id
    logger.debug(f"Fasta header to identifier is: {description_to_id}")

    # exit = False
    # open the input BED file for reading
    with open(input_bed, "r") as input_file:
        # open the output BED file for writing
        with open(output_bed, "w") as output_file:
            # loop through each line in the input BED file
            for line in input_file:
                # split the line into columns
                fields = line.rstrip().split("\t")
                if len(fields) != 4:
                    logger.error(fields)
                    logger.error(
                        f"Number of columns is {len(fields)}! This is not 4! Please make it 4 columns!"
                    )
                    exit = True

                fields = [item.strip() for item in fields]
                # check if the chromosome in the line matches a key in the dictionary
                if fields[0] in description_to_id:
                    # exit = True
                    # logger.error(
                        # f"{fields[0]} found in bed file! Please use {description_to_id[fields[0]]} instead! Reformatted bedfile outputted here {output_bed}. Please rerun callmemobile with that bed file instead!"
                    # )
                    # if there is a match, replace the chromosome value with the dictionary value
                    fields[0] = description_to_id[fields[0]]

                if int(fields[1]) > int(fields[2]):
                    logger.warning(f"{fields[1]} is greater than {fields[2]}. Switching these")
                    start = fields[2]
                    end = fields[1]
                    fields[1] = start
                    fields[2] = end

                # write the updated line to the output BED file
                output_file.write("\t".join(fields) + "\n")
    return description_to_id


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        "-i",
        help="[REQUIRED] Path to input bed file",
        type=os.path.abspath,
        required=True,
    )
    parser.add_argument(
        "--fasta",
        "-f",
        help="[REQUIRED] Path to fasta file being analyzed",
        type=os.path.abspath,
        required=True,
    )
    parser.add_argument(
        "--output",
        "-o",
        help="[REQUIRED] Path to the output reformatted bed file",
        type=os.path.abspath,
        required=True,
    )
    args = parser.parse_args()

    description_to_id = check_input(args.fasta, args.input, args.output)
    logger.success("Completed reformatting")


if __name__ == "__main__":
    main()
