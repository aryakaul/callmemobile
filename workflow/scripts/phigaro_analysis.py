#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
from loguru import logger
import pandas as pd


def read_tsv(phigaro_tsv):
    """
    Reads a Phigaro TSV file and returns a DataFrame for further processing.
    """
    try:
        prophage_data = pd.read_csv(phigaro_tsv, sep='\t')
        logger.debug("Successfully read Phigaro TSV file.")
        return prophage_data
    except Exception as e:
        logger.error(f"Error reading Phigaro TSV file: {e}")
        sys.exit(1)


def format_phigaro_output(prophage_data, analysis_outputdir):
    """
    Formats the Phigaro TSV data into a BED-like file.
    """
    phigaro_bed = os.path.join(analysis_outputdir, "phigaro_out.bed")
    with open(phigaro_bed, "w") as bedfile:
        for _, row in prophage_data.iterrows():
            bedfile.write(
                f"{row['scaffold']}\t{row['begin']}\t{row['end']}\t{row['id']}\t{row['taxonomy']}\n"
            )
    logger.success(f"Phigaro TSV formatted to BED at {phigaro_bed}")
    return phigaro_bed


def classify_phigaro(inputbed, phigaro_bed, maxdist):
    bed = os.path.dirname(phigaro_bed)
    output_bed_unsorted = os.path.join(
        bed, "input-phigaro_out-intersect.bed.unsorted")
    output_bed = os.path.join(bed, "input-phigaro_out-intersect.sorted.bed")

    classified_elements = set()
    custom_delim = '\x1f'  # Using ASCII Unit Separator as a custom delimiter

    # First, find elements inside prophages
    cmd_inside = f"bedmap --echo --echo-map-id-uniq --fraction-ref 1.0 --delim '{custom_delim}' {inputbed} {phigaro_bed}"
    output = subprocess.run(cmd_inside, shell=True, capture_output=True)
    if output.returncode != 0:
        logger.error(
            "Error in classifying elements inside prophages using bedmap!")
        logger.error(output.stderr.decode())
        sys.exit(2)

    # Process the output for elements inside prophages
    with open(output_bed_unsorted, "w") as f:
        for line in output.stdout.decode().splitlines():
            if not line.strip():
                continue
            if line.endswith(custom_delim):
                # No overlap
                continue
            else:
                # Overlap exists
                parts = line.strip().split(custom_delim)
                bed_entry = parts[0]
                map_ids = parts[1] if len(parts) > 1 else ''
                fields = bed_entry.strip().split('\t')
                element_key = '\t'.join(fields[:3])  # chrom, start, end
                if element_key not in classified_elements:
                    classified_elements.add(element_key)
                    # Append "|inside-prophage" to the name field
                    if len(fields) >= 4:
                        fields[3] += "|inside-prophage"
                    else:
                        fields.append("inside-prophage")
                    f.write('\t'.join(fields) + '\n')

    # Now find elements within maxdist distance from prophages
    cmd_within = f"bedmap --echo --echo-map-id --range {maxdist} --delim '{custom_delim}' {inputbed} {phigaro_bed}"
    output = subprocess.run(cmd_within, shell=True, capture_output=True)
    if output.returncode != 0:
        logger.error(
            "Error in identifying elements near prophages using bedmap!")
        logger.error(output.stderr.decode())
        sys.exit(2)

    # Process the output for elements near prophages
    with open(output_bed_unsorted, "a") as f:
        for line in output.stdout.decode().splitlines():
            if not line.strip():
                continue
            if line.endswith(custom_delim):
                # No overlap
                continue
            else:
                # Overlap exists
                parts = line.strip().split(custom_delim)
                bed_entry = parts[0]
                map_ids = parts[1] if len(parts) > 1 else ''
                fields = bed_entry.strip().split('\t')
                element_key = '\t'.join(fields[:3])
                if element_key not in classified_elements:
                    classified_elements.add(element_key)
                    # Append "|near-prophage" to the name field
                    if len(fields) >= 4:
                        fields[3] += "|near-prophage"
                    else:
                        fields.append("near-prophage")
                    f.write('\t'.join(fields) + '\n')

    # Sort the output BED file
    cmd_sort = f"sort-bed {output_bed_unsorted} > {output_bed}"
    output = subprocess.run(cmd_sort, shell=True, capture_output=True)
    if output.returncode != 0:
        logger.error("Error in sorting the classified BED file!")
        logger.error(output.stderr.decode())
        sys.exit(2)

    logger.success("Completed classifying elements relative to Phigaro output")
    return output_bed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",
                        "-i",
                        help="[REQUIRED] Path to Phigaro TSV file",
                        type=os.path.abspath,
                        required=True)
    parser.add_argument("--bed",
                        "-b",
                        help="[REQUIRED] Path to the input BED file",
                        type=os.path.abspath,
                        required=True)
    parser.add_argument(
        "--output",
        "-o",
        help="Directory for output. Default is ./phigaro_analysis_out/",
        type=os.path.abspath,
        default="./phigaro_analysis_out/")
    parser.add_argument("--maxdist",
                        "-m",
                        help="Max distance to consider near a prophage",
                        type=int,
                        default=10000)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    prophage_data = read_tsv(args.input)
    prophage_bed = format_phigaro_output(prophage_data, args.output)
    classify_phigaro(args.bed, prophage_bed, args.maxdist)


if __name__ == "__main__":
    main()
