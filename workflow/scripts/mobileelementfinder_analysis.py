#!/usr/bin/env python

import argparse
from loguru import logger
import os
import subprocess
from Bio import SeqIO
import sys


def read_fa(input_fasta):
    fasta_sequences = SeqIO.parse(open(input_fasta), "fasta")

    # define a dictionary of chromosome keys and values
    description_to_id = {}
    for i in fasta_sequences:
        description_to_id[i.description] = i.id
    logger.debug(f"Fasta header to identifier is: {description_to_id}")
    return description_to_id


def bedformat_mobileelementfinder(mge_outputcsv, analysis_outputdir,
                                  description_to_id):
    mge_outputbed = os.path.join(analysis_outputdir, "mge_out.sorted.bed")
    output = subprocess.run(
        [
            f"grep -v '^#' {mge_outputcsv} \
            | csvtk cut -f contig,start,end,name,type -T \
            | sed 1d"
        ],
        capture_output=True,
        shell=True,
    )
    if output.returncode != 0:
        logger.error("Error in formatting mge output as bedfile!")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        with open(f"{mge_outputbed}.unsorted", "w") as f:
            for lines in output.stdout.decode().split("\n"):
                linelist = lines.split("\t")
                if len(linelist) == 0:
                    continue
                mge_id = linelist[0]
                if mge_id in description_to_id:
                    f.write(description_to_id[mge_id] + "\t" +
                            "\t".join(linelist[1:]) + "\n")
                else:
                    f.write("\t".join(linelist) + "\n")
        output = subprocess.run(
            [f"sort-bed {mge_outputbed}.unsorted > {mge_outputbed}"],
            check=True,
            shell=True,
        )
        logger.success(f"Completed reformatting mge output to {mge_outputbed}")
        return mge_outputbed


def classify_mobileelementfinder(inputbed, is_elements_bed, maxdist):
    bed = os.path.dirname(is_elements_bed)
    output_bed = os.path.join(bed, "input-mge_out-intersect.sorted.bed")

    logger.info(f"Classifying elements for TN/IS association with max distance {maxdist} bp.")

    # Ensure input BED files are sorted
    inputbed_sorted = f"{inputbed}.sorted"
    subprocess.run(
        [f"sort-bed {inputbed} > {inputbed_sorted}"],
        shell=True,
        check=True,
    )

    # Test if the element overlaps with an IS element using bedmap
    overlap_output = subprocess.run(
        [
            f"bedmap --echo --skip-unmapped {inputbed_sorted} {is_elements_bed}"
        ],
        capture_output=True,
        shell=True,
        text=True,
    )
    if overlap_output.returncode != 0:
        logger.error("Error in classifying elements overlapping IS elements! bedmap")
        logger.error(overlap_output.stdout)
        logger.error(overlap_output.stderr)
        sys.exit(2)
    else:
        overlap_lines = overlap_output.stdout.strip().split('\n')
        with open(f"{output_bed}.unsorted", "w") as f:
            for line in overlap_lines:
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    element_name = fields[3] + "|overlaps-IS"
                    f.write('\t'.join([fields[0], fields[1], fields[2], element_name]) + '\n')

    # Use closest-features to find the nearest features on both sides
    cmd_closest = f"closest-features --dist --delim '\t' --no-overlaps {inputbed_sorted} {is_elements_bed}"

    # Run closest-features
    output_closest = subprocess.run(
        [cmd_closest],
        capture_output=True,
        shell=True,
        text=True,
    )

    if output_closest.returncode != 0:
        logger.error("Error in finding closest IS elements! closest-features")
        logger.error(output_closest.stderr)
        sys.exit(2)
    else:
        lines = output_closest.stdout.strip().split('\n')
        flanked_elements = {}
        for line in lines:
            fields = line.strip().split('\t')
            if len(fields) < 16:
                continue  # Not enough fields, skip

            # Input element fields
            element_chrom = fields[0]
            element_start = fields[1]
            element_end = fields[2]
            element_name = fields[3]

            # Left (upstream) feature fields
            left_fields = fields[4:9]  # chrom, start, end, name, type
            left_distance = int(fields[9])

            # Right (downstream) feature fields
            right_fields = fields[10:15]  # chrom, start, end, name, type
            right_distance = int(fields[15])

            # Check if closest features were found
            if left_fields[0] == 'NA' or right_fields[0] == 'NA':
                continue  # Skip if no upstream or downstream IS element

            # Skip if distances exceed maxdist
            if abs(left_distance) > maxdist or abs(right_distance) > maxdist:
                continue

            # Extract IS types from the name fields
            left_is_name = left_fields[3]
            left_is_type = left_fields[4]
            right_is_name = right_fields[3]
            right_is_type = right_fields[4]

            if left_is_type == right_is_type:
                element_key = (element_chrom, element_start, element_end, element_name)
                element_name_flanked = element_name + f"|flanked-IS-{left_is_type}"
                flanked_elements[element_key] = element_name_flanked

        with open(f"{output_bed}.unsorted", "a") as f:
            for key, name in flanked_elements.items():
                f.write('\t'.join([key[0], key[1], key[2], name]) + '\n')

    # Sort the output bed file
    subprocess.run(
        [f"sort-bed {output_bed}.unsorted > {output_bed}"],
        check=True,
        shell=True,
    )
    logger.success("Completed classifying elements for TN/IS association")
    return output_bed




#def classify_mobileelementfinder(inputbed, bedmge, maxdist):
#    bed = os.path.dirname(bedmge)
#    output_bed = os.path.join(bed, "input-mge_out-intersect.sorted.bed")
#
#    # test if the element is nested inside anything in ME finder
#    output = subprocess.run(
#        [
#            f"bedmap --echo --echo-map-id-uniq --fraction-ref 1.0 {inputbed} {bedmge}\
#            | grep -v '|$' \
#            | awk -F'\\t' '{{gsub(/\\|/, \"|nested-\", $4); print}}'"
#        ],
#        capture_output=True,
#        shell=True,
#    )
#    if output.returncode != 0:
#        logger.error(
#            "Error in classifying nested elements in mge's results! bedmap")
#        logger.error(output.stdout.decode())
#        logger.error(output.stderr.decode())
#        sys.exit(2)
#    else:
#        with open(f"{output_bed}.unsorted", "w") as f:
#            f.write(str(output.stdout.decode()))
#
#    output_bed = os.path.join(bed, "input-mge_out-intersect.sorted.bed")
#
#    # test if the element is w/in max distance of two ME of the same type
#    output = subprocess.run(
#        [
#            f'bedmap --echo --echo-map-id --range {maxdist} {inputbed} {bedmge} \
#            | awk -F\'\\t\' \'{{split($4,a,"|"); split(a[2], b, ";"); for(i in b){{if(++count[b[i]] > 1){{$4=a[1]"|sandwiched-"b[i]; print; break}}}}; delete count}}\''
#        ],
#        capture_output=True,
#        shell=True,
#    )
#    if output.returncode != 0:
#        logger.error(
#            "Error in identifying sandwiched elements in mge's results! bedops"
#        )
#        logger.error(output.stdout.decode())
#        logger.error(output.stderr.decode())
#        sys.exit(2)
#    else:
#        with open(f"{output_bed}.unsorted", "a") as f:
#            f.write(str(output.stdout.decode()))
#    logger.debug(output.stderr.decode())
#    output = subprocess.run(
#        [f"sort-bed {output_bed}.unsorted > {output_bed}"],
#        check=True,
#        shell=True,
#    )
#    logger.success("Completed classifying mge output")
#    return output_bed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        "-i",
        help="[REQUIRED] Path to mobileelementfinder's output csv",
        type=os.path.abspath,
        required=True,
    )
    parser.add_argument(
        "--bed",
        "-b",
        help="[REQUIRED] Path to the input bed file to query for putative\
 mobility",
        type=os.path.abspath,
        required=True,
    )
    parser.add_argument(
        "--fasta",
        "-f",
        help="[REQUIRED] Path to the input fasta file",
        type=os.path.abspath,
        required=True,
    )
    parser.add_argument(
        "--maxdist",
        "-m",
        help=
        "Max. dist. to consider as mobile when sandwiched between two MGEs of the same type",
        type=int,
        default=10000,
        required=False,
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Directory for the output. Default is ./maybemobile_out/",
        type=os.path.abspath,
        default="./maybemobile_out/",
        required=False,
    )
    args = parser.parse_args()
    description_to_id = read_fa(args.fasta)
    # print(description_to_id)
    # print(args.input)
    bedmgeout = bedformat_mobileelementfinder(args.input, args.output,
                                              description_to_id)
    classify_mobileelementfinder(args.bed, bedmgeout, args.maxdist)


if __name__ == "__main__":
    main()
