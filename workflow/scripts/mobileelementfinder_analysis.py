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


def classify_mobileelementfinder(inputbed, bedmge, maxdist):
    bed = os.path.dirname(bedmge)
    output_bed = os.path.join(bed, "input-mge_out-intersect.sorted.bed")

    # test if the element is nested inside anything in ME finder
    output = subprocess.run(
        [
            f"bedmap --echo --echo-map-id-uniq --fraction-ref 1.0 {inputbed} {bedmge}\
            | grep -v '|$' \
            | awk -F'\\t' '{{gsub(/\\|/, \"|nested-\", $4); print}}'"
        ],
        capture_output=True,
        shell=True,
    )
    if output.returncode != 0:
        logger.error(
            "Error in classifying nested elements in mge's results! bedmap")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        with open(f"{output_bed}.unsorted", "w") as f:
            f.write(str(output.stdout.decode()))

    output_bed = os.path.join(bed, "input-mge_out-intersect.sorted.bed")

    # test if the element is w/in max distance of two ME of the same type
    output = subprocess.run(
        [
            f'bedmap --echo --echo-map-id --range {maxdist} {inputbed} {bedmge} \
            | awk -F\'\\t\' \'{{split($4,a,"|"); split(a[2], b, ";"); for(i in b){{if(++count[b[i]] > 1){{$4=a[1]"|sandwiched-"b[i]; print; break}}}}; delete count}}\''
        ],
        capture_output=True,
        shell=True,
    )
    if output.returncode != 0:
        logger.error(
            "Error in identifying sandwiched elements in mge's results! bedops"
        )
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        with open(f"{output_bed}.unsorted", "a") as f:
            f.write(str(output.stdout.decode()))
    logger.debug(output.stderr.decode())
    output = subprocess.run(
        [f"sort-bed {output_bed}.unsorted > {output_bed}"],
        check=True,
        shell=True,
    )
    logger.success("Completed classifying mge output")
    return output_bed


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
        type=os.path.abspath,
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
