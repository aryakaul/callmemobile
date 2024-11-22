#!/usr/bin/env python

import argparse
from loguru import logger
import os
import subprocess
import sys


def format_integronfinderout(integronfinder_outdir, analysis_outdir):
    integron_outputs = os.path.join(integronfinder_outdir, "*.integrons")
    allintegrons = os.path.join(analysis_outdir,
                                "integronfinder_out.sorted.bed")
    with open(os.path.join(f"{allintegrons}"), 'w') as f:
        output = subprocess.run([
            f"sed 1d {integron_outputs} \
                | grep -v '^#' \
                | grep -v 'ID_integron' \
                | cut -f 2,4,5,11 \
                | awk '{{$2=($2<0)?0:$2; print}}' \
                | sort-bed -"
        ],
                                shell=True,
                                stdout=f)
        if output.returncode != 0:
            logger.error(
                "Error in formatting Integron Finder output as bedfile!")
            logger.error(output)
            logger.error(output.stdout.decode())
            logger.error(output.stderr.decode())
            sys.exit(2)
        else:
            logger.success(
                f"Completed reformatting Integron Finder output to {allintegrons}"
            )
            return allintegrons


def classify_integronfinder(input_bed, integron_bed, output_dir, bedolap):
    output_bed = os.path.join(output_dir,
                              "input-ifinder_out-intersect.sorted.bed")
    logger.info(
        f"Classifying IntegronFinder results. Need {bedolap} overlap of gene regions w/ input regions to classify."
    )
    output = subprocess.run(
        # [
        # f"bedmap --echo --echo-map-id-uniq --fraction-ref {bedolap} {integron_bed} {input_bed} \
        # | grep -v '|\s*$'"
        # ],
        [
            f"bedmap --echo --echo-map-id-uniq --fraction-map {bedolap} {input_bed} {integron_bed} \
            | grep -v '|\\s*$'"
        ],
        capture_output=True,
        shell=True,
    )
    logger.debug(output)
    with open(f"{output_bed}", "w") as f:
        f.write(str(output.stdout.decode()))
    logger.debug(output.stderr.decode())
    logger.success("Completed classifying Integron Finder output")
    return output_bed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        "-i",
        help="[REQUIRED] Path to IntegronFinder's output directory",
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
        "--overlap",
        "-bo",
        help="Percentage of input bed to overlap with Integron in order to be\
        classified as mobile",
        type=float,
        required=False,
        default=0.95)
    parser.add_argument(
        "--output",
        "-o",
        help="Directory for the output. Default is ./maybemobile_out/",
        type=os.path.abspath,
        default="./maybemobile_out/",
        required=False,
    )
    args = parser.parse_args()
    integron_bedformat = format_integronfinderout(args.input, args.output)
    classify_integronfinder(args.bed, integron_bedformat, args.output,
                            args.overlap)


if __name__ == "__main__":
    main()
