#!/usr/bin/env python3

from loguru import logger
import os
import subprocess
import sys


def run_integronfinder(input_fasta, output_path, threads):
    integronfinder_output = os.path.join(
        "/".join([output_path, "integronfinder_out"])
    )
    if not os.path.exists(integronfinder_output):
        os.makedirs(integronfinder_output)
        logger.info(
            f"Created {integronfinder_output} for output of Integron-Finder"
        )
    input_fasta_basename = os.path.basename(input_fasta)
    input_fasta_basename = os.path.splitext(input_fasta_basename)[0]
    output = subprocess.run(
        [
            "integron_finder",
            "--local-max",
            "--cpu",
            f"{threads}",
            "--circ",
            "--outdir",
            f"{integronfinder_output}",
            f"{input_fasta}",
        ],
        capture_output=True,
    )
    if output.returncode != 0:
        logger.error("Error in IntegronFinder!")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        logger.debug(output.stdout.decode())
        # logger.debug(output.stderr.decode())
        logger.success("Completed IntegronFinder")
    integronfinder_output_integrons = os.path.join(
        "/".join(
            [
                integronfinder_output,
                "Results_Integron_Finder_" + input_fasta_basename,
                input_fasta_basename + ".integrons",
            ]
        )
    )

    return integronfinder_output_integrons


def bedformat_integronfinder(ifinder_out):
    ifinder_outputbed = os.path.dirname(ifinder_out)
    ifinder_outputbed = os.path.join(
        "/".join([ifinder_outputbed, "integronfinder_out.sorted.bed"])
    )
    output = subprocess.run(
        [
            f"grep -v '^#' {ifinder_out} \
            | sed 1d \
            | cut -f2,4,5,11 \
            | sort-bed -",
        ],
        capture_output=True,
        shell=True,
    )
    if output.returncode != 0:
        logger.error("Error in formatting Integron Finder output as bedfile!")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        with open(f"{ifinder_outputbed}", "w") as f:
            f.write(str(output.stdout.decode()))
        # logger.debug(output.stderr.decode())
        logger.success(
            f"Completed reformatting Integron Finder output to {ifinder_outputbed}"
        )
        return ifinder_outputbed


def classify_integronfinder(input_bed, bedifinder):
    bed = os.path.dirname(bedifinder)
    bedolap = "0.9"
    output_bed = os.path.join(bed, "input-ifinder_out-intersect.sorted.bed")
    logger.info(
        f"Classifying IntegronFinder results. Need {bedolap} overlap of gene regions w/ input regions to classify."
    )
    output = subprocess.run(
        [
            f"bedmap --echo --echo-map-id-uniq --fraction-ref {bedolap} {input_bed} {bedifinder} \
            | grep -v '|$'"
        ],
        capture_output=True,
        shell=True,
    )
    if output.returncode != 0:
        logger.error(
            "Error in classifying Integron Finder's results! bedmap"
        )
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        with open(f"{output_bed}", "w") as f:
            f.write(str(output.stdout.decode()))
        logger.debug(output.stderr.decode())
    logger.success("Completed classifying Integron Finder output")
    return output_bed
