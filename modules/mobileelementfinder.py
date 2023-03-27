import os
import subprocess
from loguru import logger
import sys


def run_mobileelementfinder(input_fasta, output_path, threads):
    mefinder_output = os.path.join(
        "/".join([output_path, "mobileElementFinder_out"])
    )
    if not os.path.exists(mefinder_output):
        os.makedirs(mefinder_output)
        logger.info(
            f"Created {mefinder_output} for output of mobileElementFinder"
        )
    mefinder_output = os.path.join("/".join([mefinder_output, "out"]))
    output = subprocess.run(
        [
            "mefinder",
            "find",
            "--contig",
            f"{input_fasta}",
            "--threads",
            f"{threads}",
            f"{mefinder_output}",
        ],
        capture_output=True,
    )
    if output.returncode != 0:
        logger.error("Error in mobileelementfinder!")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        logger.debug(output.stdout.decode())
        # logger.debug(output.stderr.decode())
        logger.success("Completed mobileelementfinder")
    return "".join([mefinder_output, ".csv"])


def bedformat_mobileelementfinder(mge_outputcsv, description_to_id):
    mge_outputbed = os.path.dirname(mge_outputcsv)
    mge_outputbed = os.path.join(
        "/".join([mge_outputbed, "mge_out.sorted.bed"])
    )
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
                    f.write(
                        description_to_id[mge_id]
                        + "\t"
                        + "\t".join(linelist[1:])
                        + "\n"
                    )
                else:
                    f.write("\t".join(linelist) + "\n")
        output = subprocess.run(
            [f"sort-bed {mge_outputbed}.unsorted > {mge_outputbed}"],
            check=True,
            shell=True,
        )
        logger.success(f"Completed reformatting mge output to {mge_outputbed}")
        return mge_outputbed


def classify_mobileelementfinder(inputbed, bedmge):
    maxdist = 20000

    bed = os.path.dirname(bedmge)
    output_bed = os.path.join(
        "/".join([bed, "input-mge_out-intersect.sorted.bed"])
    )

    # test if the element is inside anything in ME finder
    output = subprocess.run(
        [
            f"bedmap --echo --echo-map {inputbed} {bedmge} \
            | grep -v '|$' \
            | sed 's/|/\\t/' \
            | awk -F'\\t' '{{printf \"%s\\t%s\\t%s\\t%s|%s\\n\", $1, $2, $3, $7, $8}}'"
        ],
        # f"bedmap --echo --echo-map {inputbed} {bedmge} \
        # | grep -v '|$' \
        # | sed 's/|/\\t//' \
        # | cut -f 1,2,3,7-"
        capture_output=True,
        shell=True,
    )
    if output.returncode != 0:
        logger.error(
            "Error in classifying nested elements in mge's results! bedmap"
        )
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        with open(f"{output_bed}", "w") as f:
            f.write(str(output.stdout.decode()))

    # test if the element is w/in max distance of two ME of the same type
    output = subprocess.run(
        [
            f"bedops --range {maxdist} --element-of {bedmge} {inputbed} \
                | awk '{{print $4}}' \
                | sort \
                | uniq -d \
                | xargs -I{{}} awk '$4==\"{{}}\"' {bedmge} \
                | bedops --range {maxdist} -e 1 {inputbed} - \
                | awk 'FNR==NR{{a[$0];next}}$0 in a{{print $0\"\\t\"v}}' - v=\"$(bedops --range {maxdist} --element-of 1 {bedmge} {inputbed} | awk -F'\\t' '{{print $4,$5}}' | sort | uniq -d)\" {inputbed}"
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
        with open(f"{output_bed}", "a") as f:
            f.write(str(output.stdout.decode()))
    logger.debug(output.stderr.decode())
    logger.success("Completed classifying mge output")
    return output_bed

    # output = subprocess.run(
    # [
    # f"closest-features --dist --closest {inputbed} {bedmge} \
    # | awk -F'\\t' '{{split($7, a, \"|\");if(sqrt(a[2]^2) < {maxdist}){{print $0}}}}'"
    # ],
    # capture_output=True,
    # shell=True,
    # )
    # if output.returncode != 0:
    # logger.error("Error in classifying mge's results! closest-features")
    # logger.error(output.stdout.decode())
    # logger.error(output.stderr.decode())
    # sys.exit(2)
    # else:
    # with open(f"{output_bed}", "w") as f:
    # f.write(str(output.stdout.decode()))
