import subprocess
import os
import sys
from loguru import logger
from Bio import SeqIO


def run_plasmidfinder(input_fasta, output_path):
    plasmidfinder_output = os.path.join(output_path, "plasmidfinder_out")
    if not os.path.exists(plasmidfinder_output):
        os.makedirs(plasmidfinder_output)
        logger.info(
            f"Created {plasmidfinder_output} for output of plasmidfinder"
        )
    output = subprocess.run(
        [
            "plasmidfinder.py",
            "-i",
            f"{input_fasta}",
            "-o",
            f"{plasmidfinder_output}",
            "-x",
        ],
        capture_output=True,
    )
    if output.returncode != 0:
        logger.error("Error in plasmidfinder!")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        logger.debug(output.stdout.decode())
        # logger.debug(output.stderr.decode())
        logger.success("Completed plasmidfinder")
    return os.path.join(plasmidfinder_output, "results_tab.tsv")


def bedformat_plasmidfinder(pfinderout, description_to_id):
    pfinder_outputbed = os.path.dirname(pfinderout)
    pfinder_outputbed = os.path.join(
        pfinder_outputbed, "plasmidfinder_out.sorted.bed"
    )
    output = subprocess.run(
        [
            f"csvtk cut -t -f Contig,'Position in contig',Plasmid,Note {pfinderout} \
            | sed 1d"
        ],
        capture_output=True,
        shell=True,
    )
    if output.returncode != 0:
        logger.error("Error in formatting plasmidfinder output as bedfile!")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        with open(f"{pfinder_outputbed}.unsorted", "w") as f:
            for lines in output.stdout.decode().split("\n"):
                linelist = lines.split("\t")
                if len(linelist) <= 1:
                    continue
                pfinder_id = linelist[0]
                print(linelist[1])
                pfinder_coords = linelist[1].split("..")
                if pfinder_id in description_to_id:
                    f.write(
                        description_to_id[pfinder_id]
                        + "\t"
                        + "\t".join(pfinder_coords)
                        + "\t"
                        + "\t".join(linelist[2:])
                        + "\n"
                    )
                else:
                    f.write(
                        pfinder_id
                        + "\t"
                        + "\t".join(pfinder_coords)
                        + "\t"
                        + "\t".join(linelist[2:])
                        + "\n"
                    )
        output = subprocess.run(
            [f"sort-bed {pfinder_outputbed}.unsorted > {pfinder_outputbed}"],
            check=True,
            shell=True,
        )
        logger.success(
            f"Completed reformatting plasmidfinder output to {pfinder_outputbed}"
        )
        return pfinder_outputbed


def classify_plasmidfinder(inputfasta, inputbed, bedpfinder):
    max_plasmidlen = 200000
    maxdist = 20000

    bed = os.path.dirname(bedpfinder)
    output_bed = os.path.join(bed, "input-pfinder_out-intersect.sorted.bed")
    with open(f"{output_bed}", "w") as f:
        pass

    for contig in SeqIO.parse(inputfasta, "fasta"):
        output = subprocess.run(
            ["grep", contig.id, bedpfinder], stdout=subprocess.PIPE
        ).stdout.decode("utf-8")
        if not output:
            logger.info(
                f"No plasmid elements found on {contig.id}. Continuing"
            )
            continue
        output = subprocess.run(
            ["grep", contig.id, inputbed], stdout=subprocess.PIPE
        ).stdout.decode("utf-8")
        if not output:
            logger.info(f"No tested elements are on {contig.id}. Continuing")
            continue
        contig_len = len(contig)
        if contig_len > max_plasmidlen:
            logger.info(
                f"{contig.id} is {contig_len:,} bp which is greater than max \
plasmid len: {max_plasmidlen:,} bp."
            )
            logger.info(f"treating {contig.id} as a chromosome.")
            logger.info(
                f"searching for input elements within {maxdist} of plasmid elements"
            )
            output = subprocess.run(
                [
                    f"closest-features --chrom {contig.id} --dist --closest {inputbed} {bedpfinder} \
                    | grep -v '|NA|NA$' \
                    | awk -F'\\t' '{{split($7, a, \"|\");if(sqrt(a[2]^2) < {maxdist}){{print $0}}}}'"
                ],
                capture_output=True,
                shell=True,
            )
            if output.returncode != 0:
                logger.error(
                    f"Error in classifying PlasmidFinder's results! closest-features {contig.id}"
                )
                logger.error(output.stdout.decode())
                logger.error(output.stderr.decode())
                sys.exit(2)
            else:
                with open(f"{output_bed}", "w") as f:
                    f.write(str(output.stdout.decode()))
                # logger.debug(output.stderr.decode())
        else:
            logger.info(f"{contig.id} is {contig_len:,} bp")
            logger.info(
                f"treating {contig.id} as a plasmid. classifying all input elements on this contig as maybe mobile"
            )
            output = subprocess.run(
                [
                    # f"grep {contig.id} {bedpfinder} | cut -f4"
                    f'grep {contig.id} {inputbed} | sed "s~$~\\tplasmid-contig|$(printf \'%q\' "$(grep {contig.id} {bedpfinder} | cut -f4)")~"'
# | sed "s/$/\\tplasmid-contig|$(grep {contig.id} {bedpfinder} | cut -f4 | sed \'s/$/-/\')/"'
# | sed "s/$/\tplasmid-contig|$(grep {contig.id} {bedpfinder} | cut -f4 | tr "\\n" "-" | sed \'s/-$/{chr(92)}n/\')/"'
                ],
                capture_output=True,
                shell=True,
            )
            logger.debug(output)
            # sys.exit(2)
            if output.returncode != 0:
                logger.error(
                    f"Error in classifying PlasmidFinder's results! grep {contig.id}"
                )
                logger.error(output)
                logger.error(output.stdout.decode())
                logger.error(output.stderr.decode())
                sys.exit(2)
            else:
                with open(f"{output_bed}", "a") as f:
                    f.write(str(output.stdout.decode()))
    logger.success("Completed classifying PlasmidFinder output")
    return output_bed
