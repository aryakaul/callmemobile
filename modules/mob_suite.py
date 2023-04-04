import sys
import os
import subprocess
from loguru import logger
from Bio import SeqIO


def run_mobrecon(ifasta):
    mobrecon_output = os.path.join("/".join([output_path, "mobrecon_out"]))
    if not os.path.exists(mobrecon_output):
        os.makedirs(mobrecon_output)
        logger.info(f"Created {mobrecon_output} for output of mob_recon")
    output = subprocess.run(
        [
            "mob_recon",
            "--infile",
            f"{ifasta}",
            "--num_threads",
            f"{threads}",
            "--outdir",
            f"{mobrecon_output}",
            "--force",
        ],
        capture_output=True,
    )
    if output.returncode != 0:
        logger.error("Error in mob_recon!")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        logger.debug(output.stdout.decode())
        logger.success("Completed mob_recon")
    return os.path.join("/".join([mobrecon_output, "contig_report.txt"]))


def classify_mobrecon(input_fasta, input_bed, mobrecondir, outputdir):
    if not os.path.exists(os.path.join(outputdir, "mobrecon_out")):
        os.mkdir(os.path.join(outputdir, "mobrecon_out"))
    contigreport = os.path.join(mobrecondir, "contig_report.txt")
    mobtyper = os.path.join(mobrecondir, "mobtyper_results.txt")
    if not os.path.exists(mobtyper):
        logger.info(
            "No plasmids identified by mob_recon from the provided assembly sequence"
        )
        return ""
    else:
        output_bed = os.path.join(
            outputdir,
            "mobrecon_out",
            "input-mobrecon_out-intersect.sorted.bed",
        )
        with open(f"{output_bed}", "w") as f:
            pass
        for contig in SeqIO.parse(input_fasta, "fasta"):
            plasmidid = (
                subprocess.run(
                    [
                        f"grep \"{contig.description}\" {contigreport} \
                    | grep 'plasmid' \
                    | cut -f3"
                    ],
                    capture_output=True,
                    shell=True,
                )
                .stdout.decode()
                .strip()
            )
            if plasmidid:
                output = subprocess.run(
                    [
                        f'grep "{contig.id}" {input_bed} \
                        | sed "s/$/\\tplasmid-contig:{plasmidid}|"$(grep {plasmidid} {mobtyper} | cut -f14)"/"'
                    ],
                    capture_output=True,
                    shell=True,
                )
                print(output)
                if output.returncode != 0:
                    logger.error(
                        f"Error in classifying mobrecon's results! grep {contig.id} & {plasmidid}"
                    )
                    logger.error(output.stdout.decode())
                    logger.error(output.stderr.decode())
                    sys.exit(2)
                else:
                    with open(f"{output_bed}", "a") as f:
                        f.write(str(output.stdout.decode()))
    logger.success("Completed classifying mobrecon's output")
    return output_bed
