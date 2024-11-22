intermediate_dir = config["intermediate_dir"]


rule bakta:
    output:
        f"{intermediate_dir}/{{batch}}/bakta/sequence_{{seqnum}}/{{sample}}.gff3",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
    params:
        db=config["bakta_db"],
    conda:
        "../envs/bakta.yml"
    shell:
        """
        bakta --db {params.db} --output $(dirname {output}) --force {input.fa}
        """


rule phispy:
    output:
        f"{intermediate_dir}/{{batch}}/PhiSpy/raw/sequence_{{seqnum}}/{{sample}}.phispy",
    input:
        gbk=f"{intermediate_dir}/{{batch}}/bakta/sequence_{{seqnum}}/{{sample}}.gff3",
    conda:
        "../envs/phispy.yml"
    shell:
        """
        PhiSpy.py {input} $(dirname {output})
        """


rule phigaro:
    output:
        f"{intermediate_dir}/{{batch}}/phigaro/raw/sequence_{{seqnum}}/{{sample}}.phigaro.tsv",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
    threads: config["threads"]
    params:
        pvogdb=config["pvog_db"],
    conda:
        "../envs/phigaro.yml"
    shell:
        """
        set -euo pipefail  # Exit on error, undefined variables, and errors in pipes

        echo "Processing {input.fa}"
        
        if [ -s {input.fa} ]; then
            echo "Input file is non-empty."
            
            if [ ! -d {params.pvogdb} ]; then
                echo "PVOG database directory not found. Setting up PVOG database at {params.pvogdb}."
                phigaro-setup --no-updatedb -p {params.pvogdb}
            else
                echo "PVOG database directory exists at {params.pvogdb}."
            fi
            
            sum_len=$(seqkit stats -Ta {input.fa} | csvtk cut -t -f "sum_len" | csvtk del-header)
            max_len=$(seqkit stats -Ta {input.fa} | csvtk cut -t -f "max_len" | csvtk del-header)
            echo "Total sum of sequence lengths: $sum_len"
            echo "Max sequence length: $max_len"
            
            if [ "$max_len" -gt 20000 ]; then
                echo "Found at least one sequence longer than 20kb."
                phigaro -f {input.fa} -t {threads} -o $(dirname {output}) -d -e tsv
            else
                echo "Max sequence length ($max_len) does not exceed 20kb. Creating empty output."
                touch {output}
            fi
        else
            echo "Input file is empty. Creating empty output."
            touch {output}
        fi
        """


rule integron_finder:
    output:
        f"{intermediate_dir}/{{batch}}/IntegronFinder/raw/sequence_{{seqnum}}/Results_Integron_Finder_{{sample}}/{{sample}}.summary",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
    conda:
        "../envs/integronfinder.yml"
    threads: config["threads"]
    shell:
        """
        OUTDIR=$(dirname $(dirname {output}))
        if [ -s {input.fa} ]; then
            integron_finder --local-max --cpu {threads} --outdir $OUTDIR {input.fa}
        else
            touch {output}
        fi
        """


rule plasmidfinder:
    output:
        f"{intermediate_dir}/{{batch}}/PlasmidFinder/raw/sequence_{{seqnum}}/{{sample}}/data.json",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
    params:
        mincov=config["plasmidfinder_mincov"],
        threshold=config["plasmidfinder_threshold"],
    conda:
        "../envs/plasmidfinder.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        if [ -s {input.fa} ]; then
            plasmidfinder.py -i {input.fa} -o $(dirname {output}) --mincov {params.mincov} --threshold {params.threshold}
        else
            touch {output}
        fi
        """


rule mob_recon:
    output:
        f"{intermediate_dir}/{{batch}}/mob-suite/raw/sequence_{{seqnum}}/{{sample}}/contig_report.txt",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
    conda:
        "../envs/mobrecon.yml"
    shell:
        """
        if [ -s {input.fa} ]; then
            mob_recon -i {input} -o $(dirname {output}) --force
        else
            touch {output}
        fi
        """


rule mobileelementfinder:
    output:
        f"{intermediate_dir}/{{batch}}/mobileelementfinder/raw/sequence_{{seqnum}}/{{sample}}/mge_results.csv",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
    threads: config["threads"]
    conda:
        "../envs/mobileelementfinder.yml"
    shell:
        """
        if [ -s {input.fa} ]; then
            mefinder find --contig {input} $(dirname {output})/mge_results -t {threads} --temp-dir $(dirname {output})/tmp-mge
        else
            touch {output}
        fi
        """
