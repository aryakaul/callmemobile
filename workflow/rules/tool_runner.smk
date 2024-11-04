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
    threads: 4
    params:
        pvogdb=config["pvog_db"],
    conda:
        "../envs/phigaro.yml"
    shell:
        """
        if [ ! -d {params.pvogdb} ]; then
            phigaro-setup --no-updatedb -p {params.pvogdb}
        fi
        phigaro -f {input} -t {threads} -o $(dirname {output}) -d -e tsv
        """


rule integron_finder:
    output:
        f"{intermediate_dir}/{{batch}}/IntegronFinder/raw/sequence_{{seqnum}}/Results_Integron_Finder_{{sample}}/{{sample}}.summary",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
    conda:
        "../envs/integronfinder.yml"
    shell:
        """
        OUTDIR=$(dirname $(dirname {output}))
        integron_finder --local-max --cpu 4 --outdir $OUTDIR {input.fa}
        """


rule plasmidfinder:
    output:
        f"{intermediate_dir}/{{batch}}/PlasmidFinder/raw/sequence_{{seqnum}}/{{sample}}/data.json",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
    conda:
        "../envs/plasmidfinder.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        plasmidfinder.py -i {input.fa} -o $(dirname {output}) --mincov 0.60 --threshold 0.90
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
        mob_recon -i {input} -o $(dirname {output}) --force
        """


rule mobileelementfinder:
    output:
        f"{intermediate_dir}/{{batch}}/mobileelementfinder/raw/sequence_{{seqnum}}/{{sample}}/mge_results.csv",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
    threads: 4
    conda:
        "../envs/mobileelementfinder.yml"
    shell:
        """
        mefinder find --contig {input} $(dirname {output})/mge_results -t {threads} --temp-dir $(dirname {output})/tmp-mge
        """
