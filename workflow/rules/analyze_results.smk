intermediate_dir = config["intermediate_dir"]
output_dir = config["output_dir"]


rule plasmidfinder_bed:
    output:
        f"{intermediate_dir}/{{batch}}/PlasmidFinder/processed/sequence_{{seqnum}}/{{sample}}/plasmidfinder_out.sorted.bed",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
        bed=lambda wildcards: fn_cleanbed(wildcards.batch, wildcards.seqnum),
        plasmidfinderout=f"{intermediate_dir}/{{batch}}/PlasmidFinder/raw/sequence_{{seqnum}}/{{sample}}/data.json",
    params:
        script=Path(workflow.basedir) / "scripts/plasmidfinder_analysis.py",
    conda:
        "../envs/callmemobile.yml"
    shell:
        """
        {params.script} \\
                -i {input.plasmidfinderout} \\
                -b {input.bed} \\
                -f {input.fa} \\
                -o {output}
        """


rule mobsuite_bed:
    output:
        f"{intermediate_dir}/{{batch}}/mob-suite/processed/sequence_{{seqnum}}/{{sample}}/mobrecon_out/input-mobrecon_out-intersect.sorted.bed",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
        bed=lambda wildcards: fn_cleanbed(wildcards.batch, wildcards.seqnum),
        mobout=f"{intermediate_dir}/{{batch}}/mob-suite/raw/sequence_{{seqnum}}/{{sample}}/contig_report.txt",
    params:
        script=Path(workflow.basedir) / "scripts/mobsuite_analysis.py",
    conda:
        "../envs/callmemobile.yml"
    shell:
        """
        {params.script} \\
                -i $(dirname {input.mobout}) \\
                -b {input.bed} \\
                -f {input.fa} \\
                -o $(dirname {output})
        """


rule integronfinder_bed:
    output:
        f"{intermediate_dir}/{{batch}}/IntegronFinder/processed/sequence_{{seqnum}}/{{sample}}/input-ifinder_out-intersect.sorted.bed",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
        bed=lambda wildcards: fn_cleanbed(wildcards.batch, wildcards.seqnum),
        ifinderout=f"{intermediate_dir}/{{batch}}/IntegronFinder/raw/sequence_{{seqnum}}/Results_Integron_Finder_{{sample}}/{{sample}}.summary",
    params:
        script=Path(workflow.basedir) / "scripts/integronfinder_analysis.py",
        overlap=0.9,
    conda:
        "../envs/callmemobile.yml"
    shell:
        """
        {params.script} \\
                -i $(dirname {input.ifinderout}) \\
                -b {input.bed} \\
                -bo {params.overlap} \\
                -o $(dirname {output})
        """


rule mobileelementfinder_bed:
    output:
        f"{intermediate_dir}/{{batch}}/mobileelementfinder/processed/sequence_{{seqnum}}/{{sample}}/input-mge_out-intersect.sorted.bed",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
        bed=lambda wildcards: fn_cleanbed(wildcards.batch, wildcards.seqnum),
        mefinderout=f"{intermediate_dir}/{{batch}}/mobileelementfinder/raw/sequence_{{seqnum}}/{{sample}}/mge_results.csv",
    params:
        script=Path(workflow.basedir) / "scripts/mobileelementfinder_analysis.py",
    conda:
        "../envs/callmemobile.yml"
    shell:
        """
        {params.script} \\
                -i {input.mefinderout} \\
                -b {input.bed} \\
                -f {input.fa} \\
                -o $(dirname {output})
        """


rule phigaro_bed:
    output:
        f"{intermediate_dir}/{{batch}}/phigaro/processed/sequence_{{seqnum}}/{{sample}}/input-phigaro_out-intersect.sorted.bed",
    input:
        bed=lambda wildcards: fn_cleanbed(wildcards.batch, wildcards.seqnum),
        phigaro=f"{intermediate_dir}/{{batch}}/phigaro/raw/sequence_{{seqnum}}/{{sample}}.phigaro.tsv",
    params:
        script=Path(workflow.basedir) / "scripts/phigaro_analysis.py",
        maxdist=10000,
    conda:
        "../envs/callmemobile.yml"
    shell:
        """
        {params.script} \\
                -i {input.phigaro} \\
                -b {input.bed} \\
                -o $(dirname {output}) \\
                -m {params.maxdist}
        """

rule aggregate_output:
    output:
        f"{output_dir}/{{batch}}/sequence_{{seqnum}}-{{sample}}-callmemobile.tsv",
    input:
        bed=lambda wildcards: fn_cleanbed(wildcards.batch, wildcards.seqnum),
        phigaro=lambda wildcards: fn_phigaroprocessed(wildcards.batch, wildcards.seqnum, wildcards.sample), 
        mobrecon=lambda wildcards: fn_mobprocessed(wildcards.batch, wildcards.seqnum, wildcards.sample), 
        plasmidfinder=lambda wildcards: fn_plasmidfinderprocessed(wildcards.batch, wildcards.seqnum, wildcards.sample), 
        mefinder=lambda wildcards: fn_mefinderprocessed(wildcards.batch, wildcards.seqnum, wildcards.sample), 
        integronfinder=lambda wildcards: fn_integronfinderprocessed(wildcards.batch, wildcards.seqnum, wildcards.sample), 
    params:
        script=Path(workflow.basedir) / "scripts/callmemobile.py",
    conda:
        "../envs/callmemobile.yml"
    shell:
        """
        {params.script} \\
                --integronfinder {input.integronfinder} \\
                --plasmidfinder {input.plasmidfinder} \\
                --mob_suite {input.mobrecon} \\
                --phigaro {input.phigaro} \\
                --mobileelementfinder {input.mefinder} \\
                --bed {input.bed} \\
                -o {output} 
        """
