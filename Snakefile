import pandas as pd
from pathlib import Path

configfile: "config/config.yaml"

SAMPLES = pd.read_csv(config["samples_tsv"], sep="\t")
SAMPLE_NAMES = sorted(SAMPLES["sample"].unique().tolist())

def fastqs_for_sample(wc):
    return SAMPLES.loc[SAMPLES["sample"] == wc.sample, "fastq"].tolist()

MEDAKA_ON   = bool(config["params"]["medaka"]["enabled"])
PROKKA_ON   = bool(config["params"]["prokka"]["enabled"])
CHECKM_ON   = bool(config["params"]["checkm"]["enabled"])
IDGENES_ON  = bool(config["params"]["id_genes"]["enabled"])
AMR_ON      = bool(config["params"]["amrfinder"]["enabled"])

# Final targets
rule all:
    input:
        expand("results/{sample}/reads/{sample}.filtered.fastq.gz", sample=SAMPLE_NAMES),
        expand("results/{sample}/assembly/flye/assembly.fasta", sample=SAMPLE_NAMES),
        expand("results/{sample}/polish/medaka/polished_consensus.fasta", sample=SAMPLE_NAMES) if MEDAKA_ON else [],
        expand("results/{sample}/annotation/prokka/prokka_annotation.tsv", sample=SAMPLE_NAMES) if PROKKA_ON else [],
        expand("results/{sample}/qc/checkm/summary.tsv", sample=SAMPLE_NAMES) if CHECKM_ON else [],
        expand("results/{sample}/id_genes/16S_rRNA.fasta", sample=SAMPLE_NAMES) if IDGENES_ON else [],
        expand("results/{sample}/amr/amrfinder_pro_results.tsv", sample=SAMPLE_NAMES) if AMR_ON else [],

# Merge FASTQs per sample
rule merge_fastq:
    input:
        fastqs=fastqs_for_sample
    output:
        "results/{sample}/reads/{sample}.merged.fastq"
    threads: 1
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"

        # Portable across macOS + Linux:
        # - gzip -cd prints decompressed data to stdout
        # - cat for uncompressed files
        for f in {input.fastqs}; do
          case "$f" in
            *.gz) gzip -cd "$f" ;;
            *)    cat "$f" ;;
          esac
        done > "{output}"
        """

# Filter with chopper
rule chopper_filter:
    conda:
        "envs/chopper.yaml"
    input:
        "results/{sample}/reads/{sample}.merged.fastq"
    output:
        "results/{sample}/reads/{sample}.filtered.fastq.gz"
    threads: 4
    params:
        q=lambda wc: config["params"]["chopper"]["q"],
        minlength=lambda wc: config["params"]["chopper"]["minlength"]
    shell:
        r"""
        mkdir -p $(dirname {output})
        chopper -q {params.q} --minlength {params.minlength} --threads {threads} -i {input} \
          | pigz --fast -p {threads} > {output}
        """

# Flye assembly
rule flye_assembly:
    conda:
        "envs/flye.yaml"
    input:
        reads="results/{sample}/reads/{sample}.filtered.fastq.gz"
    output:
        fasta="results/{sample}/assembly/flye/assembly.fasta",
        info="results/{sample}/assembly/flye/assembly_info.txt",
        gfa="results/{sample}/assembly/flye/assembly_graph.gfa"
    threads: 24
    params:
        mode=lambda wc: config["params"]["flye"]["mode"]
    shell:
        r"""
        set -euo pipefail
        outdir=results/{wildcards.sample}/assembly/flye
        mkdir -p "$outdir"
        flye {params.mode} {input.reads} --out-dir "$outdir" --threads {threads}
        """

# Medaka polish (optional)
rule medaka_polish:
    conda:
        "envs/medaka.yaml"
    input:
        reads="results/{sample}/reads/{sample}.filtered.fastq.gz",
        draft="results/{sample}/assembly/flye/assembly.fasta"
    output:
        polished="results/{sample}/polish/medaka/polished_consensus.fasta"
    threads: 32
    params:
        model=lambda wc: config["params"]["medaka"]["model"]
    shell:
        r"""
        outdir=results/{wildcards.sample}/polish/medaka
        mkdir -p "$outdir"
        medaka_consensus -i {input.reads} -d {input.draft} -o "$outdir" -t {threads} -m {params.model}
        medaka sequence "$outdir"/*.hdf {input.draft} {output.polished}
        """

# Prokka annotate (optional; input is medaka if enabled else flye)
rule prokka:
    input:
        assembly=lambda wc: (
            f"results/{wc.sample}/polish/medaka/polished_consensus.fasta"
            if MEDAKA_ON else
            f"results/{wc.sample}/assembly/flye/assembly.fasta"
        )
    output:
        tsv="results/{sample}/annotation/prokka/prokka_annotation.tsv",
        faa="results/{sample}/annotation/prokka/prokka_annotation.faa",
        ffn="results/{sample}/annotation/prokka/prokka_annotation.ffn",
        gff="results/{sample}/annotation/prokka/prokka_annotation.gff"
    threads: 16
    conda: "envs/prokka.yaml"
    shell:
        r"""
        outdir=results/{wildcards.sample}/annotation/prokka
        mkdir -p "$outdir"
        prokka --outdir "$outdir" {input.assembly} --prefix prokka_annotation --cpus {threads} --force
        """

# CheckM (optional)
rule checkm:
    conda:
        "envs/checkm.yaml"
    input:
        asm=lambda wc: (
            f"results/{wc.sample}/polish/medaka/polished_consensus.fasta"
            if MEDAKA_ON else
            f"results/{wc.sample}/assembly/flye/assembly.fasta"
        )
    output:
        "results/{sample}/qc/checkm/summary.tsv"
    threads: 16
    shell:
        r"""
        outdir=results/{wildcards.sample}/qc/checkm
        workdir=$outdir/work
        mkdir -p "$workdir"

        # Stage the assembly as a .fasta in a directory (what your -x fasta workflow expects)
        cp {input.asm} "$workdir/{wildcards.sample}.fasta"

        # Run CheckM
        checkm lineage_wf --reduced_tree -t {threads} -x fasta "$workdir" "$outdir"

        # Write QA summary table to the declared output
        checkm qa -o 2 --tab_table -t {threads} -f {output} "$outdir/lineage.ms" "$outdir"
        """

# Identifier genes (optional)
rule identifier_genes:
    input:
        ffn="results/{sample}/annotation/prokka/prokka_annotation.ffn"
    output:
        rrna="results/{sample}/id_genes/16S_rRNA.fasta",
        dnaa="results/{sample}/id_genes/dnaA.fasta",
        rpob="results/{sample}/id_genes/rpoB.fasta"
    threads: 1
    shell:
        r"""
        outdir=results/{wildcards.sample}/id_genes
        mkdir -p "$outdir"
        echo '>16S ribosomal RNA' > {output.rrna}
        awk '/16S ribosomal RNA/{{flag=1;next}}/^>/{{flag=0}}flag' {input.ffn} >> {output.rrna}

        echo '>dnaA' > {output.dnaa}
        awk '/dnaA/{{flag=1;next}}/^>/{{flag=0}}flag' {input.ffn} >> {output.dnaa}

        echo '>rpoB' > {output.rpob}
        awk '/rpoB/{{flag=1;next}}/^>/{{flag=0}}flag' {input.ffn} >> {output.rpob}
        """

# AMRFinder database setup
rule amrfinder_db:
    output:
        touch("db/amrfinder_db_initialized.txt")
    conda:
        "envs/amrfinder.yaml"
    shell:
        r"""
        set -euo pipefail

        # Update the default internal DB location
        amrfinder --force_update

        # create marker file so Snakemake knows it's done
        mkdir -p db
        touch {output}
        """

# AMRFinder execution (optional)
rule amrfinder:
    input:
        faa="results/{sample}/annotation/prokka/prokka_annotation.faa",
        gff="results/{sample}/annotation/prokka/prokka_annotation.gff",
        db_init="db/amrfinder_db_initialized.txt"
    output:
        tsv="results/{sample}/amr/amrfinder_pro_results.tsv",
        proteins="results/{sample}/amr/amrfinder_pro.fasta"
    threads: 12
    conda: "envs/amrfinder.yaml"
    log:
        "logs/{sample}/amrfinder.log"
    shell:
        r"""
        set -euo pipefail
        outdir=results/{wildcards.sample}/amr
        mkdir -p "$outdir" "$(dirname "{log}")"

        amrfinder \
          -a prokka \
          --threads {threads} \
          -p {input.faa} \
          -g {input.gff} \
          --print_node --plus \
          --protein_output {output.proteins} \
          > {output.tsv} 2> "{log}"
        """