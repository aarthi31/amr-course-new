import json
from os.path import join, basename, dirname

# globals ----------------------------

configfile: 'config/config.yaml'

# Full path to a folder where intermediate output files will be created.
OUT_DIR = config['OUT_DIR']
print(OUT_DIR)
FILES = json.load(open(config['SAMPLES_JSON']))
print(FILES)
SAMPLES = sorted(FILES.keys())
print(SAMPLES)
FASTQC_DIR = config['FASTQC_DIR']

QC_DIR = config['QC_DIR']

#CONSENSUS_DIR = config['CONSENSUS_DIR']

#RESULTS_DIR = config['RESULTS_DIR']

#FASTQ_CONVERTED = config['FASTQ_CONVERTED']

# if not os.path.exists(FASTQ_CONVERTED):
#     os.makedirs(FASTQ_CONVERTED)


# if not os.path.exists(FASTQC_DIR):
#     os.makedirs(FASTQC_DIR)

# if not os.path.exists(CONSENSUS_DIR):
#     os.makedirs(CONSENSUS_DIR)

# if not os.path.exists(RESULTS_DIR):
#     os.makedirs(RESULTS_DIR)

reads = ['1', '2']

rule all:
    input:
        expand("results/fastqc/{sample}_1_fastqc.html", sample=SAMPLES),
        expand("results/fastqc/{sample}_2_fastqc.html", sample=SAMPLES),
        expand("results/fastp/{sample}_filtered_1.fastq.gz", sample=SAMPLES),
        expand("results/fastp/{sample}_filtered_2.fastq.gz", sample=SAMPLES),
        expand("results/fastqcrerun/{sample}_filtered_1_fastqc.html", sample=SAMPLES),
        expand("results/fastqcrerun/{sample}_filtered_2_fastqc.html", sample=SAMPLES),
        expand("results/metaspades_output/{sample}/", sample=SAMPLES)
        # expand("results/metabat_output/{sample}/", sample=SAMPLES),
        # expand("results/checkm_output/{sample}/", sample=SAMPLES),
        # expand("results/kraken_folder/{sample}_kraken_report.txt", sample=SAMPLES),
        # expand("results/kraken_folder/{sample}_kraken_output.txt", sample=SAMPLES)
rule fastqc:
    input:
        r1="data/{sample}_1.fastq.gz",
        r2="data/{sample}_2.fastq.gz"
    output:
        r1_qc="results/fastqc/{sample}_1_fastqc.html",
        r2_qc="results/fastqc/{sample}_2_fastqc.html",
        r1_qc_zip="results/fastqc/{sample}_1_fastqc.zip",
        r2_qc_zip="results/fastqc/{sample}_2_fastqc.zip"
    params:
        fastqc_opts="--threads 4"
    shell:
        "fastqc {params.fastqc_opts} -o results/fastqc {input.r1} {input.r2}"

rule fastp:
    input:
        r1="data/{sample}_1.fastq.gz",
        r2="data/{sample}_2.fastq.gz"
    output:
        r1_filtered="results/fastp/{sample}_filtered_1.fastq.gz",
        r2_filtered="results/fastp/{sample}_filtered_2.fastq.gz",
        json="results/fastp/{sample}_fastp.json",
        html="results/fastp/{sample}_fastp.html"
    params:
        fastp_opts="--thread 4 --detect_adapter_for_pe"
    shell:
        "fastp {params.fastp_opts} --in1 {input.r1} --in2 {input.r2} \
        --out1 {output.r1_filtered} --out2 {output.r2_filtered} \
        --json {output.json} --html {output.html}"

rule fastqcrerun:
    input:
        r1="results/fastp/{sample}_filtered_1.fastq.gz",
        r2="results/fastp/{sample}_filtered_2.fastq.gz"
    output:
        r1_qc="results/fastqcrerun/{sample}_filtered_1_fastqc.html",
        r2_qc="results/fastqcrerun/{sample}_filtered_2_fastqc.html",
        r1_qc_zip="results/fastqcrerun/{sample}_filtered_1_fastqc.zip",
        r2_qc_zip="results/fastqcrerun/{sample}_filtered_2_fastqc.zip"
    params:
        fastqc_opts="--threads 4"
    shell:
        "fastqc {params.fastqc_opts} -o results/fastqcrerun {input.r1} {input.r2}"

rule assembly:
    input:
        r1="results/fastp/{sample}_filtered_1.fastq.gz",
        r2="results/fastp/{sample}_filtered_2.fastq.gz"
    output:
        assembly="results/metaspades_output/{sample}/"
    shell:
        "metaspades.py -1 {input.r1} -2 {input.r2} -o {output.assembly}"

rule binning:
    input:
        contig = "results/metaspades_output/{sample}/contigs.fasta"
    output:
        bins="results/metabat_output/{sample}/"
    shell:
        "metabat2 -i {input.contig} -o {output.bins} -m 1500"

# rule qualitycontrol:
#     input:
#         bins="results/metabat_output/{sample}/"
#     output:
#         cm = "results/checkm_output/{sample}/"
#     shell:
#         "checkm lineage_wf {input.bins} {output.cm}"

# rule kraken:
#     input:
#         contig = "results/metaspades_output/{sample}/contigs.fasta"
#     output:
#         krreport = "results/kraken_folder/{sample}_kraken_report.txt"
#         kroutput = "results/kraken_folder/{sample}_kraken_output.txt"
#     params:
#         db= config["kraken_db"]
#     shell:
#         "kraken2 --db {params.db} --threads 8 --output {output.kroutput} --report {output.krreport} {input.contig}"



# #include:"rules/convert_to_fastq.smk"
# #include:"rules/reverse_complement.smk"
# include:"rules/fastqc.smk"
# #include:"rules/adapter_lowqual.smk"
# #include:"rules/make_consensus.smk"
# #include:"rules/blast_rules.smk"
# #include:"rules/blast_top_hits.smk"
# #include:"rules/merge_all_res.smk"

# rule all:
#     input:
#        expand(FASTQ_CONVERTED + "/{sample}_1.fastq", sample=SAMPLES),
#        expand(FASTQ_CONVERTED + "/{sample}_2.fastq", sample=SAMPLES),
#        expand(FASTQ_CONVERTED + "/{sample}_rc.fastq", sample=SAMPLES),
#        expand(FASTQC_DIR + "/{sample}_{read}_fastqc.html", sample=SAMPLES, read=reads),
#        expand(FASTQC_DIR + "/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=reads),
#        #expand(QC_DIR + "/{sample}_{read}_QC.fastq", sample=SAMPLES, read=reads),
#        #expand(CONSENSUS_DIR + "/{sample}_consensus.fasta", sample=SAMPLES),
#        #expand(CONSENSUS_DIR + "/{sample}_consensus.aln", sample=SAMPLES),
#        #expand(OUT_DIR+"/blast/{sample}.tsv", sample=SAMPLES),
#        #expand(OUT_DIR+"/blast-tophits/{sample}-tophits.tsv", sample=SAMPLES),
#        #os.path.join(RESULTS_DIR, 'final_results.csv')

