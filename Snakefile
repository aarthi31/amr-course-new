import json
from os.path import join, basename, dirname

# globals ----------------------------

configfile: 'config/config.yaml'
print(os.getcwd())
# Full path to a folder where intermediate output files will be created.
OUT_DIR = config['OUT_DIR']
#print(OUT_DIR)
FILES = json.load(open(config['SAMPLES_JSON']))
#print(FILES)
SAMPLES = sorted(FILES.keys())
print(SAMPLES)
FASTQC_DIR = config['FASTQC_DIR']
QC_DIR = config['QC_DIR']
DATA_DIR = os.getcwd() + '/'+ config['DATA_DIR']
#CONSENSUS_DIR = config['CONSENSUS_DIR']

#RESULTS_DIR = config['RESULTS_DIR']

FASTQ_CONVERTED = config['FASTQ_CONVERTED']

if not os.path.exists(FASTQ_CONVERTED):
    os.makedirs(FASTQ_CONVERTED)

rule debug:
    input:
        r1 = lambda wc: [config['SAMPLES_JSON'][sample][0] for sample in config['SAMPLES_JSON']]
        r1 = lambda wc: config['SAMPLES_JSON'][wc.sample][0],
        r2 = lambda wc: config['SAMPLES_JSON'][wc.sample][1],
        #r1=expand("data/{sample}_1.fastq.gz", sample=SAMPLES),
        #r2=expand("data/{sample}_R2.fastq.gz", sample=SAMPLES)
    output:
       # expand(FASTQ_CONVERTED + "/{sample}.txt", sample=SAMPLES)
       o1 = "lsres/{sample}.txt"
    shell:
        "ls -lh {input.r1}"

rule grep:
    input:
        lsres = rules.debug.output["o1"]
    output:
        grepres = "grepres/{sample}.txt"

#expand(FASTQ_CONVERTED + "/{sample}_R.fastq", sample=SAMPLES)
# rule all:
#     input:
#         # Debug raw data files
#         expand(DATA_DIR + "/{sample}_R1.fastq.gz", sample=SAMPLES),
#         #expand("data/{sample}_R2.fastq.gz", sample=SAMPLES),
# #     input:
# #         # expand("results/fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
# #         # expand("results/fastqc/{sample}_R2_fastqc.html", sample=SAMPLES),
# #         # expand("results/fastp/{sample}_filtered_R1.fastq.gz", sample=SAMPLES),
# #         # expand("results/fastp/{sample}_filtered_R2.fastq.gz", sample=SAMPLES),
# #         # Debug inputs
# #         expand(DATA_DIR+"/{sample}_R1.fastq.gz", sample=SAMPLES),
#         #expand("data/{sample}_R2.fastq.gz", sample=SAMPLES)
# rule debug:
#     input:
#         r1=expand(DATA_DIR+ "/{sample}_R1.fastq.gz", sample=SAMPLES)
#         #r2=expand(os.path.join("data", "{sample}_R2.fastq.gz"), sample=SAMPLES)
#     shell:
#         "echo {input.r1}"
# rule fastqc:
#     input:
#         r1=os.path.join("data","{sample}_R1.fastq.gz"),
#         r2=os.path.join("data","{sample}_R2.fastq.gz")
#     output:
#         r1_qc="results/fastqc/{sample}_R1_fastqc.html",
#         r2_qc="results/fastqc/{sample}_R2_fastqc.html",
#         r1_qc_zip="results/fastqc/{sample}_R1_fastqc.zip",
#         r2_qc_zip="results/fastqc/{sample}_R2_fastqc.zip"
#     params:
#         fastqc_opts="--threads 4"
#     shell:
#         "fastqc {params.fastqc_opts} -o results/fastqc {input.r1} {input.r2}"

# rule fastp:
#     input:
#         r1="data/{sample}_R1.fastq.gz",
#         r2="data/{sample}_R2.fastq.gz"
#     output:
#         r1_filtered="results/fastp/{sample}_filtered_R1.fastq.gz",
#         r2_filtered="results/fastp/{sample}_filtered_R2.fastq.gz",
#         json="results/fastp/{sample}_fastp.json",
#         html="results/fastp/{sample}_fastp.html"
#     params:
#         fastp_opts="--thread 4 --detect_adapter_for_pe"
#     shell:
#         "fastp {params.fastp_opts} --in1 {input.r1} --in2 {input.r2} \
#         --out1 {output.r1_filtered} --out2 {output.r2_filtered} \
#         --json {output.json} --html {output.html}"


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

