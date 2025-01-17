import json
import os
from glob import glob

configfile: 'config/config.yaml'
FILES = json.load(open(config['SAMPLES_JSON']))
#print(FILES)
SAMPLES = sorted(FILES.keys())
#print(SAMPLES)

# Full path to a folder where input is and intermediate output files will be created.
INPUT_DATA_DIR = config['INPUT_DIR']
OUT_DIR = config['OUT_DIR']
reads = ['1', '2']
ABRICATE_OUT_DIR = f"{OUT_DIR}/abricate/"
BIN_DIR = f"{OUT_DIR}/metabat_output/"

rule all:
    input:
        expand(f"{OUT_DIR}/fastqc/{{sample}}_{{read}}_fastqc.html", sample=SAMPLES, read=reads),
        expand(f"{OUT_DIR}/fastqc/{{sample}}_{{read}}_fastqc.zip", sample=SAMPLES, read=reads),
        expand(f"{OUT_DIR}/fastp/{{sample}}_filtered_1.fastq.gz", sample=SAMPLES),
        expand(f"{OUT_DIR}/fastp/{{sample}}_filtered_2.fastq.gz", sample=SAMPLES),
        expand(f"{OUT_DIR}/fastqcrerun/{{sample}}_filtered_{{read}}_fastqc.html", sample=SAMPLES, read=reads),
        expand(f"{OUT_DIR}/fastqcrerun/{{sample}}_filtered_{{read}}_fastqc.zip", sample=SAMPLES, read=reads),
        expand(f"{OUT_DIR}/metaspades_output/{{sample}}/contigs.fasta", sample=SAMPLES),
        expand(f"{OUT_DIR}/metabat_output/{{sample}}/check.empty", sample=SAMPLES),
        expand(f"{OUT_DIR}/checkm_output/{{sample}}/lineage.ms", sample=SAMPLES),
        expand(f"{OUT_DIR}/kraken_folder/{{sample}}_kraken_report.txt", sample=SAMPLES),
        expand(f"{OUT_DIR}/kraken_folder/{{sample}}_kraken_output.txt", sample=SAMPLES),
        expand(f"{OUT_DIR}/kraken_reads/{{sample}}_kraken_report.txt", sample=SAMPLES),
        expand(f"{OUT_DIR}/kraken_reads/{{sample}}_kraken_output.txt", sample=SAMPLES),
        expand(f"{OUT_DIR}/prokka/{{sample}}_prokka.log", sample=SAMPLES),
        expand(f"{OUT_DIR}/abricate/{{sample}}/abricate_res.txt", sample=SAMPLES)

rule fastqc:
    input:
        r1=lambda wildcards: FILES[wildcards.sample]["1"][0],  # Access the 1st read for the sample
        r2=lambda wildcards: FILES[wildcards.sample]["2"][0] 
        #r1=f"{INPUT_DATA_DIR}/{{sample}}_1.fastq.gz",
        #r2=f"{INPUT_DATA_DIR}/{{sample}}_2.fastq.gz"
    output:
        r1_qc=f"{OUT_DIR}/fastqc/{{sample}}_{{reads}}_fastqc.html",
        r2_qc=f"{OUT_DIR}/fastqc/{{sample}}_{{reads}}_fastqc.zip"
    params:
        fastqc_opts="--threads 4"
    shell:
        "fastqc {params.fastqc_opts} {input.r1} {input.r2} --outdir=results/fastqc "

rule fastp:
    input:
        r1=lambda wildcards: FILES[wildcards.sample]["1"][0],  # Access the 1st read for the sample
        r2=lambda wildcards: FILES[wildcards.sample]["2"][0] 
        #r1=f"{INPUT_DATA_DIR}/{{sample}}_1.fastq.gz",
        #r2=f"{INPUT_DATA_DIR}/{{sample}}_2.fastq.gz"
    output:
        r1_filtered=f"{OUT_DIR}/fastp/{{sample}}_filtered_1.fastq.gz",
        r2_filtered=f"{OUT_DIR}/fastp/{{sample}}_filtered_2.fastq.gz",
        json=f"{OUT_DIR}/fastp/{{sample}}_fastp.json",
        html=f"{OUT_DIR}/fastp/{{sample}}_fastp.html"
    params:
        fastp_opts="--thread 4 --detect_adapter_for_pe"
    shell:
        "fastp {params.fastp_opts} --in1 {input.r1} --in2 {input.r2} \
        --out1 {output.r1_filtered} --out2 {output.r2_filtered} \
        --json {output.json} --html {output.html}"

rule fastqcrerun:
    input:
        r1=f"{OUT_DIR}/fastp/{{sample}}_filtered_1.fastq.gz",
        r2=f"{OUT_DIR}/fastp/{{sample}}_filtered_2.fastq.gz"
    output:
        r1_qc=f"{OUT_DIR}/fastqcrerun/{{sample}}_filtered_{{reads}}_fastqc.html",
        r2_qc=f"{OUT_DIR}/fastqcrerun/{{sample}}_filtered_{{reads}}_fastqc.zip"
    params:
        fastqc_opts="--threads 4",
        outdir=f"{OUT_DIR}/fastqcrerun"
    shell:
        "fastqc {params.fastqc_opts} --outdir={params.outdir} {input.r1} {input.r2}"

rule assembly:
    input:
        r1=f"{OUT_DIR}/fastp/{{sample}}_filtered_1.fastq.gz",
        r2=f"{OUT_DIR}/fastp/{{sample}}_filtered_2.fastq.gz"
    output:
        assembly=f"{OUT_DIR}/metaspades_output/{{sample}}/contigs.fasta"
    params:
        f"{OUT_DIR}/metaspades_output/{{sample}}"
    shell:
        "metaspades.py -1 {input.r1} -2 {input.r2} -o {params}"

rule binning:
    input:
        contig =f"{OUT_DIR}/metaspades_output/{{sample}}/contigs.fasta"
    output:
        f"{OUT_DIR}/metabat_output/{{sample}}/check.empty"
    params:
        f"{OUT_DIR}/metabat_output/{{sample}}/bin"
    shell:
        "metabat2 -i {input.contig} -o {params} -m 1500 && touch {output}"

rule qualitycheck:
    input:
        bins=f"{OUT_DIR}/metabat_output/{{sample}}/check.empty"
    output:
        cm = f"{OUT_DIR}/checkm_output/{{sample}}/lineage.ms"
    params:
        cmop = f"{OUT_DIR}/checkm_output/{{sample}}/",
        sample=f"{OUT_DIR}/metabat_output/{{sample}}/"
    shell:
        "checkm lineage_wf -x fa -t 8 {params.sample} {params.cmop} --pplacer_threads 8"

rule kraken:
     input:
         contig = f"{OUT_DIR}/metaspades_output/{{sample}}/contigs.fasta"
     output:
         krreport = f"{OUT_DIR}/kraken_folder/{{sample}}_kraken_report.txt",
         kroutput = f"{OUT_DIR}/kraken_folder/{{sample}}_kraken_output.txt"
     log:
        f"{OUT_DIR}/kraken_folder/{{sample}}.log"
     params:
         db= config["kraken_db"]
     shell:
         "kraken2 --db {params.db} --threads 8 --output {output.kroutput} --report {output.krreport} {input.contig} > {log} 2>&1"

rule kraken_reads:
     input:
        r1=f"{OUT_DIR}/fastp/{{sample}}_filtered_1.fastq.gz",
        r2=f"{OUT_DIR}/fastp/{{sample}}_filtered_2.fastq.gz"
     output:
         krreport = f"{OUT_DIR}/kraken_reads/{{sample}}_kraken_report.txt",
         kroutput = f"{OUT_DIR}/kraken_reads/{{sample}}_kraken_output.txt"
     log:
        f"{OUT_DIR}/kraken_reads/{{sample}}.log"
     params:
         db= config["kraken_db"]
     shell:
         "kraken2 --db {params.db} --threads 8 --output {output.kroutput} --report {output.krreport} {input.r1} {input.r2} > {log} 2>&1"

rule run_prokka:
     input:
        fasta = f"{OUT_DIR}/metaspades_output/{{sample}}/contigs.fasta"
     output:
        f"{OUT_DIR}/prokka/{{sample}}_prokka.log"  # Prokka generates a directory for each genome but need to specify file for dag
     params:
        outdir = f"{OUT_DIR}/prokka/{{sample}}",  # Output directory for Prokka
        prefix = "{sample}"  # Prefix for output files
     log:
        f"{OUT_DIR}/prokka/{{sample}}_prokka.log"
     shell:
        """
        #touch {output}
        prokka --outdir {params.outdir} --centre X --compliant --prefix {params.prefix} {input.fasta} > {log} 2>&1
        """
#print( expand(f"{OUT_DIR}/abricate/{{sample}}/{{bin}}/abricate_check.empty", sample=SAMPLES,
             #  bin=[bin for bin in get_all_bins()]))

rule run_abricate:
    input:
        f"{OUT_DIR}/metabat_output/{{sample}}/check.empty"
        #fasta=f"{OUT_DIR}/metabat_output/{{sample}}/{{bin}}.faa"
    output:
        f"{OUT_DIR}/abricate/{{sample}}/abricate_res.txt"
    params:
        #txt=f"{OUT_DIR}/abricate/{{sample}}/{{bin}}_abricate.txt",
        fasta=f"{OUT_DIR}/metabat_output/{{sample}}/"
    #log:
    #    f"{OUT_DIR}/abricate/{{sample}}/{{bin}}.log"
    shell:
        """
        abricate --db resfinder {params.fasta}/*.fa > {output}
        """
