rule fastqc:
    input:
       os.path.join('data', '{sample}_{read}.fastq.gz')
    output:
        os.path.join(FASTQC_DIR, '{sample}_{read}_fastqc.html'),
        os.path.join(FASTQC_DIR, '{sample}_{read}_fastqc.zip')
    # conda:
    #     "../envs/fastqc.yaml"
    shell:
        "fastqc {input} --outdir={FASTQC_DIR}"
