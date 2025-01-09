

rule link:
    input:
        get_fqs
    output:
        link_fq1="link/{raw}_1.fastq.gz",
        link_fq2="link/{raw}_2.fastq.gz"
    threads: 8
    shell:
        """
        ln -s {input[0]}  {output.link_fq1}
        ln -s {input[1]}  {output.link_fq2}
        """



rule trim_adapters:
    input:
        link_fq1="link/{raw}_1.fastq.gz",
        link_fq2="link/{raw}_2.fastq.gz"
    output:
        trimmed_fq1="trimmed/{raw}_1.trimmed.fastq.gz",
        trimmed_fq2="trimmed/{raw}_2.trimmed.fastq.gz",
        fastqc="qc/fastqc/{raw}_1_fastqc.html",
        tfastqc="qc/trimgalore/{raw}_1.trimmed_fastqc.html"
    threads: 8
    shell:
        """
        fastqc \
            --quiet \
            --threads {threads} \
            {input.link_fq1} {input.link_fq2} \
            -o qc/fastqc/

        cd trimmed

        trim_galore \
            --fastqc \
            --cores {threads} \
            --paired \
            --gzip \
            ../{input.link_fq1} ../{input.link_fq2}

        mv {wildcards.raw}_1_val_1.fq.gz {wildcards.raw}_1.trimmed.fastq.gz
        mv {wildcards.raw}_2_val_2.fq.gz {wildcards.raw}_2.trimmed.fastq.gz

        mv {wildcards.raw}_1_val_1_fastqc.html ../qc/trimgalore/{wildcards.raw}_1.trimmed_fastqc.html
        mv {wildcards.raw}_2_val_2_fastqc.html ../qc/trimgalore/{wildcards.raw}_2.trimmed_fastqc.html

        mv {wildcards.raw}_1_val_1_fastqc.zip ../qc/trimgalore/{wildcards.raw}_1.trimmed_fastqc.zip
        mv {wildcards.raw}_2_val_2_fastqc.zip ../qc/trimgalore/{wildcards.raw}_2.trimmed_fastqc.zip
        
        mv {wildcards.raw}_1.fastq.gz_trimming_report.txt ../qc/trimgalore/{wildcards.raw}_1.fastq.gz_trimming_report.txt
        mv {wildcards.raw}_2.fastq.gz_trimming_report.txt ../qc/trimgalore/{wildcards.raw}_2.fastq.gz_trimming_report.txt
        
        """