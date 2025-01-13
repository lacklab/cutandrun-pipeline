rule align_bowtie2:
    input:
        fq1="trimmed/{raw}_1.trimmed.fastq.gz",
        fq2="trimmed/{raw}_2.trimmed.fastq.gz"
    output:
        bam="results/mapping/{raw}.raw.bam",
        log="qc/bowtie2/{raw}.bowtie2.log"
    params:
        idx=lambda wildcards: config["REFERENCES"][ref]["BOWTIE2_IDX"]
    threads: 32
    shell:
        """
        bowtie2 \
            -x {params.idx} \
            -1 {input.fq1} -2 {input.fq2} \
            --threads {threads} \
            --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 --minins 10 --maxins 700 --dovetail \
            2> {output.log} \
            | samtools view --threads {threads} -o {output.bam} -

        """


rule bam_process:
    input:
        "results/mapping/{raw}.raw.bam"
    output:
        bam="results/mapping/{raw}.target.sorted.bam"
    threads: 12
    params:
        fa=lambda wildcards: config["REFERENCES"][ref]["FA"]
    shell:
        """
        samtools \
            sort \
            -@ {threads} \
            -o {output.bam} \
            -T {wildcards.raw}.target.sorted \
            {input}

        samtools \
            index \
            -@ 1 \
            {output.bam}
        """




rule bam_filter:
    input:
        "results/mapping/{raw}.target.sorted.bam"
    output:
        bam="results/mapping/{raw}.target.filtered.sorted.bam",
        interbam=temp("results/mapping/{raw}.target.filtered.bam")
    threads: 12
    params:
        fa=lambda wildcards: config["REFERENCES"][ref]["FA"]
    shell:
        """
        samtools \
            view \
            --threads {threads} \
            -b -q 20 -F 0x004 -F 0x0008 -f 0x001 \
            {input} \
            > {output.interbam}

        samtools \
            sort \
            -@ {threads} \
            -o {output.bam} \
            -T {wildcards.raw}.target.filtered.sorted \
            {output.interbam}

        samtools \
            index \
            -@ 1 \
            {output.bam}
        """



rule bam_markdup:
    input:
        "results/mapping/{raw}.target.filtered.sorted.bam"
    output:
        bam="results/mapping/{raw}.target.markdup.sorted.bam",
        rg=temp("results/mapping/{raw}.target.filtered.rg.bam"),
        interbam=temp("results/mapping/{raw}.target.markdup.bam"),
        metric="qc/picard/{raw}.target.markdup.MarkDuplicates.metrics.txt"
    params:
        fa=lambda wildcards: config["REFERENCES"][ref]["FA"]
    shell:
        """
        picard AddOrReplaceReadGroups \
            I={input} \
            O={output.rg} \
            RGID=1 \
            RGLB=lib1 \
            RGPL=ILLUMINA \
            RGPU=unit1 \
            RGSM={wildcards.raw}

        picard \
            -Xmx29491M \
            MarkDuplicates \
            --ASSUME_SORT_ORDER coordinate --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT --TMP_DIR tmp \
            --INPUT {output.rg} \
            --OUTPUT {output.interbam} \
            --REFERENCE_SEQUENCE {params.fa} \
            --METRICS_FILE {output.metric}

        samtools \
            sort \
            -@ {threads} \
            -o {output.bam} \
            -T {wildcards.raw}.target.markdup.sorted \
            {output.interbam}

        samtools \
            index \
            -@ 1 \
            {output.bam}
        """


rule bam_dedup:
    input:
        "results/mapping/{raw}.target.markdup.sorted.bam"
    output:
        bam="results/mapping/{raw}.target.dedup.sorted.bam",
        interbam=temp("results/mapping/{raw}.target.dedup.bam"),
        metric="qc/picard/{raw}.target.dedup.MarkDuplicates.metrics.txt"
    threads: 12
    params:
        fa=lambda wildcards: config["REFERENCES"][ref]["FA"]
    shell:
        """
        picard \
            -Xmx29491M \
            MarkDuplicates \
            --ASSUME_SORT_ORDER coordinate --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY LENIENT --TMP_DIR tmp \
            --INPUT {input} \
            --OUTPUT {output.interbam} \
            --REFERENCE_SEQUENCE {params.fa} \
            --METRICS_FILE {output.metric}

        samtools \
            sort \
            -@ {threads} \
            -o {output.bam} \
            -T {wildcards.raw}.target.dedup.sorted \
            {output.interbam}

        samtools \
            index \
            -@ 1 \
            {output.bam}
        """



