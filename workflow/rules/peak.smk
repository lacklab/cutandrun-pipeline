


rule macs:
    input:
        get_macs_i
    output:
        narrowpeak="results/peaks/{raw}_peaks.narrowPeak",
        qc="qc/macs/{raw}_peaks.xls"
    threads:
        16
    params:
        get_macs_p
    shell:
        """
        macs3 callpeak \
            {params}  \
            -n results/peaks/{wildcards.raw} \
            -q 0.1 -g hs -f BAMPE

        mv results/peaks/{wildcards.raw}_peaks.xls {output.qc}
        """



rule genomecov:
    input:
        bam="results/mapping/{raw}.target.dedup.sorted.bam"
    output:
        bg=temp("results/bedgraph/{raw}.temp.bg"),
        clip="results/bedgraph/{raw}.clip.bg",
        bw="results/bigwig/{raw}.bw"
    params:
        chrSizes=config[f"REFERENCES"][ref]["CHROM_SIZES"]
    threads:
        16
    run:
        bam = pysam.AlignmentFile(input["bam"])
        scale = f" -scale {10**6 / bam.mapped}"
        
        shell("""
            bedtools genomecov \
                -ibam {input.bam} \
                -bg -pc {scale} \
                | sort -k1,1 -k2,2n --parallel={threads} > {output.bg}

            bedClip \
                {output.bg} \
                {params.chrSizes} \
                {output.clip}

            bedGraphToBigWig \
                {output.bg} \
                {params.chrSizes} \
                {output.bw}
            """)



rule seacr_stringent:
    input:
        get_seacr_i
    output:
        "results/peaks/{raw}.seacr.peaks.stringent.bed"
    params:
        mode=lambda wildcards: config["OUTPUT"]['SAECR_MODE']
    shell:
        """
        ~/apps/SEACR/SEACR_1.3.sh \
            {input} \
            non stringent \
            results/peaks/{wildcards.raw}.seacr.peaks
        """


rule seacr_relaxed:
    input:
        get_seacr_i
    output:
        "results/peaks/{raw}.seacr.peaks.relaxed.bed"
    params:
        mode=lambda wildcards: config["OUTPUT"]['SAECR_MODE']
    shell:
        """
        ~/apps/SEACR/SEACR_1.3.sh \
            {input} \
            non relaxed \
            results/peaks/{wildcards.raw}.seacr.peaks
        """









