


rule bam_qc:
    input:
        "results/mapping/{bam}.bam"
    output:
        idxstats="qc/samtools/idxstats/{bam}.idxstats",
        flagstat="qc/samtools/flagstat/{bam}.flagstat",
        stats="qc/samtools/stats/{bam}.stats"
    params:
        fa=lambda wildcards: config["REFERENCES"][ref]["FA"]
    shell:
        """
        samtools \
            idxstats \
            --threads 0 \
            {input} \
            > {output.idxstats}

        samtools \
            flagstat \
            --threads 1 \
            {input} \
            > {output.flagstat}

        samtools \
            stats \
            --threads 1 \
            --reference {params.fa} \
            {input} \
            > {output.stats}
        """



################################################################################


from collections import Counter
rule annotatepeaks_qc:
	input:
		"results/homer/{raw}_annotatepeaks.txt"
	output:
		"qc/homer/{raw}_summary_mqc.txt"
	run:
		header = ["INTERGENIC", "INTRON ", "PROMOTER-TSS ", "EXON ", "3' UTR ", "5' UTR ", "TTS ", "NON-CODING "]
		with open(output[0], "w") as f:
			f.write(assets["annotatepeaks"])
			tmp = pd.read_table(input[0])
			if tmp.shape[0] == 0:
				nAnnot = dict(zip(header, [8]*0))
			else:
				tmp["shortAnn"] = tmp["Annotation"].str.split("(", expand=True)[0].str.upper()
				nAnnot = Counter(tmp["shortAnn"])
			for k in header:
				f.write(f"{k}\t{nAnnot[k]}\n")





rule frip:
    input:
        bams=expand("results/mapping/{raw}.target.dedup.sorted.bam", raw=samples["Raw"].tolist()),   # Fetch BAM files for FRIP calculation
        peak=expand("results/peaks/{raw}_peaks.narrowPeak", raw=samples["Raw"].tolist())    # Fetch peak files for FRIP calculation
    output:
        "qc/frip_mqc.tsv"
    run:
        import deeptools.countReadsPerBin as crpb
        import pysam
        import numpy as np

        # Prepare the MultiQC-compatible header
        with open(output[0], "w") as f:
            f.write("# plot_type: 'generalstats'\n")
            f.write("Sample Name\tFRiP\tNumber of Peaks\tMedian Fragment Length\n")

            # Loop through BAM and Peak pairs
            for b, p in zip(input.bams, input.peak):
                
                # Calculate FRIP using deepTools
                cr = crpb.CountReadsPerBin([b], bedFile=[p], numberOfProcessors=10)
                reads_at_peaks = cr.run()
                total_reads_at_peaks = reads_at_peaks.sum(axis=0)

                # Calculate total mapped reads using pysam
                bam = pysam.AlignmentFile(b)
                total_mapped_reads = bam.mapped

                # Calculate number of peaks
                with open(p, 'r') as peak_file:
                    num_peaks = sum(1 for _ in peak_file)

                # Calculate median fragment length using pysam
                fragment_lengths = [
                    abs(read.template_length) for read in bam.fetch() if read.is_proper_pair
                ]
                median_fragment_length = np.median(fragment_lengths)

                # Calculate FRIP score
                sample_name = p.split("/")[-1].split("_peaks")[0]
                frip_score = float(total_reads_at_peaks[0]) / total_mapped_reads

                # Write results into a MultiQC-compatible TSV file
                f.write(f"{sample_name}\t{frip_score:.4f}\t{num_peaks}\t{median_fragment_length:.2f}\n")


rule multiqc:
    input:
        get_multiqc
    output:
        "qc/multiqc_report.html"
    shell:
        """
        cd qc/ && multiqc . --filename
        """