

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
        cd qc/ && multiqc .
        """