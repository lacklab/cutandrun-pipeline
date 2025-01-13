
import pandas as pd
import urllib.request
import os
import pysam

samples = pd.read_table(config["SAMPLES"])
samples["Raw"] = samples["Name"] + "." + samples["Unit"].astype(str)


assets = {}
with open("assets/annotatepeaks.asset", "r") as f:
	assets["annotatepeaks"]	= ""
	for line in f.readlines():
		assets["annotatepeaks"]	+= line



##########

def get_fqs(wildcards):

	fq1 = samples.loc[samples["Raw"] == wildcards.raw, "Fastq1"].unique()[0]
	fq2 = samples.loc[samples["Raw"] == wildcards.raw, "Fastq2"].unique()[0]
	
	return fq1, fq2


def get_control(wildcards):
	return samples.loc[samples["Raw"] == wildcards.raw, "Control"].unique()[0]

# Peak-calling functions
def get_macs_p(wildcards):
    control = get_control(wildcards)
    param = f"-t results/mapping/{wildcards.raw}.target.dedup.sorted.bam "
    if control != "-":
        param += f"-c results/mapping/{control}.target.dedup.sorted.bam "
    return param

def get_macs_i(wildcards):
    control = get_control(wildcards)
    inputs = [f"results/mapping/{wildcards.raw}.target.dedup.sorted.bam"]
    if control != "-":
        inputs.append(f"results/mapping/{control}.target.dedup.sorted.bam")
    return inputs



def get_seacr_i(wildcards):
    control = get_control(wildcards)
    inputs = [f"results/bedgraph/{wildcards.raw}.clip.bg"]
    if control != "-":
        inputs.append(f"results/bedgraph/{control}.clip.bg")
    return inputs



def get_group_fingerprint(wildcards):
	if wildcards.group == 'all_bam':
		return expand(
			"results/mapping/{raw}.target.markdup.sorted.bam", 
			raw=samples["Raw"].tolist()
		)
	else:
		return expand(
			"results/mapping/{raw}.target.markdup.sorted.bam", 
			raw=samples.loc[samples['Name'] == wildcards.group, "Raw"].tolist()
		)


def get_group_multibam(wildcards):
	if wildcards.group == 'all_bam':
		return expand(
			"results/mapping/{raw}.target.dedup.sorted.bam", 
			raw=samples["Raw"].tolist()
		)
	else:

		return expand(
			"results/mapping/{raw}.target.dedup.sorted.bam", 
			raw=samples.loc[samples['Name'] == wildcards.group, "Raw"].tolist()
		)

def get_multiqc(wildcards):
    out = []
    
    # List of common quality control file types and tools
    qc_tools = {
        "fastqc": [
            "{raw}_1_fastqc.html", 
            "{raw}_2_fastqc.html"
        ],
        "trimgalore": [
            "{raw}_1.fastq.gz_trimming_report.txt",
            "{raw}_2.fastq.gz_trimming_report.txt",
            "{raw}_1.trimmed_fastqc.html",
            "{raw}_2.trimmed_fastqc.html"
        ],
        "samtools": [
            "flagstat/{raw}.target.sorted.flagstat",
            "flagstat/{raw}.target.filtered.sorted.flagstat",
            "flagstat/{raw}.target.markdup.sorted.flagstat",
            "flagstat/{raw}.target.dedup.sorted.flagstat",
            "idxstats/{raw}.target.sorted.idxstats",
            "idxstats/{raw}.target.filtered.sorted.idxstats",
            "idxstats/{raw}.target.markdup.sorted.idxstats",
            "idxstats/{raw}.target.dedup.sorted.idxstats",
            "stats/{raw}.target.sorted.stats",
            "stats/{raw}.target.filtered.sorted.stats",
            "stats/{raw}.target.markdup.sorted.stats",
            "stats/{raw}.target.dedup.sorted.stats"
        ],
        "picard": [
            "{raw}.target.markdup.MarkDuplicates.metrics.txt",
            "{raw}.target.dedup.MarkDuplicates.metrics.txt"
        ],
        "bowtie2": [
            "{raw}.bowtie2.log"
        ],
        "macs": [
            "{raw}_peaks.xls"
        ],
        "deeptools": [
            "all_bam.bamSummary.npz",
            "all_bam.plotCorrelation.mat.tab",
            "all_bam.plotFingerprint.qcmetrics.txt",
            "all_bam.plotFingerprint.raw.txt",
            "all_bam.plotPCA.tab"
        ],
        "homer":
        [
            "{raw}_annotatepeaks.txt"
        ]
    }

    # Iterate through each sample and append all files based on the defined templates
    for _, row in samples.iterrows():
        raw = row['Raw']

        is_control = row['Control'] == '-'
        
        # Generate output paths for each tool and file pattern
        for tool, patterns in qc_tools.items():
            for pattern in patterns:
                if (tool in ['macs', 'homer']) and (not is_control):
                    continue
                out.append(f"qc/{tool}/{pattern.format(raw=raw)}")

    # Add FRIP score file (outside the loop as a single file)
    out.append("qc/frip_mqc.tsv")
    
    return expand(out)



##########

ref = config["OUTPUT"]["REF"]


outputs = []

if config["OUTPUT"]["RUN"]["QC"]:
	outputs += ["qc/multiqc_report.html"]


if config["OUTPUT"]["RUN"]["PEAKS"]:
	outputs += [
		f"results/peaks/{row['Raw']}.seacr.peaks.{mode}.bed"
		for i, row in samples.iterrows()
		for mode in config['OUTPUT']['SAECR_MODE']
		if row['Control'] != '-'
	]

	outputs += [
		f"results/peaks/{row['Raw']}_peaks.narrowPeak"
		for i, row in samples.iterrows()
		if row['Control'] != '-'
	]

##########

if config["OUTPUT"]["RUN"]["BWS"]:
	outputs += [
		f"results/bigwig/{row['Raw']}.bw"
		for i, row in samples.iterrows()
	]


# <<< OUTPUTS <<<