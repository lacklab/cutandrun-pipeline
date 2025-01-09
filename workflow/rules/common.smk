
import pandas as pd
import urllib.request
import os
import pysam

samples = pd.read_table(config["SAMPLES"])
samples["Raw"] = samples["Name"] + "." + samples["Unit"].astype(str)



def get_fqs(wildcards):

	fq1 = samples.loc[samples["Raw"] == wildcards.raw, "Fastq1"].unique()[0]
	fq2 = samples.loc[samples["Raw"] == wildcards.raw, "Fastq2"].unique()[0]
	
	return fq1, fq2


##########


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

if config["OUTPUT"]["RUN"]["FINGERPRINTS"]:
	outputs += [
		f"qc/deeptools/{name}.plotFingerprint.pdf"
		for name in samples['Name'].unique()
	]

if config["OUTPUT"]["RUN"]["BAMSUMMARY"]:
	outputs += [
		f"qc/deeptools/all_bam.bamSummary.npz"
		for name in samples['Name'].unique()
	]

# <<< OUTPUTS <<<