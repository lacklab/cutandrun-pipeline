rule homer_annotatepeaks:
	input:
		"results/peaks/{raw}_peaks.narrowPeak"
	output:
		"results/homer/{raw}_annotatepeaks.txt"
	params:
		config["OUTPUT"]["REF"]
	shell:
		"""
		annotatePeaks.pl {input} {params} > {output}
		"""