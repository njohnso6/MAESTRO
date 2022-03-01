rule scatac_countpeak:
    input:
        finalpeak = "Result/ATAC/Analysis/%s_final_peaks.bed" %(config["outprefix"]),
        validbarcode = "Result/ATAC/QC/%s_scATAC_validcells.txt" %(config["outprefix"]),
        frag = "Result/ATAC/Mapping/fragments_corrected_dedup_count.tsv"
    output:
        counts = "Result/ATAC/Analysis/" + config["outprefix"] + "_peak_count.h5"
    params:
        peakcount = temp("Result/ATAC/Analysis/%s_sorted_peakcount.tsv" %(config["outprefix"])),
        species = config["species"],
        outdir = "Result/ATAC/Analysis",
        outpre = config["outprefix"]
    threads:
        config["cores"]
    benchmark:
        "Result/Benchmark/%s_scATAC_PeakCount.benchmark" %(config["outprefix"])
    shell:
        """
        LC_ALL=C fgrep -f {input.validbarcode} -w {input.frag} | bedtools intersect -sorted -wa -wb -a {input.finalpeak} -b - | sort -k1,1 -k2,2 -k3,3 -k7,7 - | bedtools groupby -g 1,2,3,7 -c 7 -o count > {params.peakcount};
        MAESTRO scatac-peakcount --binary --filtered_count {params.peakcount} --species {params.species} --directory {params.outdir} --outprefix {params.outpre}
        """

rule scatac_qcfilter:
    input:
        counts = "Result/ATAC/Analysis/%s_peak_count.h5" %(config["outprefix"]),
    output:
        filtercount = "Result/ATAC/QC/%s_filtered_peak_count.h5" %(config["outprefix"]),
    params:
        outdir = "Result/ATAC/QC",
        outpre = config["outprefix"],
        peak = config["cutoff"]["atac_peak"],
        cell = config["cutoff"]["atac_cell"],
    benchmark:
        "Result/Benchmark/%s_scATAC_QCFilter.benchmark" %(config["outprefix"])
    shell:
        "MAESTRO scatac-qc --format h5 --peakcount {input.counts} --peak-cutoff {params.peak} --cell-cutoff {params.cell} "
        "--directory {params.outdir} --outprefix {params.outpre}"
