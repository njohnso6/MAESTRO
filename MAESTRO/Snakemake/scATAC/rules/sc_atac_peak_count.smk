_countpeak_threads = 4

rule scatac_countpeak:
    input:
        finalpeak = "Result/Analysis/{sample}/{sample}_final_peaks.bed",
        validbarcode = "Result/QC/{sample}/{sample}_scATAC_validcells.txt",
        frag = "Result/Mapping/{sample}/fragments_corrected_dedup_count.tsv"
    output:
        counts = "Result/Analysis/{sample}/{sample}_peak_count.h5"
    params:
        peakcount = temp("Result/Analysis/{sample}/{sample}_sorted_peakcount.tsv"),
        species = config["species"],
        outdir = "Result/Analysis/{sample}",
        outpre = "{sample}"
    threads:
        _countpeak_threads
    benchmark:
        "Result/Benchmark/{sample}_PeakCount.benchmark"
    shell:
        """
        LC_ALL=C fgrep -f {input.validbarcode} -w {input.frag} | bedtools intersect -sorted -wa -wb -a {input.finalpeak} -b - | sort -k1,1 -k2,2 -k3,3 -k7,7 - | bedtools groupby -g 1,2,3,7 -c 7 -o count > {params.peakcount};
        MAESTRO scatac-peakcount --binary --filtered_count {params.peakcount} --species {params.species} --directory {params.outdir} --outprefix {params.outpre}
        """

rule scatac_qcfilter:
    input:
        counts = "Result/Analysis/{sample}/{sample}_peak_count.h5"
    output:
        filtercount = "Result/QC/{sample}/{sample}_filtered_peak_count.h5"
    params:
        outdir = "Result/QC/{sample}",
        outpre = "{sample}",
        peak = config["cutoff"]["peak"],
        cell = config["cutoff"]["cell"]
    benchmark:
        "Result/Benchmark/{sample}_QCFilter.benchmark"
    shell:
        """
        MAESTRO scatac-qc --format h5 --peakcount {input.counts} --peak-cutoff {params.peak} --cell-cutoff {params.cell} \
        --directory {params.outdir} --outprefix {params.outpre}
        """

rule scatac_genescore:
    input:
        filtercount = "Result/QC/{sample}/{sample}_filtered_peak_count.h5",
        genebed = "%s/annotations/%s_ensembl.bed" %(SCRIPT_PATH, config["species"]),
    output:
        genescore = "Result/Analysis/{sample}/{sample}_gene_score.h5"
    params:
        genedistance = config["genedistance"],
        species = config["species"],
        outdir = "Result/Analysis/{sample}",
        outpre = "{sample}",
        rpmodel = config["rpmodel"]
    benchmark:
        "Result/Benchmark/{sample}_GeneScore.benchmark"
    shell:
        "MAESTRO scatac-genescore --format h5 --peakcount {input.filtercount} --species {params.species} --directory {params.outdir} --outprefix {params.outpre} --model {params.rpmodel}"
