_shortfragment_threads = 2

macs2_genome = "hs" if config["species"] == "GRCh38" else "mm"

rule scatac_allpeakcall:
    input:
        frag_dedup = "Result/ATAC/Mapping/fragments_corrected_dedup_count.tsv"
    output:
        peak = "Result/ATAC/Analysis/%s_all_peaks.narrowPeak" %(config["outprefix"]),
        bdg = "Result/ATAC/Analysis/%s_all_treat_pileup.bdg" %(config["outprefix"]),
    params:
        name = "%s_all" %(config["outprefix"]),
        genome = macs2_genome
    log:
        "Result/ATAC/Log/%s_macs2_allpeak.log" %(config["outprefix"])
    benchmark:
        "Result/Benchmark/%s_scATAC_AllPeakCall.benchmark" %(config["outprefix"])
    shell:
        "macs2 callpeak -f BEDPE -g {params.genome} --outdir Result/ATAC/Analysis/ -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t {input.frag_dedup}"

if config["shortpeaks"] and not config["custompeaks"] :
    rule scatac_shortfragment:
        input:
            frag_dedup = "Result/ATAC/Mapping/fragments_corrected_dedup_count.tsv"
        output:
            frag_short = "Result/ATAC/Mapping/fragments_corrected_150bp.tsv"
        threads:
            _shortfragment_threads
        benchmark:
            "Result/Benchmark/%s_scATAC_ShortFrag.benchmark" %(config["outprefix"])
        shell:
            "awk -F'\\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{if (abs($3-$2)<=150) print}}' {input.frag_dedup} > {output.frag_short}"

    rule scatac_shortpeakcall:
        input:
            frag_short = "Result/ATAC/Mapping/fragments_corrected_150bp.tsv"
        output:
            bed = "Result/ATAC/Analysis/%s_150bp_peaks.narrowPeak" %(config["outprefix"])
        params:
            name = "%s_150bp" %(config["outprefix"]),
            genome = macs2_genome
        log:
            "Result/ATAC/Log/%s_macs2_shortpeak.log" %(config["outprefix"])
        benchmark:
            "Result/Benchmark/%s_scATAC_ShortPeakCall.benchmark" %(config["outprefix"])
        shell:
            "macs2 callpeak -f BEDPE -g {params.genome} --outdir Result/ATAC/Analysis/{wildcards.sample} -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t {input.frag_short}"

    rule scatac_mergepeak:
        input:
            allpeak = "Result/ATAC/Analysis/%s_all_peaks.narrowPeak" %(config["outprefix"]),
            shortpeak = "Result/ATAC/Analysis/%s_150bp_peaks.narrowPeak" %(config["outprefix"]),
        output:
            finalpeak = "Result/ATAC/Analysis/%s_final_peaks.bed" %(config["outprefix"])
        params:
            catpeaksort = "Result/ATAC/Analysis/%s_cat_peaks.bed" %(config["outprefix"])
        benchmark:
            "Result/Benchmark/%s_scATAC_PeakMerge.benchmark" %(config["outprefix"])
        shell:
            """
            cat {input.allpeak} {input.shortpeak} \
            | sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort}

            mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak}

            rm {params.catpeaksort}
            """

elif config["custompeaks"] and config["shortpeaks"]:
    rule scatac_shortfragment:
        input:
            frag_dedup = "Result/ATAC/Mapping/fragments_corrected_dedup_count.tsv"
        output:
            frag_short = "Result/ATAC/Mapping/fragments_corrected_150bp.tsv"
        threads:
            _shortfragment_threads
        benchmark:
            "Result/Benchmark/%s_scATAC_ShortFrag.benchmark" %(config["outprefix"])
        shell:
            "awk -F'\\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{if (abs($3-$2)<=150) print}}' {input.frag_dedup} > {output.frag_short}"

    rule scatac_shortpeakcall:
        input:
            frag_short = "Result/ATAC/Mapping/fragments_corrected_150bp.tsv"
        output:
            bed = "Result/ATAC/Analysis/%s_150bp_peaks.narrowPeak" %(config["outprefix"])
        params:
            name = "%s_150bp" %(config["outprefix"]),
            genome = macs2_genome
        log:
            "Result/ATAC/Log/%s_macs2_shortpeak.log" %(config["outprefix"])
        benchmark:
            "Result/Benchmark/%s_scATAC_ShortPeakCall.benchmark" %(config["outprefix"])
        shell:
            "macs2 callpeak -f BEDPE -g {params.genome} --outdir Result/ATAC/Analysis/{wildcards.sample} -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t {input.frag_short}"

    rule scatac_mergepeak:
        input:
            allpeak = "Result/ATAC/Analysis/%s_all_peaks.narrowPeak" %(config["outprefix"]),
            shortpeak = "Result/ATAC/Analysis/%s_150bp_peaks.narrowPeak" %(config["outprefix"]),
            custompeak = config["custompeaksloc"],
        output:
            finalpeak = "Result/ATAC/Analysis/%s_final_peaks.bed" %(config["outprefix"])
        params:
            catpeaksort = "Result/ATAC/Analysis/%s_cat_peaks.bed" %(config["outprefix"])
        benchmark:
            "Result/Benchmark/%s_scATAC_PeakMerge.benchmark" %(config["outprefix"])
        shell:
            """
            cat {input.allpeak} {input.shortpeak} {input.custompeak} \
            | sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort}

            mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak}

            rm {params.catpeaksort}
            """

elif config["custompeaks"] and not config["shortpeaks"]:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/ATAC/Analysis/%s_all_peaks.narrowPeak" %(config["outprefix"]),
            custompeaks = config["custompeaksloc"]
        output:
            finalpeak = "Result/ATAC/Analysis/%s_final_peaks.bed" %(config["outprefix"])
        params:
            catpeaksort = "Result/ATAC/Analysis/%s_cat_peaks.bed" %(config["outprefix"])
        benchmark:
            "Result/Benchmark/%s_scATAC_PeakMerge.benchmark" %(config["outprefix"])
        shell:
            """
            cat {input.allpeak} {input.custompeaks} \
            | sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort}

            mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak}

            rm {params.catpeaksort}
            """

else:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/ATAC/Analysis/%s_all_peaks.narrowPeak" %(config["outprefix"])
        output:
            finalpeak = "Result/ATAC/Analysis/%s_final_peaks.bed" %(config["outprefix"])
        params:
            catpeaksort = "Result/ATAC/Analysis/%s_cat_peaks.bed" %(config["outprefix"])
        benchmark:
            "Result/Benchmark/%s_scATAC_PeakMerge.benchmark" %(config["outprefix"])
        shell:
            """
            cat {input.allpeak} \
            | sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort}

            mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak}

            rm {params.catpeaksort}
            """

rule scatac_bdg2bw:
    input:
        bdg = "Result/ATAC/Analysis/%s_all_treat_pileup.bdg" %(config["outprefix"])
    output:
        bw = "Result/ATAC/Analysis/%s.bw" %(config["outprefix"])
    benchmark:
        "Result/Benchmark/%s_scATAC_Bdg2Bw.benchmark" %(config["outprefix"])
    params:
        chrom_len = "%s/annotations/%s_chrom_len.txt" %(SCRIPT_PATH, config["species"]),
        sort_bdg = "Result/ATAC/Analysis/%s_sort.bdg" %(config["outprefix"]),
        clip_bdg = "Result/ATAC/Analysis/%s_sort.clip" %(config["outprefix"]),
    shell:
        # https://gist.github.com/taoliu/2469050  bdg file generated from MACS2 can exceed chromosome limits
        "bedtools slop -i {input.bdg} -g {params.chrom_len} -b 0 | " + SCRIPT_PATH + "/utils/bedClip stdin {params.chrom_len} {params.clip_bdg};"
        "sort -k1,1 -k2,2n {params.clip_bdg} > {params.sort_bdg};"
        "" + SCRIPT_PATH + "/utils/bedGraphToBigWig {params.sort_bdg} {params.chrom_len} {output.bw};"
        "rm {params.sort_bdg};"
        "rm {params.clip_bdg}"
