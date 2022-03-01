_samtools_thead = 4
_qcplot_threads = 4

rule scatac_qcstat_mapped:
    input:
        frag_count = "Result/ATAC/Mapping/fragments_corrected_dedup_count.tsv",
    output:
        mapped_stat = temp("Result/ATAC/QC/singlecell_mapped.txt"),
    params:
        frag_count_sort = "Result/ATAC/Mapping/fragments_corrected_count_sortedbybarcode.tsv"
    benchmark:
        "Result/Benchmark/%s_scATAC_SingleQCMappability.benchmark" %(config["outprefix"])
    shell:
        """
        sort -k4,4 -V {input.frag_count} > {params.frag_count_sort}

        bedtools groupby -i {params.frag_count_sort} -g 4 -c 5 -o sum > {output.mapped_stat}
        """


rule scatac_qcstat_promoter:
    input:
        frag_count = "Result/ATAC/Mapping/fragments_corrected_dedup_count.tsv",
        promoter = "%s/annotations/%s_promoter.bed" %(SCRIPT_PATH, config["species"])
    output:
        promoter_stat = temp("Result/ATAC/QC/singlecell_promoter.txt"),
    params:
        frag_promoter = "Result/ATAC/Mapping/fragments_promoter.tsv",
        frag_promoter_sort = "Result/ATAC/Mapping/fragments_promoter_sortbybarcode.tsv",
    benchmark:
        "Result/Benchmark/%s_scATAC_SingleQCPromoter.benchmark" %(config["outprefix"])
    shell:
        """
        bedtools intersect -wa -a {input.frag_count} -b {input.promoter} -u > {params.frag_promoter}

        sort -k4,4 -V {params.frag_promoter} > {params.frag_promoter_sort}

        bedtools groupby -i {params.frag_promoter_sort} -g 4 -c 5 -o sum > {output.promoter_stat}
        """


rule scatac_qcstat_singlecell:
    input:
        mapped_stat = "Result/ATAC/QC/singlecell_mapped.txt",
        promoter_stat = "Result/ATAC/QC/singlecell_promoter.txt",
    output:
        single_stat = "Result/ATAC/QC/singlecell.txt"
    params:
        mapped_stat_sort = "Result/ATAC/QC/singlecell_mapped_sortbybarcode.txt",
        promoter_stat_sort = "Result/ATAC/QC/singlecell_promoter_sortbybarcode.txt"
    benchmark:
        "Result/Benchmark/%s_scATAC_SingleQCStat.benchmark" %(config["outprefix"])
    shell:
        """
        sort -k1,1 {input.mapped_stat} > {params.mapped_stat_sort}

        sort -k1,1 {input.promoter_stat} > {params.promoter_stat_sort}

        join --nocheck-order -t $'\t' -a1 -e'0' -o'1.1 1.2 2.2' -1 1 -2 1 {params.mapped_stat_sort} {params.promoter_stat_sort} > {output.single_stat}
        """


if config["format"] == "fastq" and config["mapping"] == "minimap2":
    rule scatac_qcstat_bulk:
        input:
            bam = "Result/ATAC/Mapping/%s.sortedByPos.rmdp.CBadded.bam" %(config["outprefix"]),
            promoter = "%s/annotations/%s_promoter.bed" %(SCRIPT_PATH, config["species"]),
            peak = "Result/ATAC/Analysis/%s_all_peaks.narrowPeak" %(config["outprefix"])
        output:
            bulk_stat = "Result/ATAC/QC/flagstat.txt",
            bam = "Result/ATAC/Mapping/%s.sortedByPos.rmdp.CBadded.unique.bam" %(config["outprefix"]),
            bed = "Result/ATAC/Mapping/%s.sortedByPos.rmdp.CBadded.unique.bed" %(config["outprefix"]),
        threads:
            _samtools_thead
        benchmark:
            "Result/Benchmark/%s_scATAC_BulkQCStat.benchmark" %(config["outprefix"])
        shell:
            """
            samtools flagstat --threads {threads} {input.bam} > {output.bulk_stat}
            samtools view -F 2316 -f 0x2 -q 30 -b -o {output.bam} {input.bam}
            samtools view {output.bam} -c >> {output.bulk_stat}
            bedtools bamtobed -i {output.bam} > {output.bed}
            grep 'chrM' {output.bed} -c >> {output.bulk_stat} || true
            grep -v 'chrM' {output.bed} | bedtools intersect -wa -a - -b {input.promoter} -u | wc -l >> {output.bulk_stat} || true
            grep -v 'chrM' {output.bed} | bedtools intersect -wa -a - -b {input.peak} -u | wc -l >> {output.bulk_stat} || true
            """

elif config["format"] == "fastq" and config["mapping"] == "chromap":
    rule scatac_qcstat_bulk:
        input:
            frag_dedup = "Result/ATAC/Mapping/fragments_corrected_dedup_count.tsv",
            promoter = "%s/annotations/%s_promoter.bed" %(SCRIPT_PATH, config["species"]),
            peak = "Result/ATAC/Analysis/%s_all_peaks.narrowPeak" %(config["outprefix"]),
            log = "Result/ATAC/Log/%s_chromap.log" %(config["outprefix"]),
        output:
            bulk_stat = "Result/ATAC/QC/flagstat.txt",
            fragbed = "Result/ATAC/QC/%s_frag.bed" %(config["outprefix"])
        threads:
            _samtools_thead
        benchmark:
            "Result/Benchmark/%s_scATAC_BulkQCStat.benchmark" %(config["outprefix"])
        shell:
            """
            awk '{{print $3-$2}}' {input.frag_dedup} > {output.fragbed}
            grep 'Number of reads' {input.log} | head -1 | cut -f4 -d' ' | cut -f1 -d'.' > {output.bulk_stat}
            grep 'Number of uniquely mapped reads' {input.log} | cut -f6 -d' ' | cut -f1 -d'.' >> {output.bulk_stat}
            grep 'chrM' {input.frag_dedup} -c >> {output.bulk_stat} || true
            grep -v 'chrM' {input.frag_dedup} | bedtools intersect -wa -a - -b {input.promoter} -u | wc -l >> {output.bulk_stat} || true
            grep -v 'chrM' {input.frag_dedup} | bedtools intersect -wa -a - -b {input.peak} -u | wc -l >> {output.bulk_stat} || true
            """

elif config["format"] == "fragments":
    rule scatac_qcstat_bulk:
        input:
            frag_dedup = "Result/ATAC/Mapping/fragments_corrected_dedup_count.tsv",
            promoter = "%s/annotations/%s_promoter.bed" %(SCRIPT_PATH, config["species"]),
            peak = "Result/ATAC/Analysis/%s_all_peaks.narrowPeak" %(config["outprefix"]),
        output:
            bulk_stat = "Result/ATAC/QC/flagstat.txt",
            fragbed = "Result/ATAC/QC/%s_frag.bed" %(config["outprefix"]),
        threads:
            _samtools_thead
        benchmark:
            "Result/Benchmark/%s_scATAC_BulkQCStat.benchmark" %(config["outprefix"])
        shell:
            """
            awk '{{print $3-$2}}' {input.frag_dedup} > {output.fragbed}
            grep 'chrM' {input.frag_dedup} -c >> {output.bulk_stat} || true
            grep -v 'chrM' {input.frag_dedup} | bedtools intersect -wa -a - -b {input.promoter} -u | wc -l >> {output.bulk_stat} || true
            grep -v 'chrM' {input.frag_dedup} | bedtools intersect -wa -a - -b {input.peak} -u | wc -l >> {output.bulk_stat} || true
            """

if config["format"] == "fastq":
    rule scatac_qcplot:
        input:
            fragbed = "Result/ATAC/QC/%s_frag.bed" %(config["outprefix"]),
            single_stat = "Result/ATAC/QC/singlecell.txt",
            bulk_stat = "Result/ATAC/QC/flagstat.txt"
        output:
            readdistr = "Result/ATAC/QC/%s_scATAC_read_distr.png" %(config["outprefix"]),
            qcfrag = "Result/ATAC/QC/%s_scATAC_fragment_size.png" %(config["outprefix"]),
            qcfrip = "Result/ATAC/QC/%s_scATAC_cell_filtering.png" %(config["outprefix"]),
            validbarcode = "Result/ATAC/QC/%s_scATAC_validcells.txt" %(config["outprefix"]),
        params:
            outdir = "Result/ATAC/QC/",
            outpre = config["outprefix"],
            fragbed = "%s_frag.bed" %(config["outprefix"]),
            single_stat = "singlecell.txt",
            bulk_stat = "flagstat.txt",
            counts = config["cutoff"]["atac_count"],
            frip = config["cutoff"]["atac_frip"],
            mapping = config["mapping"],
        threads:
            _qcplot_threads
        benchmark:
            "Result/Benchmark/%s_scATAC_QCPlot.benchmark" %(config["outprefix"])
        shell:
            "Rscript " + RSCRIPT_PATH + "/scATACseq_qc.R --bulkstat {params.bulk_stat} --fragment {params.fragbed} --mapping {params.mapping} --singlestat {params.single_stat} "
            "--countcutoff {params.counts} --fripcutoff {params.frip} --prefix {params.outpre} --outdir {params.outdir}"

elif config["format"] == "fragments":
    rule scatac_qcplot:
        input:
            fragbed = "Result/ATAC/QC/%s_frag.bed" %(config["outprefix"]),
            single_stat = "Result/ATAC/QC/singlecell.txt",
        output:
            qcfrag = "Result/ATAC/QC/%s_scATAC_fragment_size.png" %(config["outprefix"]),
            qcfrip = "Result/ATAC/QC/%s_scATAC_cell_filtering.png" %(config["outprefix"]),
            validbarcode = "Result/ATAC/QC/%s_scATAC_validcells.txt" %(config["outprefix"]),
        params:
            outdir = "Result/ATAC/QC/",
            outpre = config["outprefix"],
            fragbed = "%s_frag.bed" %(config["outprefix"]),
            single_stat = "singlecell.txt",
            counts = config["cutoff"]["count"],
            frip = config["cutoff"]["frip"]
        threads:
            _qcplot_threads
        benchmark:
            "Result/Benchmark/%s_scATAC_QCPlot.benchmark" %(config["outprefix"]),
        shell:
            "Rscript " + RSCRIPT_PATH + "/scATACseq_qc.R --fragment {params.fragbed} --singlestat {params.single_stat} "
            "--countcutoff {params.counts} --fripcutoff {params.frip} --prefix {params.outpre} --outdir {params.outdir}"
