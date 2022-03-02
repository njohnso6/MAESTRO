
if config["barcode"]["atac_whitelist"]:
    rule scatac_barcodecorrect:
        input:
            r2 = "Result/ATAC/Tmp/%s_R2.fastq" %(config["outprefix"]),
            whitelist_atac = config["barcode"]["atac_whitelist"],
            whitelist_rna = config["barcode"]["rna_whitelist"]
        output:
            bc_correct = "Result/ATAC/Mapping/barcode_correct.txt",
            bc_correct_uniq = "Result/ATAC/Mapping/barcode_correct_uniq.txt"
        params:
            outdir = "Result/ATAC/Mapping/"
        benchmark:
            "Result/Benchmark/%s_scATAC_BarcodeCorrect.benchmark" %(config["outprefix"])
        shell:
            "python " + SCRIPT_PATH + "/scATAC_10x_BarcodeCorrect.py -b {input.r2} -B {input.whitelist_atac} --barcodelib-rna {input.whitelist_rna} -O {params.outdir};"
            "sort -k1,1 -k3,3 {output.bc_correct} | uniq > {output.bc_correct_uniq}"


else:
    rule scatac_barcodecorrect:
        input:
            r2 = "Result/ATAC/Tmp/%s_R2.fastq" %(config["outprefix"]),
        output:
            bc_correct = "Result/ATAC/Mapping/barcode_correct.txt",
            bc_correct_uniq = "Result/ATAC/Mapping/barcode_correct_uniq.txt"
        params:
            outdir = "Result/ATAC/Mapping/%s"
        benchmark:
            "Result/Benchmark/%s_scATAC_BarcodeCorrect.benchmark" %(config["outprefix"])
        shell:
            "python " + SCRIPT_PATH + "/scATAC_10x_BarcodeCorrect.py -b {input.r2} -O {params.outdir};"
            "sort -k1,1 -k3,3 {output.bc_correct} | uniq > {output.bc_correct_uniq}"


rule scatac_fragmentcorrect:
    input:
        fragments = "Result/ATAC/Mapping/fragments.tsv",
        bc_correct = "Result/ATAC/Mapping/barcode_correct.txt",
    output:
        frag_count = "Result/ATAC/Mapping/fragments_corrected_count.tsv",
        frag_sort = temp("Result/ATAC/Mapping/fragments_corrected_sorted.tsv"),
    params:
        outdir = "Result/ATAC/Mapping/",
        frag_correct = "Result/ATAC/Mapping/fragments_corrected.tsv",
    benchmark:
        "Result/Benchmark/%s_scATAC_FragCorrect.benchmark" %(config["outprefix"])
    shell:
        "python " + SCRIPT_PATH + "/scATAC_FragmentCorrect.py -F {input.fragments} -C {input.bc_correct} -O {params.outdir};"
        "sort -k1,1 -k2,2 -k3,3 -k4,4 -V {params.frag_correct} > {output.frag_sort};"
        "bedtools groupby -i {output.frag_sort} -g 1,2,3,4 -c 4 -o count | sort -k1,1 -k2,2n > {output.frag_count};"
        "rm {params.frag_correct}"
