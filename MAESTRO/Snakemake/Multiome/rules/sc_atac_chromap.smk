_chromap_threads = 8

rule scatac_chromap:
    input:
        r1 = "Result/ATAC/Tmp/%s_R1.fastq" %(config["outprefix"]),
        r3 = "Result/ATAC/Tmp/%s_R3.fastq" %(config["outprefix"]),
        barcode = "Result/ATAC/Tmp/%s_R2.fastq" %(config["outprefix"]),
    output:
        temp("Result/ATAC/Mapping/fragments_dedup_count.tsv")
    log:
        "Result/ATAC/Log/%s_chromap.log" %(config["outprefix"])
    benchmark:
        "Result/Benchmark/%s_scATAC_chromap.benchmark" %(config["outprefix"])
    threads:
        _chromap_threads
    params:
        index = config["genome"]["atac_mapindex"],
        fasta = config["genome"]["atac_fasta"],
        whitelist = config["barcode"]["atac_whitelist"],
    shell:
        "chromap --preset atac -x {params.index} -r {params.fasta} -1 {input.r1} -2 {input.r3} -o {output} -b {input.barcode} -t {threads} --barcode-whitelist {params.whitelist} --chr-order " + SCRIPT_PATH + "/utils/chr_list_sorted.tsv 2> {log}"

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
        fragments = "Result/ATAC/Mapping/fragments_dedup_count.tsv",
        bc_correct = "Result/ATAC/Mapping/barcode_correct.txt",
    output:
        "Result/ATAC/Mapping/fragments_pre_corrected_dedup_count.tsv",
    params:
        outdir = "Result/ATAC/Mapping/",
        suffix = "fragments_pre_corrected_dedup_count.tsv"
    benchmark:
        "Result/Benchmark/%s_scATAC_FragCorrect.benchmark" %(config["outprefix"])
    shell:
        "python " + SCRIPT_PATH + "/scATAC_FragmentCorrect.py -F {input.fragments} -C {input.bc_correct} -O {params.outdir} -S {params.suffix};"


rule scatac_process_frag:
    input:
        "Result/ATAC/Mapping/fragments_pre_corrected_dedup_count.tsv"
    output:
        frag = "Result/ATAC/Mapping/fragments_corrected_dedup_count.tsv",
        fraggz = "Result/ATAC/Mapping/fragments_corrected_dedup_count.tsv.gz"
    benchmark:
        "Result/Benchmark/%s_scATAC_fragReshape.benchmark" %(config["outprefix"])
    shell:
        "python " + SCRIPT_PATH + "/scATAC_FragmentReshape.py -F {input} -O {output.frag};"
        "bgzip -c {output.frag} > {output.fraggz};"
        "tabix -p bed {output.fraggz}"
