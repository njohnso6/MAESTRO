_chromap_threads = 8

rule scatac_chromap:
    input:
        r1 = "Result/Tmp/{sample}/{sample}_R1.fastq",
        r3 = "Result/Tmp/{sample}/{sample}_R3.fastq",
        barcode = "Result/Tmp/{sample}/{sample}_R2.fastq",
    output:
        temp("Result/Mapping/{sample}/fragments_pre_corrected_dedup_count.tsv")
    log:
        "Result/Log/{sample}_chromap.log"
    benchmark:
        "Result/Benchmark/{sample}_chromap.benchmark"
    threads:
        _chromap_threads
    params:
        index = config["genome"]["index"],
        fasta = config["genome"]["fasta"],
        whitelist = config["whitelist"],
    shell:
        "chromap --preset atac -x {params.index} -r {params.fasta} -1 {input.r1} -2 {input.r3} -o {output} -b {input.barcode} -t {threads} --barcode-whitelist {params.whitelist} --chr-order " + SCRIPT_PATH + "/utils/chr_list_sorted.tsv 2> {log}"


rule scatac_process_frag:
    input:
        "Result/Mapping/{sample}/fragments_pre_corrected_dedup_count.tsv"
    output:
        frag = "Result/Mapping/{sample}/fragments_corrected_dedup_count.tsv",
        fraggz = "Result/Mapping/{sample}/fragments_corrected_dedup_count.tsv.gz"
    benchmark:
        "Result/Benchmark/{sample}_fragReshape.benchmark"
    shell:
        "python " + SCRIPT_PATH + "/scATAC_FragmentReshape.py -F {input} -O {output.frag};"
        "bgzip -c {output.frag} > {output.fraggz};"
        "tabix -p bed {output.fraggz}"
