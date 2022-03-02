_cat_threads= 2

rule scatac_preprocess:
    input:
        r1 = lambda wildcards: FILES["ATAC"]["R1"],
        r2 = lambda wildcards: FILES["ATAC"]["R2"],
        r3 = lambda wildcards: FILES["ATAC"]["R3"],
    output:
        r1cat = temp("Result/ATAC/Tmp/%s_R1.fastq" %(config["outprefix"])),
        r2cat = temp("Result/ATAC/Tmp/%s_R2.fastq" %(config["outprefix"])),
        r3cat = temp("Result/ATAC/Tmp/%s_R3.fastq" %(config["outprefix"])),
    log:
        "Result/ATAC/Log/%s_preprocess.log" %(config["outprefix"])
    benchmark:
        "Result/Benchmark/%s_scATAC_Preprocess.benchmark" %(config["outprefix"])
    threads: _cat_threads
    shell:
        """
        gunzip -c {input.r1} > {output.r1cat} 2> {log}
        gunzip -c {input.r2} > {output.r2cat} 2>> {log}
        gunzip -c {input.r3} > {output.r3cat} 2>> {log}

        """

rule scatac_fqaddbarcode:
    input:
        r1 = "Result/ATAC/Tmp/%s_R1.fastq" %(config["outprefix"]),
        r2 = "Result/ATAC/Tmp/%s_R2.fastq" %(config["outprefix"]),
        r3 = "Result/ATAC/Tmp/%s_R3.fastq" %(config["outprefix"]),
    output:
        r1 = temp("Result/ATAC/Tmp/%s_R1.barcoded.fastq"),
        r3 = temp("Result/ATAC/Tmp/%s_R3.barcoded.fastq"),
    benchmark:
        "Result/Benchmark/%s_scATAC_FqAddbarcode.benchmark" %(config["outprefix"])
    shell:
        """
        base=`head -n 2 {input.r2} | tail -n 1 | wc -L`

        sinto barcode --barcode_fastq {input.r2} --read1 {input.r1} --read2 {input.r3} -b $base
        """
