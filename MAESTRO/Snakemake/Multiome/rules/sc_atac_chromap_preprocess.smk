_cat_threads= 2

if config["gzip"]:
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
else:
    rule scatac_link_fastq:
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
        shell:
            """
            cat {input.r1} > {output.r1cat}
            cat {input.r2} > {output.r2cat}
            cat {input.r3} > {output.r3cat}
            """
