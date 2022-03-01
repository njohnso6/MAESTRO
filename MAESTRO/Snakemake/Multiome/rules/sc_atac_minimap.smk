_minmap_threads = 8
_picard_threads = 4
_bamAddCB_threads = 4
_bamindex_threads = 4

rule scatac_minimap:
    input:
        fasta = config["genome"]["atac_fasta"],
        r1 = "Result/ATAC/Tmp/%s_R1.barcoded.fastq" %(config["outprefix"]),
        r3 = "Result/ATAC/Tmp/%s_R3.barcoded.fastq" %(config["outprefix"]),
    output:
        bam = temp("Result/ATAC/Mapping/%s.sortedByPos.bam" %(config["outprefix"]))
    threads:
        _minmap_threads
    benchmark:
        "Result/Benchmark/%s_scATAC_Minimap2.benchmark" %(config["outprefix"])
    shell:
        """
        minimap2 -ax sr -t {threads} {input.fasta} {input.r1} {input.r3} \
        | samtools view --threads {threads} -b \
        | samtools sort --threads {threads} -o {output.bam}
        """

rule scatac_fragmentgenerate:
    input:
        bam = "Result/ATAC/Mapping/%s.sortedByPos.bam" %(config["outprefix"])
    output:
        fragments = "Result/ATAC/Mapping/fragments.tsv",
        bam = "Result/ATAC/Mapping/%s.sortedByPos.CRadded.bam" %(config["outprefix"]),
    params:
        outdir = "Result/ATAC/Mapping/
    benchmark:
        "Result/Benchmark/%s_scATAC_FragGenerate.benchmark" %(config["outprefix"])
    shell:
        "python " + SCRIPT_PATH + "/scATAC_FragmentGenerate.py -B {input.bam} -O {params.outdir} --addtag CR"

rule scatac_rmdp:
    input:
        bam = "Result/ATAC/Mapping/%s.sortedByPos.CRadded.bam" %(config["outprefix"]),
    output:
        bam = "Result/ATAC/Mapping/%s.sortedByPos.CRadded.rmdp.bam" %(config["outprefix"]),
        metric = "Result/ATAC/Mapping/%s.rmdp.txt" %(config["outprefix"]),
        fragbed = "Result/ATAC/QC/%s_frag.bed" %(config["outprefix"]),
    params:
        sam = "Result/ATAC/Mapping/%s.sortedByPos.CRadded.rmdp.sample.sam" %(config["outprefix"]),
    threads:
        _picard_threads
    benchmark:
        "Result/Benchmark/%s_scATAC_Rmdp.benchmark" %(config["outprefix"])
    shell:
        """
        picard MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metric} TMP_DIR=Result/ATAC/Tmp

        samtools view -@ {threads} -s 0.01 -o {params.sam} {input.bam}

        awk '{{if ($9>0) print $9}}' {params.sam} > {output.fragbed}
        """

rule scatac_bamaddCB:
        input:
            bam = "Result/ATAC/Mapping/%s.sortedByPos.CRadded.rmdp.bam" %(config["outprefix"]),
            bc_correct = "Result/ATAC/Mapping/barcode_correct_uniq.txt"
        output:
            bam = "Result/ATAC/Mapping/%s.sortedByPos.rmdp.CBadded.bam" %(config["outprefix"])
        params:
            outdir = "Result/ATAC/Mapping/,
            outprefix = "%s.sortedByPos.rmdp.CBadded" %(config["outprefix"])
        threads:
            _bamAddCB_threads
        benchmark:
            "Result/Benchmark/%s_scATAC_BamAddCB.benchmark" %(config["outprefix"])
        shell:
            "python " +  SCRIPT_PATH + "/scATAC_BamAddTag.py -B {input.bam} -T {input.bc_correct} -C CR "
            "-O {params.outdir} -P {params.outprefix}"


rule scatac_bamindex:
    input:
        bam = "Result/ATAC/Mapping/%s.sortedByPos.rmdp.CBadded.bam" %(config["outprefix"]),
    output:
        bai = "Result/ATAC/Mapping/%s.sortedByPos.rmdp.CBadded.bam.bai" %(config["outprefix"]),
    threads:
        _bamindex_threads
    benchmark:
        "Result/Benchmark/%s_scATAC_BamIndex.benchmark" %(config["outprefix"])
    shell:
        "samtools index -@ {threads} {input.bam}"
