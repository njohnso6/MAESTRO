

STAR_VERSION = subprocess.check_output("STAR --version", shell=True)

rule STARsolo:
    """
    STARsolo to align single cell or single nuceli data
    specify --soloFeatures Gene for single-cell  data
    specify --soloFeatures GeneFull for single-nuclei data
    specify --soloFeatures Gene GeneFull for getting both counts in exons level and exon + intron level (velocity)
    """
    input:
        mapindex = config["genome"]["rna_mapindex"],
        whitelist = config["barcode"]["rna_whitelist"]
    output:
        bam = "Result/RNA/STAR/%sAligned.sortedByCoord.out.bam" %(config["outprefix"]),
        bai = "Result/RNA/STAR/%sAligned.sortedByCoord.out.bam.bai" %(config["outprefix"]),
        rawmtx = "Result/RNA/STAR/%sSolo.out/%s/raw/matrix.mtx" %(config["outprefix"],config["STARsolo_Features"].split(" ")[0]),
        feature = "Result/RNA/STAR/%sSolo.out/%s/raw/features.tsv" %(config["outprefix"],config["STARsolo_Features"].split(" ")[0]),
        barcode = "Result/RNA/STAR/%sSolo.out/%s/raw/barcodes.tsv" %(config["outprefix"],config["STARsolo_Features"].split(" ")[0])
    params:
        star_custom = config["STARsolo_Features"],
        outprefix = "Result/RNA/STAR/%s"%(config["outprefix"]),
        transcript = lambda wildcards: ','.join(FILES["RNA"]["R2"]),
        barcode = lambda wildcards: ','.join(FILES["RNA"]["R1"]),
        barcodestart = config["barcode"]["barcodestart"],
        barcodelength = config["barcode"]["barcodelength"],
        umistart = config["barcode"]["umistart"],
        umilength = config["barcode"]["umilength"]
    version: STAR_VERSION
    log:
        "Result/RNA/Log/%s_STAR.log" %(config["outprefix"])
    benchmark:
        "Result/Benchmark/%s_scRNA_STAR.benchmark" %(config["outprefix"])
    threads:
        config["cores"]
    shell:
        """
        STAR \
            --runMode alignReads \
            --genomeDir {input.mapindex} \
            --runThreadN {threads} \
            --outFileNamePrefix {params.outprefix} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
            --soloType CB_UMI_Simple \
            --soloFeatures {params.star_custom} \
            --soloCBwhitelist {input.whitelist} \
            --soloCBstart {params.barcodestart} \
            --soloCBlen {params.barcodelength} \
            --soloUMIstart {params.umistart} \
            --soloUMIlen {params.umilength} \
            --soloCBmatchWLtype 1MM_multi_pseudocounts \
            --soloUMIfiltering MultiGeneUMI \
            --readFilesIn {params.transcript} {params.barcode} \
			--readFilesCommand zcat \
            --genomeSAindexNbases 2 \
            --limitOutSJcollapsed 5000000 \
            > {log} 2>&1

        samtools index -b -@ {threads} {output.bam} >> {log} 2>&1
        """
