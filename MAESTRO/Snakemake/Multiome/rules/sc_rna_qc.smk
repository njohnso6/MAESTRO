rule scrna_qc:
    input:
        rawmtx = "Result/RNA/STAR/%sSolo.out/%s/raw/matrix.mtx" %(config["outprefix"],config["STARsolo_Features"].split(" ")[0]),
        feature = "Result/RNA/STAR/%sSolo.out/%s/raw/features.tsv" %(config["outprefix"],config["STARsolo_Features"].split(" ")[0]),
        barcode = "Result/RNA/STAR/%sSolo.out/%s/raw/barcodes.tsv" %(config["outprefix"],config["STARsolo_Features"].split(" ")[0]),
    output:
        countgene = "Result/RNA/QC/%s_count_gene_stat.txt" %(config["outprefix"]),
        filtermatrix = "Result/RNA/QC/%s_filtered_gene_count.h5" %(config["outprefix"]),
        rnafilterplot = "Result/RNA/QC/%s_scRNA_cell_filtering.png" %(config["outprefix"]),
    params:
        counts = config["cutoff"]["rna_count"],
        gene = config["cutoff"]["rna_gene"],
        outpre = config["outprefix"],
        outdir = "Result/RNA/QC",
        species = config["species"]
    benchmark:
        "Result/Benchmark/%s_scRNA_QC.benchmark" %(config["outprefix"])
    shell:
        """
        MAESTRO scrna-qc \
            --format mtx \
            --matrix {input.rawmtx} \
            --feature {input.feature} \
            --barcode {input.barcode} \
            --species {params.species} \
            --count-cutoff {params.counts} \
            --gene-cutoff {params.gene} \
            --directory {params.outdir} \
            --outprefix {params.outpre}
        """
