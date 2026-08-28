velocyto run -b ./CRC1/result_bam/barcodes.tsv.gz \
        -o ./CRC1/result_bam/velocyto_loom/ \
        -m ./ref/RhapRef_Mouse_WTA/mm39_rmsk.gtf \
        ./CRC1/result_bam/UB_Combined_WTA_CRC1.BAM \
        ./ref/RhapRef_Mouse_WTA/gencode.vM31.primary_assembly.annotation-filtered.gtf

velocyto run -b ./CRC2/result_bam/barcodes.tsv.gz \
        -o ./CRC2/result_bam/velocyto_loom/ \
        -m ./ref/RhapRef_Mouse_WTA/mm39_rmsk.gtf \
        ./CRC2/result_bam/UB_Combined_WTA_CRC2.BAM \
        ./ref/RhapRef_Mouse_WTA/gencode.vM31.primary_assembly.annotation-filtered.gtf

velocyto run -b ./CRC3/result_bam/barcodes.tsv.gz \
        -o ./CRC3/result_bam/velocyto_loom/ \
        -m ./ref/RhapRef_Mouse_WTA/mm39_rmsk.gtf \
        ./CRC3/result_bam/UB_Combined_WTA_CRC3.BAM \
        ./ref/RhapRef_Mouse_WTA/gencode.vM31.primary_assembly.annotation-filtered.gtf

velocyto run -b ./CRC4/result_bam/barcodes.tsv.gz \
        -o ./CRC4/result_bam/velocyto_loom/ \
        -m ./ref/RhapRef_Mouse_WTA/mm39_rmsk.gtf \
        ./CRC4/result_bam/UB_Combined_WTA_CRC4.BAM \
        ./ref/RhapRef_Mouse_WTA/gencode.vM31.primary_assembly.annotation-filtered.gtf

velocyto run -b ./CRC5/result_bam/barcodes.tsv.gz \
        -o ./CRC5/result_bam/velocyto_loom/ \
        -m ./ref/RhapRef_Mouse_WTA/mm39_rmsk.gtf \
        ./CRC5/result_bam/UB_Combined_WTA_CRC5.BAM \
        ./ref/RhapRef_Mouse_WTA/gencode.vM31.primary_assembly.annotation-filtered.gtf               