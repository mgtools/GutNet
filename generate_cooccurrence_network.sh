#!/bin/bash

LIB=$1
WD=/N/project/BacInteraction/DiseaseNetwork_w_phyloseq/


declare -a arr=( "AA" "ACVD" "CD" "CRC" "WT" "OB" "OW" "RA" "T2D" "UC" )

# 1) quantification - download SRA FASTQ and run bracken 
while read -r $ERR <&3; do
    ${LIB}/download_sra_quantify_bracken.fixed.sh ${ERR}
done

# 2) get id to taxa annotations
Rscript lib/get_taxa_annotations.R

# 3) merge bracken outputs into phyloseq object
kraken-biom ../bracken_quantification/*/*s.breport -o sung_bracken.biom --fmt hdf5

# 4) filter bracken results 
Rscript lib/get_bracken_by_disease.R -i 0_mk_phyloseq/sung_bracken.biom -m metadata.filtered.csv -o 3_abundance_matrix/ --ignore 2_samples2remove/ignore_these_sra.txt

# 5) calculate spiec-easi - species wise
for i in "${arr[@]}";do
    Rscript ${WD}/lib/run_spiec_easi.mb.R ${WD}/3_abundance_matrix/${i}_bracken_filtered.s.csv ${WD}/4_spiec_easi_out/${i}_spiec.s.RData ${WD}/4_spiec_easi_out/${i}_spiec_mb.s.adj
done

# 6) make GML network files
for i in "${arr[@]}";do
    # build spiec-easi - species
    echo "${i}: se - s"
    python ${WD}/lib/build_network_gml.py -i ${WD}/4_spiec_easi_out/${i}_spiec_mb.s.adj -t spiec-easi -o ${WD}/6_mk_network_files/${i}_spiec-easi.filtered.annotated.gml --annotations ${WD}/0_taxa_annotations/phyloseq_id_to_annotation.csv
done
