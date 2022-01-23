#!/bin/bash

OUT=bracken_quantification/
LOG=metadata.log
WD=$(pwd)
HG38=GRCh38_noalt_as/GRCh38_noalt_as
TEMP=/tmp
QUERY=$1
MINLEN=60
THREADS=12
KRAKENDB=kraken_db/k2_standard_20201202/

echo "$QUERY"

# get accession
SEARCH=$(esearch -db sra -query "${QUERY}" | efetch -format runinfo | grep "WGS" | grep "ILLUMINA" | sort -t "," -rn -k8,8 -k6,6 -k4,4 -k17,17)
ERR=$(echo ${SEARCH} | cut -d ',' -f 1 | head -n 1)
ENDEDNESS=$(echo ${SEARCH} | cut -d ',' -f 16 | head -n 1)

echo $ERR
#make outdir
OUTDIR=${OUT}/${QUERY}/
mkdir -p ${OUTDIR}

# if single end
if [ ${ENDEDNESS} == "SINGLE" ]; then

    echo "${ERR} is single-end"
    read_type='SE'

    #download
    echo "Downloading: ${ERR}"
    fasterq-dump ${ERR} --skip-technical --outdir ${OUTDIR}/ --outfile ${QUERY}.fastq -t ${TEMP} --threads ${THREADS}
    echo "Download Complete: ${ERR}"

    # get number of reads
    raw_read_count=$(wc -l ${OUTDIR}/${QUERY}.fastq | awk '{print $1}')
    raw_read_count=$(expr "${raw_read_count}" / 4)

    # trim fastqc
    if [ -f "${OUTDIR}/${QUERY}.fastq" ]; then
        echo "Trimming: ${ERR}"
        trimmomatic SE -threads ${THREADS} -phred33 ${OUTDIR}/${QUERY}.fastq ${OUTDIR}/${QUERY}.trimmed.fastq.gz ILLUMINACLIP:/N/u/tjlam/Quartz/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:${MINLEN}
        echo "Trimming Complete: ${ERR}"
    else
        echo "FASTQ missing: ${ERR}"
        exit 1 
    fi

    # get number of reads
    trimmed_read_count=$(gzip -dc ${OUTDIR}/${QUERY}.trimmed.fastq.gz | wc -l)
    trimmed_read_count=$(expr "${trimmed_read_count}" / 4)

    # delete raw
    rm ${OUTDIR}/${QUERY}.fastq

    # map to hg38 and remove host
    if [ -f "${OUTDIR}/${QUERY}.trimmed.fastq.gz" ]; then
        echo "Running Bowtie - Removing Host: ${ERR}"
        bowtie2 -p ${THREADS} -x ${HG38} -U ${OUTDIR}/${QUERY}.trimmed.fastq.gz --un-gz ${OUTDIR}/${QUERY}_host_removed > ${OUTDIR}/${QUERY}_mapped_and_unmapped.sam
        echo "Bowtie Complete"
    else
        echo "Trimmed missing: ${ERR}"
        exit 1
    fi

    # rename host removed
    mv ${OUTDIR}/${QUERY}_host_removed ${OUTDIR}/${QUERY}_host_removed.fastq.gz

    # remove sam file
    rm ${OUTDIR}/${QUERY}_mapped_and_unmapped.sam

    # remove trimmed
    rm ${OUTDIR}/${QUERY}.trimmed.fastq.gz 

    # count host removed
    host_rm_count=$(gzip -dc ${OUTDIR}/${QUERY}_host_removed.fastq.gz | wc -l)
    host_rm_count=$(expr "${host_rm_count}" / 4)

    # run kraken
    if [ -f "${OUTDIR}/${QUERY}_host_removed.fastq.gz" ]; then
        echo "Running Kraken: ${ERR}"
        kraken2 --db ${KRAKENDB} --thread ${THREADS} --report ${OUTDIR}/${QUERY}.kreport --out ${OUTDIR}/${QUERY}.kraken ${OUTDIR}/${QUERY}_host_removed.fastq.gz
        echo "Kraken Complete"
    else
        echo "host rm missing: ${ERR}"
        exit 1
    fi

    # remove host removed fasta
    rm ${OUTDIR}/${QUERY}_host_removed.fastq.gz

# if paired-end
elif [ ${ENDEDNESS} == "PAIRED" ]; then
    echo "${ERR} is paired-end"
    read_type='PE'

    # download    
    echo "Downloading: ${ERR}"
    fasterq-dump ${ERR} --skip-technical --split-files --outdir ${OUTDIR} --outfile ${QUERY} -t ${TEMP} -e ${THREADS}
    echo "Download Complete: ${ERR}"

    # get number of reads
    raw_read_count=$(wc -l ${OUTDIR}/${QUERY}_1.fastq | awk '{print $1}')
    raw_read_count=$(expr "${raw_read_count}" / 4)

    # trim fastqc
    if [ -f "${OUTDIR}/${QUERY}_1.fastq" ]; then
        echo "Trimming: ${ERR}"
        trimmomatic PE -threads ${THREADS} -phred33 ${OUTDIR}/${QUERY}_1.fastq ${OUTDIR}/${QUERY}_2.fastq ${OUTDIR}/${QUERY}_1.trimmed.fastq.gz ${OUTDIR}/${QUERY}_1.unpaired.fastq.gz ${OUTDIR}/${QUERY}_2.trimmed.fastq.gz ${OUTDIR}/${QUERY}_2.unpaired.fastq.gz ILLUMINACLIP:/N/u/tjlam/Quartz/bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:${MINLEN}
        echo "Trimming Complete: ${ERR}"
    else
        echo "FASTQ missing: ${ERR}"
        exit 1 
    fi

    # get number of reads
    trimmed_read_count=$(gzip -dc ${OUTDIR}/${QUERY}_1.trimmed.fastq.gz | wc -l)
    trimmed_read_count=$(expr "${trimmed_read_count}" / 4)

    # delete raw
    rm ${OUTDIR}/${QUERY}_1.fastq ${OUTDIR}/${QUERY}_2.fastq

    # map to hg38 and remove host
    if [ -f "${OUTDIR}/${QUERY}_1.trimmed.fastq.gz" ]; then
        echo "Running Bowtie - Removing Host: ${ERR}"
        bowtie2 -p ${THREADS} -x ${HG38} -1 ${OUTDIR}/${QUERY}_1.trimmed.fastq.gz -2 ${OUTDIR}/${QUERY}_2.trimmed.fastq.gz -U ${OUTDIR}/${QUERY}_1.unpaired.fastq.gz,${OUTDIR}/${QUERY}_2.unpaired.fastq.gz  --un-conc-gz ${OUTDIR}/${QUERY}_host_removed > ${OUTDIR}/${QUERY}_mapped_and_unmapped.sam 
        echo "Bowtie Complete"
    else
        echo "Trimmed missing: ${ERR}"
        exit 1
    fi

    # rename host removed
    mv ${OUTDIR}/${QUERY}_host_removed.1 ${OUTDIR}/${QUERY}_host_removed_R1.fastq.gz
    mv ${OUTDIR}/${QUERY}_host_removed.2 ${OUTDIR}/${QUERY}_host_removed_R2.fastq.gz

    # remove sam file
    rm ${OUTDIR}/${QUERY}_mapped_and_unmapped.sam

    # removed trimmed and unpaired
    rm ${OUTDIR}/${QUERY}_1.trimmed.fastq.gz ${OUTDIR}/${QUERY}_2.trimmed.fastq.gz ${OUTDIR}/${QUERY}_1.unpaired.fastq.gz ${OUTDIR}/${QUERY}_2.unpaired.fastq.gz

    # count host removed
    host_rm_count=$(gzip -dc ${OUTDIR}/${QUERY}_host_removed_R1.fastq.gz | wc -l)
    host_rm_count=$(expr "${host_rm_count}" / 4)
    
    # run kraken
    if [ -f "${OUTDIR}/${QUERY}_host_removed_R1.fastq.gz" ]; then
        echo "Running Kraken: ${ERR}"
        kraken2 --db ${KRAKENDB} --thread ${THREADS} --report ${OUTDIR}/${QUERY}.kreport --out ${OUTDIR}/${QUERY}.kraken --paired ${OUTDIR}/${QUERY}_host_removed_R1.fastq.gz ${OUTDIR}/${QUERY}_host_removed_R2.fastq.gz
        echo "Kraken Complete"
    else
        echo "host rm missing: ${ERR}"
        exit 1
    fi

    # remove raw kraken out to save space
    rm ${OUTDIR}/${QUERY}.kraken

    # remove PE
    rm ${OUTDIR}/${QUERY}_host_removed_R1.fastq.gz ${OUTDIR}/${QUERY}_host_removed_R2.fastq.gz

fi

# filter kraken results to bacteria only
echo "Filtering Kraken: ${ERR}"
python ${WD}/lib/kraken_filter.py -k ${OUTDIR}/${QUERY}.kreport -o ${OUTDIR}/${QUERY}.filtered.kreport
echo "Filtering Complete"

# run bracken on filtered kraken results
if [ -f "${OUTDIR}/${QUERY}.filtered.kreport" ]; then
    echo "Running Bracken: ${ERR}"
    bracken -d ${KRAKENDB} -i ${OUTDIR}/${QUERY}.filtered.kreport -o ${OUTDIR}/${QUERY}.s.bracken -w ${OUTDIR}/${QUERY}.s.breport -l S
    bracken -d ${KRAKENDB} -i ${OUTDIR}/${QUERY}.filtered.kreport -o ${OUTDIR}/${QUERY}.g.bracken -w ${OUTDIR}/${QUERY}.g.breport -l G
    echo "Bracken Complete"
else
    echo "missing Kraken: ${ERR}"
    exit 1
fi

# log metadata
echo -e "${QUERY}\t${ERR}\t${read_type}\t${raw_read_count}\t${trimmed_read_count}\t${host_rm_count}" >> $LOG
