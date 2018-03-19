#!/bin/bash

#---
# Title: 'First real multiplexing assay with ONT, basecall and demux'
# Author: 'antoine.bridier-nahmias@inserm.fr'
# Date: '06-03-2018'
#---

if [ "$#" == 0 ]; then
	echo -e "--otoro:\t\tBasecalling and demux with albacore"
	echo -e "--petunia_pig:\t\tPorechop"
	echo -e "--ryan_the_racoon:\tSubsampling and quality filtering"
	echo -e "--howard_the_duck:\tIllumina read quality filtering with BBduk"
	echo -e "--monocycle:\t\tUnicycler"
fi

### PATH
spades_path=${HOME}/Documents/tool/SPAdes-3.11.1-Linux/bin
racon_path=${HOME}/Documents/tool/SPAdes-3.11.1-Linux/racon/build/bin
pilon_path=${HOME}/Documents/tool/pilon
bowtie2_path=${HOME}/Documents/tool/bowtie2-2.3.4.1-linux-x86_64
samtools_path=${HOME}/Documents/tool/samtools-1.7
export PATH=${PATH}:${spades_path}:${racon_path}:${pilon_path}:${bowtie2_path}:${samtools_path}
JAVA_HOME="/usr/lib/jvm/java-9-openjdk-amd64/bin/"
export JAVA_HOME

### Global vars
n_threadus=10
wanted_barcoda="08"

### Tools
porechop="$HOME/Documents/tool/Porechop/porechop-runner.py"
fastq_to_fastq="$HOME/Documents/tool/Fast5-to-Fastq/fastq_to_fastq.py"
unicycler="$HOME/Documents/tool/Unicycler/unicycler-runner.py"
bbduk="$HOME/Documents/tool/bbmap/bbduk.sh"
MinionQC="${HOME}/Documents/tool/minion_qc/MinionQC.R"
###############################################################################
###############    Basecalling and demux with albacore    #####################
###############################################################################
basecall_dir="./data/basecall"
fast5_dir="/home/nanopore/nanopore_data/20180227_1652_Multiplex_first_real_27022018/fast5/"

if [[ "$@" =~ "--otoro" ]]; then
	if [ ! -d ${basecall_dir} ]; then mkdir ${basecall_dir}; fi
	echo "Basecall and Demux"
	echo "starting: $(date)"
	time \
#	read_fast5_basecaller.py \
#		--input ${fast5_dir} \
#		--output_format "fastq" \
#		--save_path ${basecall_dir} \
#		--flowcell FLO-MIN106 \
#		--kit  SQK-LSK108 \
#		--recursive \
#		--barcoding \
#		--worker_threads ${n_threadus}
	Rscript ${MinionQC} \
		-i "./data/basecall/sequencing_summary.txt" \
		-o "./data/basecall/minionQC" \
		-p ${n_threadus}
	echo "Basecall and Demux"
	echo "ending: $(date)"

fi

###############################################################################
##################     Trimming and re-de-muxing     ##########################
###############################################################################
if [[ "$@" =~ "--petunia_pig" ]]; then
	for n_barcodus in ${wanted_barcoda}; do
		barcodus=barcode${n_barcodus}
		echo "Treating sample: ${barcodus}"
		fastq_dir="./data/fastq/${barcodus}"
		if [ -d ${fastq_dir} ]; then
			echo "Trimming and re-de-muxing"
			echo "starting $(date)"			
			time \
			cat ${fastq_dir}/*_*.fastq | gzip > ${fastq_dir}/${barcodus}.fastq.gz
			trimmed_dir="${fastq_dir}/trimmed/"
			time \
			$porechop \
				-i ${fastq_dir}/${barcodus}.fastq.gz \
				-b ${trimmed_dir} \
				--discard_middle \
				--threads ${n_threadus}
			echo "Trimming and re-de-muxing"
			echo "ending: $(date)"			
#			echo "La limpieza"
#			find ${trimmed_dir} -type f ! -name "*${n_barcodus}*" | xargs rm	
		fi
	done
fi

###############################################################################
##################    Subsampling and quality filtering      ##################
###############################################################################
if [[ "$@" =~ "--ryan_the_racoon" ]]; then
	echo "Subsampling and quality filtering"
	for n_barcodus in ${wanted_barcoda}; do
		barcodus=barcode${n_barcodus}
		curr_fastq="./data/fastq/barcode${n_barcodus}/trimmed/BC${n_barcodus}.fastq.gz"
		subsampled_dir="$(dirname ${curr_fastq})/subsampled"
		min_lengthus=2000
		target_basus=500000000
		if [ -s "${curr_fastq}" ]; then
			echo "Treating sample: ${barcodus}"
			echo "starting $(date)"
			if  [ ! -d ${subsampled_dir} ]; then mkdir ${subsampled_dir}; fi
			time \
			$fastq_to_fastq \
				--min_length ${min_lengthus} \
				--target_bases ${target_basus}\
				${curr_fastq} | gzip > "${subsampled_dir}/${barcodus}.fastq.gz"
			echo "ending $(date)"
		fi
	done
fi

###############################################################################
###################    Illumina read quality filtering      ###################
###############################################################################
if [[ "$@" =~ "--howard_the_duck" ]]; then
	echo "Illumina read quality filtering with BBduk"
	echo "starting $(date)"
	for n_barcodus in ${wanted_barcoda}; do
		barcodus=barcode${n_barcodus}
		echo "$barcodus"
		illumina_dir="./data/fastq/barcode${n_barcodus}/illumina"
		illumina_clean_dir="${illumina_dir}/clean"
		n_illumina=$(find ${illumina_dir} -type f -name "*fastq*"| wc -l)
		if [ ! -d ${illumina_clean_dir} ]; then mkdir ${illumina_clean_dir}; fi
		if [ "$n_illumina" == 1 ]; then
			echo "$n_illumina"
			illumina_in=$(find ${illumina_dir} -type f -name "*fastq*")
			${bbduk} \
					in=${illumina_in} \
					out=${illumina_clean_dir}/${barcodus}_clean.fastq \
					maq=30 \
					minlen=50 
		elif [ "$n_illumina" == 2 ]; then
			echo "$n_illumina"
			illumina_R1_in=$(find ${illumina_dir} -type f -name "*R1*fastq*")
			illumina_R2_in=$(find ${illumina_dir} -type f -name "*R2*fastq*")
			${bbduk} \
					in1=${illumina_R1_in} \
					in2=${illumina_R2_in} \
					out1=${illumina_clean_dir}/${barcodus}_R1_clean.fastq \
					out2=${illumina_clean_dir}/${barcodus}_R2_clean.fastq \
					maq=30 \
					minlen=50 
		else
			echo "Wrong number of illumina reads: $n_illumina"
			echo "Skipping ${barcodus}"
		fi

	done
	echo "ending $(date)"
fi

###############################################################################
##################    Assembling ONT and Illumina reads      ##################
###############################################################################
if [[ "$@" =~ "--monocycle" ]]; then
	echo "Assembling ONT and Illumina reads"
	for n_barcodus in ${wanted_barcoda}; do
		barcodus=barcode${n_barcodus}
		curr_fastq="./data/fastq/barcode${n_barcodus}/trimmed/subsampled/${barcodus}.fastq.gz"
		subsample_dir="$(dirname ${curr_fastq})/subsampled/${barcodus}"
		illumina_dir="./data/fastq/${barcodus}/illumina"
		assembled_dir="./data/assembled/${barcodus}"
		if [ -s "${curr_fastq}" ]; then
			echo "Treating sample: ${barcodus}"	
			echo "starting $(date)"
			if  [ ! -d ${assembled_dir} ]; then mkdir -p ${assembled_dir}; fi
			n_illumina_reads=$(ls -1 ${illumina_dir}/*.fastq* | wc -l)
			if [ "$n_illumina_reads" == 1 ]; then 
				illumina_fastq=`find ${illumina_dir} -type f -name "*clean*fastq*"`
				unicycler_input="--unpaired ${illumina_fastq}"
			fi
			if [ "$n_illumina_reads" == 2 ]; then 
				illumina_R1=$(find ${illumina_dir} -type f -name "*R1*clean*.fastq*")
				illumina_R2=`find ${illumina_dir} -type f -name "*R2*clean*.fastq*"`
				unicycler_input="-1 ${illumina_R1} -2 ${illumina_R2}"
			fi
			echo -e "Launching unicycler with reads: ${curr_fastq} and ${unicycler_input}"
			${unicycler} \
				${unicycler_input} \
				--long ${curr_fastq}\
				--verbosity 1 \
				--out ${assembled_dir} \
				--keep 0 \
				--threads ${n_threadus}
			echo "ending $(date)"
		fi
	done
fi













###############################################################################
###################                                       #####################
###############################################################################

###############################################################################
###################                                       #####################
###############################################################################
