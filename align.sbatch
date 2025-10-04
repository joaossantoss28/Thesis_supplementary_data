#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=170G
#SBATCH --job-name="vast_tools_1_align"
#SBATCH --output=vast_%j.vast_tools_1_align.out
#SBATCH --array=0-7%6

image_path="/mnt/beegfs/apptainer/images/vast-tools_latest.sif"

###########################################################################################
##################################### Input Variables #####################################
###########################################################################################

metadataFile=$1

if [[ ! -f $metadataFile ]]; then
	echo "Metadata file not found!"
        exit 1
fi

datasetName=$(cat $metadataFile | grep datasetName | awk -F , '{print $2}')
rawDataFolder=$(cat $metadataFile | grep rawDataFolder | awk -F , '{print $2}')
processingFolder=$(cat $metadataFile | grep processingFolder | awk -F , '{print $2}')
sampleNamesFile=$(cat $metadataFile | grep sampleNamesFile | awk -F , '{print $2}') 
fastqExtension=$(cat $metadataFile | grep fastqExtension | awk -F , '{print $2}')
read1Extension=$(cat $metadataFile | grep read1Extension | awk -F , '{print $2}')
read2Extension=$(cat $metadataFile | grep read2Extension | awk -F , '{print $2}')

echo "test:$sampleNamesFile"

readarray -t sampleNamesList < <(cat $sampleNamesFile)
sampleName=${sampleNamesList[$SLURM_ARRAY_TASK_ID]}

input1=$rawDataFolder/${sampleName}${read1Extension}.$fastqExtension.gz
input2=$rawDataFolder/${sampleName}${read2Extension}.$fastqExtension.gz
input1=$(readlink -f $input1)
input2=$(readlink -f $input2)

resultsFolder="/mnt/nfs/lobo/IMM-NFS/imm/files/vast/output_dir"
mkdir -p $resultsFolder

echo "Sample Name: $sampleName"
echo "Input 1 Path: $input1"
echo "Input 2 Path: $input2"
echo "Read 1 Extension: $read1Extension"
echo "Read 2 Extension: $read2Extension"

##############################################
############## Other variables ###############
##############################################

databaseDir="/mnt/beegfs/scratch/IMM/joaossantos/database_dir"

###########################################################################################
######################################### Commands ########################################
###########################################################################################

echo "********************************************************************"
echo "STARTING TO RUN BATCH FILE TO RUN VAST-TOOLS ALIGN"
date
echo "********************************************************************"
echo ""

echo "********************************************************************************"
echo "Running vast-tools align for sample $sampleName"
echo ""
echo "Inputs:"
echo "$input1"
echo "$input2"
echo ""
echo "Outputs:"
echo "$resultsFolder"
echo "********************************************************************************"
echo ""

srun apptainer exec $image_path vast-tools align $input1 $input2 \
        --sp mm10 \
        --name $sampleName \
        --dbDir $databaseDir \
        -c 8 \
        --output $resultsFolder

echo ""
echo "********************************************************************"
echo "FINISHED SCRIPT"
date
echo "********************************************************************"


###########################################################################################
#################################### Move output file  ####################################
###########################################################################################

#mkdir -p /mnt/beegfs/scratch/IMM/joaossantos/output_dir$datasetName/
#srun mv /mnt/beegfs/scratch/IMM/joaossantos/output_dir_$SLURM_JOB_ID.$SLURM_JOB_NAME.out \
#        ~/general/outFiles/$datasetName/slurm_$SLURM_JOB_ID.$SLURM_JOB_NAME.${sampleName}.out


