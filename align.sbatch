#!/bin/bash
#SBATCH --job-name=vast_zh 
#SBATCH --mail-type=END,FAIL          
#SBATCH --mail-user=joaosilvasantos@igc.gulbenkian.pt      
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=4             
#SBATCH --mem=50gb                     
#SBATCH --time=5-00:00:00               
#SBATCH --output=vast_zh%j.log   
#SBATCH --array=0-8%8
# IMPORTANT: if you need to separate error output, add a -e or --error line, specifying
# the filename for the error log (use the --output line as a template for the new line).

export SINGULARITY_BIND="/ifs"  # required for access to /ifs (to improve safety, you can restrict to specific folders)


image_path="/Users/.singularity_containers/ADA/vast.sif"

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

input1=$rawDataFolder/${sampleName}${read1Extension}.$fastqExtension
input2=$rawDataFolder/${sampleName}${read2Extension}.$fastqExtension
input1=$(readlink -f $input1)
input2=$(readlink -f $input2)

resultsFolder="/ifs/igc/folders/RRA/zhang_files/vast/output_dir"
mkdir -p $resultsFolder

echo "Sample Name: $sampleName"
echo "Input 1 Path: $input1"
echo "Input 2 Path: $input2"
echo "Read 1 Extension: $read1Extension"
echo "Read 2 Extension: $read2Extension"

##############################################
############## Other variables ###############
##############################################

databaseDir="/ifs/igc/folders/RRA/casella_files/references"

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

srun singularity exec $image_path vast-tools align $input1 $input2 \
        --sp hg38 \
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


