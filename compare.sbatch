#!/bin/bash
#SBATCH --job-name=vast3
#SBATCH --mail-type=END,FAIL          
#SBATCH --mail-user=joaosilvasantos@igc.gulbenkian.pt      
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=4             
#SBATCH --mem=50gb                     
#SBATCH --time=5-00:00:00               
#SBATCH --output=vast3%j.log   
# IMPORTANT: if you need to separate error output, add a -e or --error line, specifying
# the filename for the error log (use the --output line as a template for the new line).

export SINGULARITY_BIND="/ifs"  # required for access to /ifs (to improve safety, you can restrict to specific folders)


image_path="/Users/.singularity_containers/ADA/vast.sif"

###########################################################################################
##################################### Input Variables #####################################
###########################################################################################

## NOTE!!!! This batch takes 2 arguments as input:
## 1. Metadata file (typical, as in other batches)
## 2. Vast-tools metadata file (that will indicate the groups to be compared)
## Results will be saved in the same folder as this file.

metadataFile=$1
if [[ ! -f $metadataFile ]]; then
	echo "Metadata file not found!"
        exit 1
fi

metadata_vast_tools=$(readlink -f "$2")
if [[ ! -f $metadata_vast_tools ]]; then
	echo "This file takes another input: Metadata vast-tools (not found!). This file should indicate the groups to be compared."
        echo "And the results will be saved in the same folder as this file."
        exit 1
fi

# Example of metadata content:

# sampleName,group,groupName
# ID109M18,a,Mock18
# ID65M18,a,Mock18
# ID76M18,a,Mock18
# ID109A18,b,Aza
# ID65A18,b,Aza
# ID76A18,b,Aza


datasetName=$(cat $metadataFile | grep datasetName | awk -F , '{print $2}')
processingFolder=$(cat $metadataFile | grep processingFolder | awk -F , '{print $2}')
pairedSamples=$(cat $metadataFile | grep pairedSamples | awk -F , '{print $2}')
pairedSamples=${pairedSamples^^}
if [[ -z $pairedSamples ]]; then
        echo "Define pairedSamples variable in metadata file! (TRUE or FALSE)"
        exit 1
fi

inputFolder=${processingFolder}
resultsFolder=${processingFolder}
resultsFolderName=$(basename $resultsFolder)

echo $resultsFolder

#############################################################
########### Input and output file names and paths ###########
#############################################################

databaseDir="/ifs/igc/folders/RRA/references"
combineOutputOriginal=$resultsFolder/combineOutput_INCLUSION_LEVELS.tsv
tidyOutput=$resultsFolder/tidyOutput.tsv
combineOutputFiltered=$resultsFolder/combineOutput_INCLUSION_LEVELS_tidy.tsv
## CHOOSE COMBINE OUTPUT (FILTERED OR ORIGINAL):
combineOutput=$combineOutputOriginal

compareOutput_toPlot_fileName=compareOutput_toPlot.tsv
compareOutput_withDPSI_fileName=compareOutput_withDPSI.tsv
compareOutput_toPlot=$resultsFolder/$compareOutput_toPlot_fileName
compareOutput_withDPSI=$resultsFolder/$compareOutput_withDPSI_fileName

#############################################################
################### Groups to be compared ###################
#############################################################

## Go to metadata_vast_tools file and get samples from group A and group B:
## Get column 1 (sample name) for the lines that contain "a" and then a "b" in column 2:
group_a_samples=$(cat $metadata_vast_tools | grep ",a," | awk -F , '{print $1}')
group_b_samples=$(cat $metadata_vast_tools | grep ",b," | awk -F , '{print $1}')
## Replace spaces by commas in the sample lists:
group_a_samples=$(echo $group_a_samples | sed 's/ /,/g')
group_b_samples=$(echo $group_b_samples | sed 's/ /,/g')
nameGroupA=$(cat $metadata_vast_tools | grep -m 1 ",a," | awk -F , '{print $3}')
nameGroupB=$(cat $metadata_vast_tools | grep -m 1 ",b," | awk -F , '{print $3}')


###########################################################################################
######################################### Commands ########################################
###########################################################################################

## After align, combine and tidy, we have a file with PSI values for each sample
## (combineOutput_INCLUSION_LEVELS.tsv).
## Now it is possible to compare samples from two different groups.
## There are 2 options with vast-tools:
## 1. compare (Delta PSI = average group B - average group A)
## 2. diff (Bayesian Inference to compare PSI values between groups)

echo "********************************************************************"
echo "STARTING TO RUN VAST-TOOLS COMPARE AND DIFF BATCH FILE"
echo $resultsFolderName
date
echo "********************************************************************"

echo ""
echo "Group A samples:"
echo $group_a_samples
echo "Group B samples:"
echo $group_b_samples
echo ""

#############################################################
########################## COMPARE ##########################
#############################################################

cmd_compare="vast-tools compare \
        $combineOutput \
        -a $group_a_samples \
        -b $group_b_samples \
        -name_A $nameGroupA \
        -name_B $nameGroupB \
        --min_dPSI 0 \
        --min_range 0 \
        --noVLOW \
        --noB3 \
        --p_IR \
        --GO \
        --sp hg38 \
        --dbDir $databaseDir"

if [ "$pairedSamples" == "TRUE" ]; then
        cmd_compare="${cmd_compare} --paired"
fi

cmd_compare_to_plot="${cmd_compare} --outFile $compareOutput_toPlot_fileName"
cmd_compare_with_dPSI="${cmd_compare} --print_dPSI --outFile $compareOutput_withDPSI_fileName"


echo "********************************************************************************"
echo "Running vast-tools compare with the following command: (generate file to plot)"
echo ""
echo $cmd_compare_to_plot
echo "********************************************************************************"
echo ""
singularity exec $image_path $cmd_compare_to_plot
echo ""
echo ">> Moving output files to results folder..."
## Compare plot doesn't allow to choose output directory (only base name), so we have to move the file:
mv $inputFolder/$compareOutput_toPlot_fileName /ifs/igc/folders/RRA/zhang_files/data_quiescence_deep
echo ""


echo "********************************************************************************"
echo "Running vast-tools compare with the following command: (generate file with DPSI)"
echo ""
echo $cmd_compare_with_dPSI
echo "********************************************************************************"
echo ""
singularity exec $image_path $cmd_compare_with_dPSI
echo ""
echo ">> Moving output files to results folder..."
mv $inputFolder/$compareOutput_withDPSI_fileName /ifs/igc/folders/RRA/zhang_files/data_quiescence_deep
echo ""


echo "*****************************************************************"
echo ">> Keeping only columns of samples that are being compared..."
echo "*****************************************************************"
echo ""

file=$compareOutput_toPlot
all_samples="$group_a_samples,$group_b_samples"
columnsToPrint="1-6"
IFS=',' read -ra samples_array <<< "$all_samples"
for sample in "${samples_array[@]}"; do
        sampleColumn1=$(awk -v RS='\t' "/$sample/"'{print NR; exit}' $file)
        sampleColumn2=$(awk -v RS='\t' "/$sample-Q/"'{print NR; exit}' $file)
        columnsToPrint="${columnsToPrint},$sampleColumn1,$sampleColumn2"
done

## Extract columnsToPrint from compareOutput_toPlot:
cut -f "$columnsToPrint" "$compareOutput_toPlot" > $resultsFolder/temp.tsv
mv $resultsFolder/temp.tsv $compareOutput_toPlot

## Extract columnsToPrint from compareOutput_withDPSI:
numberColumns=$(awk '{print NF; exit}' $compareOutput_withDPSI)
columnsToPrintDPSI="${columnsToPrint},$numberColumns"
cut -f "$columnsToPrintDPSI" "$compareOutput_withDPSI" > $resultsFolder/temp.tsv
mv $resultsFolder/temp.tsv $compareOutput_withDPSI

#############################################################
#################### PLOT COMPARE RESULTS ###################
#############################################################

cmd_plot="vast-tools plot \
        -v TRUE -m 1000 \
        $compareOutput_toPlot"

if [[ -f $resultsFolder/plots_config.tsv ]]; then
        cmd_plot="$cmd_plot --config $resultsFolder/plots_config.tsv --groupMeans TRUE"
fi

echo "********************************************************************************"
echo "Running vast-tools plot with the following command:"
echo ""
echo $cmd_plot
echo "********************************************************************************"
echo ""
srun singularity exec $image_path $cmd_plot
echo ""

###########################################################################################
#################################### Move output file  ####################################
###########################################################################################

#mkdir -p ~/general/outFiles/$datasetName/
#srun mv ~/general/outFiles/slurm_$SLURM_JOB_ID.$SLURM_JOB_NAME.out \
#        ~/general/outFiles/$datasetName/slurm_$SLURM_JOB_ID.$SLURM_JOB_NAME.$resultsFolderName.out

