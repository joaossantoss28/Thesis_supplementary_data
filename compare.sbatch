#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=170G
#SBATCH --job-name="vast_tools_3_combine"
#SBATCH --output=vast3_%j.vast_tools_3_align.out

image_path="/mnt/beegfs/apptainer/images/vast-tools_latest.sif"

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

databaseDir="/mnt/nfs/lobo/IMM-NFS/imm/files/vast/database_dir"
combineOutputOriginal=$resultsFolder/combineOutput_INCLUSION_LEVELS.tsv
tidyOutput=$resultsFolder/tidyOutput.tsv
combineOutputFiltered=$resultsFolder/combineOutput_INCLUSION_LEVELS_tidy.tsv
## CHOOSE COMBINE OUTPUT (FILTERED OR ORIGINAL):
combineOutput=$combineOutputFiltered

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
        --min_dPSI 10 \
        --min_range 5 \
        --noVLOW \
        --noB3 \
        --p_IR \
        --GO \
        --sp mm10 \
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
srun apptainer exec $image_path $cmd_compare_to_plot
echo ""
echo ">> Moving output files to results folder..."
## Compare plot doesn't allow to choose output directory (only base name), so we have to move the file:
srun mv $inputFolder/$compareOutput_toPlot_fileName $resultsFolder
echo ""


echo "********************************************************************************"
echo "Running vast-tools compare with the following command: (generate file with DPSI)"
echo ""
echo $cmd_compare_with_dPSI
echo "********************************************************************************"
echo ""
srun apptainer exec $image_path $cmd_compare_with_dPSI
echo ""
echo ">> Moving output files to results folder..."
srun mv $inputFolder/$compareOutput_withDPSI_fileName $resultsFolder
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
srun apptainer exec $image_path $cmd_plot
echo ""

#############################
######### GO files ##########
#############################

echo "*****************************************************************"
echo ">> Reorganizing GO files from compare output..."
echo "*****************************************************************"
echo ""

## flag --GO will create files for posterior analysis.
## Put all those files in a folder (to organize)
mkdir -p $resultsFolder/GO_files
## Remove files from previous runs:
rm -rf $resultsFolder/GO_files/*
## Move files to GO_files folder:
mv $inputFolder/All_Ev*.txt \
$inputFolder/AltEx*.txt \
$inputFolder/BG*.txt \
$inputFolder/IR_*.txt \
$resultsFolder/GO_files

## Remove duplicate files:
rm $resultsFolder/GO_files/*toPlot.txt

## Rename files (to make it easier to understand):
mv $resultsFolder/GO_files/All* $resultsFolder/GO_files/positive_geneIDs_to_GO.txt
mv $resultsFolder/GO_files/BG* $resultsFolder/GO_files/background_geneIDs_to_GO.txt
mv $resultsFolder/GO_files/AltEx* $resultsFolder/GO_files/alternative_exon_geneIDs_to_GO.txt
mv $resultsFolder/GO_files/IR_DOWN* $resultsFolder/GO_files/intron_higher_retention_in_A_geneIDs_to_GO.txt
mv $resultsFolder/GO_files/IR_UP* $resultsFolder/GO_files/intron_higher_retention_in_B_geneIDs_to_GO.txt

## Create new file with negative genes (all genes from background minus all genes from positive)
fgrep -vwf $resultsFolder/GO_files/positive_geneIDs_to_GO.txt $resultsFolder/GO_files/background_geneIDs_to_GO.txt > $resultsFolder/GO_files/negative_geneIDs_to_GO.txt

#############################################################
########################### DIFF ############################
#############################################################

## Run diff with default parameters: (0.95 and 0.3)
diffOutputBaseName1=diffOutput_95_03
diffOutput1=$resultsFolder/${diffOutputBaseName1}.tsv
cmd_diff_1="vast-tools diff \
        -a $group_a_samples \
        -b $group_b_samples \
        --sampleNameA $nameGroupA \
        --sampleNameB $nameGroupB \
        --input $combineOutput \
        --prob 0.95 \
        --minDiff 0.1 \
        --minSamples 2 \
        --minReads 10 \
        -c 10 \
        --output $resultsFolder \
        -d $diffOutputBaseName1"

# if [ "$pairedSamples" == "TRUE" ]; then
#         cmd_diff_1="${cmd_diff_1} --paired=TRUE"
# fi

## Run diff with less strict parameters: (0.90 and 0.2)
diffOutputBaseName2=diffOutput_90_02
diffOutput2=$resultsFolder/${diffOutputBaseName2}.tsv
cmd_diff_2="vast-tools diff \
        -a $group_a_samples \
        -b $group_b_samples \
        --sampleNameA $nameGroupA \
        --sampleNameB $nameGroupB \
        --input $combineOutput \
        --prob 0.90 \
        --minDiff 0.2 \
        --minSamples 2 \
        --minReads 10 \
        -c 10 \
        --output $resultsFolder \
        -d $diffOutputBaseName2"

# if [ "$pairedSamples" == "TRUE" ]; then
#         cmd_diff_2="${cmd_diff_2} --paired=TRUE"
# fi

echo "********************************************************************************"
echo "Running vast-tools diff with the following command: (default parameters: 0.95 and 0.3)"
echo ""
echo $cmd_diff_1
echo "********************************************************************************"
echo ""
srun apptainer exec $image_path $cmd_diff_1
srun mv $resultsFolder/${diffOutputBaseName1}.tab $diffOutput1

echo ""
echo "********************************************************************************"
echo "Running vast-tools diff with the following command: (less strict parameters: 0.90 and 0.2)"
echo ""
echo $cmd_diff_2
echo "********************************************************************************"
echo ""
srun apptainer exec $image_path $cmd_diff_2
srun mv $resultsFolder/${diffOutputBaseName2}.tab $diffOutput2

echo ""
echo "*****************************************************************"
echo ">> Creating files with significant events from diffOutput files..."
echo "*****************************************************************"
echo ""

## Selects rows from diffOutput where the value in the sixth column is greater than or equal to 0.2 / 0.3:
## In short, select significant events:
awk '$6 >= 0.2' $diffOutput1 > $resultsFolder/${diffOutputBaseName1}_filtered_02.tsv
awk '$6 >= 0.3' $diffOutput1 > $resultsFolder/${diffOutputBaseName1}_filtered_03.tsv

## From vast-tools github about sixth column (or MV column):
## the minimum value (MV) for |ΔPSI| at 0.95 is 0.3, meaning that there is a 0.95
## probability that |ΔPSI| is greater than 0.3. Use this value to filter for events
## that are statistically likely to have at least a minimal difference of some magnitude
## that you deem to be biologically relevant.

#############################################################
##################### PLOT DIFF RESULTS #####################
#############################################################

## Extract events to plot:
awk 'NR==FNR{events[$2]; next} $2 in events' $resultsFolder/${diffOutputBaseName1}_filtered_02.tsv $combineOutput > $resultsFolder/${diffOutputBaseName1}_toPlot.tsv
## Extract columnsToPrint:
cut -f "$columnsToPrint" "$resultsFolder/${diffOutputBaseName1}_toPlot.tsv" > $resultsFolder/temp.tsv
mv $resultsFolder/temp.tsv $resultsFolder/${diffOutputBaseName1}_toPlot.tsv

cmd_plot="vast-tools plot \
        -v TRUE -m 1000 \
        $resultsFolder/${diffOutputBaseName1}_toPlot.tsv"

if [[ -f $resultsFolder/plots_config.tsv ]]; then
        cmd_plot="$cmd_plot --config $resultsFolder/plots_config.tsv --groupMeans TRUE"
fi

echo "********************************************************************************"
echo "Running vast-tools plot with the following command:"
echo ""
echo $cmd_plot
echo "********************************************************************************"
echo ""
srun apptainer exec $image_path $cmd_plot
echo ""


echo ""
echo "********************************************************************"
echo "FINISHED SCRIPT"
date
echo "********************************************************************"

###########################################################################################
#################################### Move output file  ####################################
###########################################################################################

#mkdir -p ~/general/outFiles/$datasetName/
#srun mv ~/general/outFiles/slurm_$SLURM_JOB_ID.$SLURM_JOB_NAME.out \
#        ~/general/outFiles/$datasetName/slurm_$SLURM_JOB_ID.$SLURM_JOB_NAME.$resultsFolderName.out

