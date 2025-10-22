#!/bin/bash
#SBATCH --job-name=vast2
#SBATCH --mail-type=END,FAIL          
#SBATCH --mail-user=joaosilvasantos@igc.gulbenkian.pt      
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=4             
#SBATCH --mem=50gb                     
#SBATCH --time=5-00:00:00               
#SBATCH --output=vast2%j.log   
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
processingFolder=$(cat $metadataFile | grep processingFolder | awk -F , '{print $2}')

resultsFolder="/ifs/igc/folders/RRA/zhang_files/vast/output_dir"

##############################################
############## Other variables ###############
##############################################

databaseDir="/ifs/igc/folders/RRA/casella_files/references"
combineOutput=$resultsFolder/combineOutput_INCLUSION_LEVELS.tsv
tidyOutput=$resultsFolder/tidyOutput.tsv
combineOutputFiltered=$resultsFolder/combineOutput_INCLUSION_LEVELS_tidy.tsv

###########################################################################################
######################################### Commands ########################################
###########################################################################################

## Align commands are run for each sample with previous batch file.
## Several files are created and saved in a folder named "to_combine"
## And now "combine" command is called to combine all those files into one.
## Note: "to_combine" folder is created automatically by vast-tools align command;
## and "combine" command will look for it automatically.

echo "********************************************************************"
echo "STARTING TO RUN BATCH FILE TO RUN VAST-TOOLS COMBINE and TIDY"
date
echo "********************************************************************"
echo ""

#############################################################
########################## COMBINE ##########################
#############################################################

cmd_combine="vast-tools combine \
        --output $resultsFolder \
        --sp hg38 \
        --dbDir $databaseDir \
        --cores 6"

echo "********************************************************************************"
echo "Running vast-tools combine with the following command:"
echo $cmd_combine
echo "********************************************************************************"
echo ""
singularity exec $image_path $cmd_combine
echo ""

## There is no argument to set output file name, so we have to move it manually
srun mv $resultsFolder/INCLUSION_*.tab $combineOutput

#############################################################
############################ TIDY ###########################
#############################################################

cmd_tidy="vast-tools tidy \
        $combineOutput \
        --min_Fr 0.5 \
        --min_SD 0 \
        -outFile $tidyOutput \
        --noVLOW \
        --p_IR"

echo "********************************************************************************"
echo "--min_Fr 0.5"
echo "This removes splicing events with low read support"
echo "********************************************************************************"

echo "--min_SD 0"
echo "This sets the minimum deviation to the integer next to it"
echo "********************************************************************************"

echo "--noVLOW"
echo "This removes the samples with very low coverage"
echo "********************************************************************************"

echo "--p_IR"
echo "Vast-tools applies a binomial test to evaluate the statistical significance of an IR event"
echo "This filters IR by the p-value of the binomial test"


echo "********************************************************************************"
echo "Running vast-tools tidy with the following command:"
echo $cmd_tidy
echo "********************************************************************************"
echo ""
srun apptainer exec $image_path $cmd_tidy
echo ""

#############################################################
################## FILTER COMBINE WITH TIDY #################
#############################################################

echo "********************************************************************************"
echo "Filtering combine table with tidy and saving in:"
echo "$combineOutputFiltered"
echo "********************************************************************************"
echo ""
srun awk -F'\t' 'NR==FNR{a[$1]; next} $2 in a' $tidyOutput $combineOutput > $combineOutputFiltered

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
#        ~/general/outFiles/$datasetName/slurm_$SLURM_JOB_ID.$SLURM_JOB_NAME.out

