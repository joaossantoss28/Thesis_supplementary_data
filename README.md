# Thesis_supplementary_data
This repository is composed of all files used for the production of the desired thesis

# Coding Languages
1.	Python 
Python is characterised as an interpreted, object-oriented programming language used to build data structures. It is one of the most popular programming languages, making it widely used in the scientific field 42. Python has an extensive ecosystem of libraries such as NumPy, Pandas and BioPython.
Pandas – Package that provides fast, flexible and expressive data structures such as data frames. It allows users to run certain tasks in data cleaning, transformation and exploration, supporting different file types, including .csv and Excel files43. 
Biopython – Collection of Python tools, designed specifically for bioinformatics and computational biology, with modules specific for reading and writing different sequence files as well as sequence alignments.44,45

2.	Bash (Shell Scripting)
Bash is a scripting language designed for Unix-based systems, for example, the LOBO HPC, which is based on this language. To navigate and work in most HPC clusters, such as LOBO, Bash is the primary scripting language used directly in the terminal and through .sh files. These scripts allow for automation and job scheduling, creating a fully functional framework to run specific tasks.46

3.	R
R is a programming language designed for statistical computing and an environment for graphics. It provides a wide variety of statistical and graphical techniques and good-quality plots 47. R can be extended via packages such as:
  •	ggplot2 – a system to create graphs with different details 48 .
  •	dplyr – a package suited for data manipulation, to help manage large datasets in a simple and efficient way 49.
  •	biomaRt – this package allows users to access Ensembl databases by connecting to the BioMart API, querying a wide range of data such as gene names, IDs, and coordinates - without manually downloading files. It is also possible to filter the results based on what information we want 50,51.
  •	readr – the goal of readr is to provide a fast way to read data from delimited files, such as CSV, TSV 52.


# Data Format
In computational biology, data is stored in standardised formats to ensure interoperability, consistency, and efficient processing. Each file format has a specific type of information, ranging from raw sequencing reads to files with processed genomic variants. 

1.	FASTQ
A FASTQ file is a text file that contains raw sequence data, obtained directly from the sequencing machines. For a single end run, only one FASTQ file is created for each sample. For a paired-end run, two FASTQ files are created for each sample (one R1 and one R2). Generally, FASTQ files are compressed and have the .fastq.gz extension.53
A FASTQ file consists of 4 lines:
1.	A sequence identifier, with the sequencing run label and the cluster information. 
2.	The sequence itself, with nucleotides (A, C, T, G and N)
3.	A separator line with the plus (+) sign
4.	The base quality scores, for the given sequence, using ASCII characters to represent numerical quality scores.


2.	FASTA
FASTA files are text files that are similar to FASTQ files, also storing nucleotide sequences, yet they do not include quality scores. Each FASTA file contains one or more entries, each of which has two parts:
1.	The identifier line, starting with a > symbol, the identifier, and some extra information, sometimes about gene function, location, species, length or coverage.
2.	The sequence line itself
The next line starts with a >, indicating that the following sequence ends and a new one starts.


3.	BAM
A BAM file is the compressed format of a SAM file, and contains information about read sequences, mapping quality, alignment position, etc. These files store aligned sequencing reads from sequencing techniques, and they are used for variant calling, visualisation, and downstream analysis, such as creating gene count data.55,56


4.	CSV and TSV
CSV (comma-separated values) files are one of the simplest and most used text file formats to store data. Each line/row has values separated by commas, and is highly compatible with Excel, and programming languages including Python and R.

5.	GTF File Format
A Gene Transfer Format file (GTF) is a tabular text format, each line of which corresponds to an annotation, or feature. Each line has 9 columns related to: reference sequence, source, method, start position, score, strand, phase, and group. This helps identify genes, transcripts, exons and coding sequences. 



# SRA TOOLKIT 

The subcommand fasterq-dump is the most recent and an upgraded version of the previous fastq-dump. While it requires more computational resources, it is faster and more efficient. The file sra.sh  is a template batch script designed to download sequencing data, using the fasterq-dump command, by taking a metadata file in CSV format with information about: datasetName, rawdataFolder, processingFolder, sampleNamesFile, fastqExtension, read1Extension, read2Extension. This structure enhances the reproducibility of the pipeline with different datasets, without directly changing the information in the script. The flag –outdir specifies the output directory for the fastq files, and the –threads sets the number of CPUs used. Increasing the number of threads reduces runtime but increases the computational load on the HPC cluster. 


# FASTQ
The file fastq.sh uses a metadata file similar to the one used in sra.sh. The output directory is specified with the flag –o and then the path of the output directory. The input files must be specified; if paired sequencing, both should be referenced.

After downloading the FASTQ files described in Deschênes, M.. et al.   , the first standard step was to perform a quality control assessment using FastQC (v0.12.1) and MultiQC (v1.20), as described in the Methods section. The results indicated overall good quality, as detailed in Figure 15. 
Despite the excellent quality, I performed a conservative trimming step using Trim-Galore (v0.6.8) 59. Default parameters were applied, including quality trimming with a Phred score cutoff of 20 and discarding reads shorter than 20 bp, with standard Illumina adapters specified for removal. The output report, multiqc_report.html, provides detailed metrics on GC content, base quality scores, number of reads and duplicates, and other relevant quality indicators


# STAR
When STAR attempts to align a read to the genome, it identifies  the longest exact match, the Maximal Mappable Prefix (MMP) that is present in the genome before encountering any mismatches). Once it finds this position, it stops the MMP there and defines whether this is a splice junction, a mismatch or a sequencing error. Considering this, STAR is very effective when detecting and analysing splicing problems, since it is aware of each junction present in the genome. 

The indexing step consists of using a specific set of flags. The –runMode specifies the STAR mode to generate the genome index. The –genomeDir to specify the directory of this process, where the index is going to be stored. The –genomeFastaFiles and the –sjdbGTFfile provide the paths of the reference genome, in this case, with the fastq and gtf format, respectively. Both reference files were retrieved from the Ensembl database.
After the genome index directory is completed, the final step is aligning the files to a reference genome. For this step, we used the script align.sh.
The parameter –runThreadN indicates the number of CPUs allocated for the analysis, while  –genomeDir defines the directory containing the previously generated genome index. The option –readFilesIn indicates the path to the input files (FASTQ files in this case), and the –readFilesCommand is used when the FASTQ files are compressed; here, the zcat command is applied to decompress the files, since STAR cannot process compressed inputs. The –sjdbGTFfile corresponds to the path of the reference GTF file, obtained from Ensembl 62 and –sjdbOverhang ensures that 100 bases of genomic sequence are added on each side of every known splice junction during genome index generation, to achieve a more accurate splice-aware alignment. Finally, the –outFileNamePrefix specifies the path where all output files will be stored, and –outSAMtype defines the format of the alignment files, which in this case are generated as BAM files sorted by genomic coordinates. These parameters collectively characterize the genomic index and alignment process. 
Alignment  

I performed read alignment with STAR solely for visualisation purposes. As described in the Methods section, the index and alignment steps with STAR were performed to generate four different BAM files, one for each biological sample. These BAM files were not used for downstream quantification but were exclusively intended for the visual inspection of splicing events using the Integrative Genomics Viewer (IGV).   





# VAST
The script align.sh  is used to execute the first step. In this script, the raw FASTQ files are provided as input “$input1 and $input2”. The –sp flag specifies the species assembly, while the –name defines the sample name. The --dbDir flag indicates the path to the database directory retrieved from VastDB, and --c sets the number of CPUs allocated for the task, and the –output flag specifies the directory where the output files will be stored.
After aligning the reads to the reference database, vast-tools combine merges the output files from each sample into a consolidated dataset. As shown in Figure 13, the combine step requires access to the results directory from the previous align step. The –output flag specifies this directory, while the remaining parameters are similar to those used in the align.sh script.


The final output file includes multiple columns, each representing distinct types of information. Each row corresponds to a specific alternative splicing event at a given genomic position. A single gene may have multiple entries if it exhibits more than one splicing event. The file format is as follows:
•	Gene – The gene symbol name (e.g SYDE1)
•	EVENT – Set of letters and numbers identifying that given event, with information about the species, the type of alternative splicing event and the ID of the event (e.g HsaINT0162432)
•	COORD – The genomic coordinate of the alternative sequence (e.g chr19:15109398-15109703, location of the retained intron for this event)
•	LENGTH – The length of the alternative sequence
•	FullCO – Full set of genomic coordinates of the alternative splicing event (e.g chr19:15109056-15109397=15109704-15110348:+, here we have the information of an IR event that occurred in chromosome 19, the coordinates of the upstream exon being 15109056-15109397, the coordinates of the downstream exon 15109704-15110348, and the strand direction in this case +)
•	COMPLEX – Type of alternative splicing event (e.g S for exon skipping, IR for intron retention, Alt3 for alternative 3’ splice and Alt5 for alternative 5’ splice)
In addition, the output includes a dedicated column for each sample, labelled with the corresponding sample ID, containing the PSI values for each splicing event. For every sample, there is also an accompanying column that provides a set of quality scores computed by vast-tools. These scores are accompanied by the number of reads supporting the inclusion of the splicing event, as well as the number of reads that support its exclusion. 


# CLUSTERPROFILER
To visualize the results, I generated a dot plot highlighting the top enriched biological processes using the clusterProfiler R package.

