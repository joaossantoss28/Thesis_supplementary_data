import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def load_genome(fasta_file):
    genome = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome[record.id] = record.seq
    return genome

def extract_sequence_from_genome(genome, chrom, start, end):
    if chrom in genome:
        return genome[chrom][start-1:end]  # Extract sequence from the genome
    else:
        raise KeyError(f"Chromosome {chrom} not found in the genome")


# Function to calculate GC content
def calculate_gc_content(sequence):
    return gc_fraction(sequence) * 100 if sequence else None

def add_gc_ratio_column(intron_data, genome):
    gc_ratios = []

    for index, row in intron_data.iterrows(): 
        chrom = str(row['CHR'])
        chrom = chrom.replace("chr", "")        # Normalize the chromosome format by removing 'chr' prefix if present
        start = int(row['Intron_Start'])
        end = int(row['Intron_End'])

        # Extract sequence
        sequence = extract_sequence_from_genome(genome, chrom, start, end)
        # Calculate GC content
        gc_content = calculate_gc_content(sequence)
        gc_ratios.append(gc_content)

    # Add the GC ratio column to the DataFrame
    intron_data['GC_Ratio'] = gc_ratios
    return intron_data


# File paths
csv_file_path = "/ifs/igc/folders/RRA/zhang_files/data_quies_senesc/data7.csv"
reference_file = "/ifs/igc/folders/RRA/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" 

# Read intron data from CSV
intron_data = pd.read_csv(csv_file_path)

print(intron_data.head())

genome = load_genome(reference_file)

# Add the GC_Ratio column
intron_data = add_gc_ratio_column(intron_data, genome)
# Save the updated DataFrame
output_file = "/ifs/igc/folders/RRA/zhang_files/data_quies_senesc/data7_withgc.csv"
intron_data.to_csv(output_file, index=False)


