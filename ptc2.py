from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

STOP_CODONS = ["TAA", "TGA", "TAG"]

def extract_sequence_from_genome(genome, chrom, start, end):
    if chrom in genome:
        return genome[chrom][start-1:end]  # Extract sequence from the genome
    else:
        raise KeyError(f"Chromosome {chrom} not found in the genome")


def check_frameshift(exon_phase, intron_length):
    """
    Determines if intron retention causes a frameshift.

    Parameters:
        exon_phase (int): The phase of the exon (0, 1, or 2).
        intron_length (int): The length of the retained intron in nucleotides.

    Returns:
        bool: True if a frameshift occurs, False otherwise.
        int: New phase after retention.
    """

    new_phase = (exon_phase + intron_length) % 3 

    # % symbol is the remainder so 5%3=2, 3%3=0


    frameshift = new_phase != exon_phase

    return frameshift, new_phase



def check_for_ptc_in_sequence(dna_sequence, phase, strand):
    """
    Check for premature stop codons in the sequence based on the translation phase.
    
    Parameters:
        dna_sequence (str): The DNA sequence to check for PTC.
        phase (int): The phase of translation (0, 1, or 2).
        strand (str): The strand ("+" or "-").
    
    Returns:
        bool: True if a PTC is found, otherwise False.
        int: The position of the PTC (relative to the translated sequence).
    """
    dna_sequence = str(dna_sequence).upper()

    # Adjust reading frame based on the phase
    if phase == 0:
        start_pos = 0
    elif phase == 1:
        start_pos = 1
    elif phase == 2:
        start_pos = 2
    else:
        raise ValueError("Invalid phase. Phase must be 0, 1, or 2.")


    # Iterate through the sequence in 3-base codons
    for i in range(start_pos, len(dna_sequence) - 2, 3):  
        codon = dna_sequence[i:i+3]
        #print(f"Codon: {codon} at position {i}")
        if codon in STOP_CODONS:
            #print(f"Found stop codon {codon} at position {i}")
            return True, i  # Return position of PTC
    
    return False, None  # No stop codon found

def get_exon_phase_from_gtf(gtf_file, chrom, exon_start, exon_end):
    phase = None

    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            columns = line.strip().split("\t")

            gtf_chrom = columns[0]
            feature = columns[2]
            gtf_start = int(columns[3])
            gtf_end = int(columns[4])
            gtf_phase = columns[7]

            if gtf_chrom == chrom and feature == "CDS" and exon_start == gtf_start and exon_end == gtf_end:
                phase = int(gtf_phase) 
                return phase 
    return None

def determine_nmd_status(gtf_file, chrom, exon1_start, exon1_end, exon2_start, exon2_end, ptc_relative_position, strand):
    """
    Determines if a PTC triggers NMD by checking the last exon-exon junction downstream.

    Parameters:
        gtf_file (str): Path to the GTF file.
        chrom (str): Chromosome of the transcript.
        exon1_start (int): Start position of exon1.
        exon1_end (int): End position of exon1.
        exon2_start (int): Start position of exon2.
        exon2_end (int): End position of exon2.
        ptc_relative_position (int): PTC position relative to modified sequence.
        strand (str): Strand of the transcript ("+" or "-").

    Returns:
        bool: True if NMD is triggered (PTC > 50 nt upstream of last downstream junction), False otherwise.
    """
    
    exon_ends = []  # Stores exon end positions

    # Convert PTC relative position to genomic coordinate
    if strand == "+":
        ptc_position = exon1_start + ptc_relative_position
    else:
        ptc_position = exon1_end - ptc_relative_position

    # Read the GTF file to extract exon junctions
    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue  # Skip comments
            columns = line.strip().split("\t")
            
            # Extract necessary fields
            gtf_chrom = columns[0]
            feature = columns[2]
            start = int(columns[3])
            end = int(columns[4])
            gtf_strand = columns[6]

            # Ensure correct chromosome and strand
            if gtf_chrom == chrom and feature == "exon" and gtf_strand == strand:
                exon_ends.append(end)

    # Sort exon ends in ascending order for "+", descending for "-"
    exon_ends.sort(reverse=(strand == "-"))

    # Find first exon junction after the PTC
    last_downstream_junction = None
    for exon_end in exon_ends:
        if (strand == "+" and exon_end > ptc_position) or (strand == "-" and exon_end < ptc_position):
            last_downstream_junction = exon_end
            break

    # If no downstream junction found, return False (NMD unlikely)
    if last_downstream_junction is None:
        return False  

    # Calculate distance from PTC to this junction
    distance = abs(last_downstream_junction - ptc_position)

    # NMD is triggered if distance is ≥ 50 nucleotides
    return distance >= 50


def process_vast_output(vast_file, genome, gtf_file):
    df_vast = pd.read_csv(vast_file, sep="\t")
    results = []  

    for index, row in df_vast.iterrows():
        if "IR" in str(row['COMPLEX']):  
            chrom, full_coords = row['FullCO'].split(":", 1)
            
            #try:
                #exon1_coords, exon2_coords = full_coords.split("=")  
            #except ValueError:
                #print(f"Skipping row {index} due to incorrect format: {full_coords}")
                #continue  
            
            #exon1_start, exon1_end = map(int, exon1_coords.split("-"))
            #exon2_coords = exon2_coords.split(":")[0]  
            #exon2_start, exon2_end = map(int, exon2_coords.split("-"))

            #intron_start, intron_end = map(int, row['COORD'].split(":")[1].split("-"))

            #intron_length=intron_end-intron_start

            # Extract sequences
            #exon1_seq = extract_sequence_from_genome(genome, chrom, exon1_start, exon1_end)
            #intron_seq = extract_sequence_from_genome(genome, chrom, intron_start, intron_end)
            #exon2_seq = extract_sequence_from_genome(genome, chrom, exon2_start, exon2_end)
            
            strand = full_coords[-1]

            # Modify the sequence based on the strand
            if strand == "+":
                # For positive strand: exon1 → intron → exon2
                exon1_coords, exon2_coords = full_coords.split("=")


                exon1_start, exon1_end = map(int, exon1_coords.split("-"))
                exon2_coords = exon2_coords.split(":")[0]  
                exon2_start, exon2_end = map(int, exon2_coords.split("-")) 

                intron_start, intron_end = map(int, row['COORD'].split(":")[1].split("-"))

                intron_length=intron_end-intron_start


                exon1_seq = extract_sequence_from_genome(genome, chrom, exon1_start, exon1_end)
                intron_seq = extract_sequence_from_genome(genome, chrom, intron_start, intron_end)
                exon2_seq = extract_sequence_from_genome(genome, chrom, exon2_start, exon2_end)


                modified_dna_sequence = exon1_seq + intron_seq + exon2_seq
                exon_phase = get_exon_phase_from_gtf(gtf_file, chrom, exon1_start, exon1_end)
            elif strand == "-":
                # For negative strand: exon2 → intron → exon1, then reverse complement
                exon1_coords, exon2_coords = full_coords.split("=") 


                exon1_start, exon1_end = map(int, exon1_coords.split("-"))
                exon2_coords = exon2_coords.split(":")[0]  
                exon2_start, exon2_end = map(int, exon1_coords.split("-")) 

                intron_start, intron_end = map(int, row['COORD'].split(":")[1].split("-"))

                intron_length=intron_end-intron_start

            
                exon1_seq = extract_sequence_from_genome(genome, chrom, exon1_start, exon1_end)
                intron_seq = extract_sequence_from_genome(genome, chrom, intron_start, intron_end)
                exon2_seq = extract_sequence_from_genome(genome, chrom, exon2_start, exon2_end)

                exon1_seq=Seq(exon1_seq).complement()               
                intron_seq=Seq(intron_seq).complement()               
                exon2_seq=Seq(exon2_seq).complement()               
                
                
                intron_seq= intron_seq[::-1]
                exon1_seq= exon1_seq[::-1]
                exon2_seq= exon2_seq[::-1]
                
                modified_dna_sequence = exon1_seq + intron_seq + exon2_seq

                exon_phase = get_exon_phase_from_gtf(gtf_file, chrom, exon1_start, exon1_end)


            
            # Validate phase (must be 0, 1, or 2)
            if exon_phase is None:
                print(f"ERROR: Invalid phase detected. Exon1 Phase: {exon_phase}")
                continue

            # Check for premature stop codons
            is_ptc, stop_codon_position = check_for_ptc_in_sequence(modified_dna_sequence, exon_phase, strand)

            frameshift, new_phase = check_frameshift(exon_phase, intron_length)


            if is_ptc:
                # Determine NMD status
                nmd_status = determine_nmd_status(gtf_file, chrom, exon1_start, exon1_end, exon2_start, exon2_end, stop_codon_position, strand)

                results.append({
                    "Gene": row['GENE'],
                    "Event": row['EVENT'],
                    "Coordinates": row['COORD'],
                    "PTC_Position": stop_codon_position,
                    "Exon_Phase": exon_phase,
                    "NMD": nmd_status,
                    "Frameshift": frameshift,  # True if it disrupts the sequence
                    "New_Phase": new_phase # True or False
                })

    return results

def load_genome(fasta_file):
    genome = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome[record.id] = record.seq
    return genome


def main():
    vast_file = "/home/joaossantos/imm/files/vast/output_dir/compareOutput_toPlot.tsv"  
    fasta_file = "/mnt/nfs/lobo/IMM-NFS/imm/files/references/mm10.fa"
    gtf_file = "/mnt/nfs/lobo/IMM-NFS/imm/files/references/mm10.ensGene.gtf"
    output_csv = "/mnt/nfs/lobo/IMM-NFS/imm/files/ptc/ptc_results.csv"  

    # Load the genome
    genome = load_genome(fasta_file)

    # Process the VAST file and generate results
    result = process_vast_output(vast_file, genome, gtf_file)

    

    # Save results to CSV
    df_results = pd.DataFrame(result)
    df_results.to_csv(output_csv, index=False)

    print(f"Results saved to {output_csv}")

