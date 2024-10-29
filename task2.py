from Bio import SeqIO
from Bio.Seq import Seq
import random


def read_bed_file(bed_file):
    """Parse BED file to extract regions for each scaffold."""
    bed_regions = {}
    with open(bed_file, 'r') as f:
        for line in f:
            if not line.strip() or len(line.strip().split()) != 3:
                continue
            scaffold, start, end = line.strip().split()
            start, end = int(start), int(end)
            if scaffold not in bed_regions:
                bed_regions[scaffold] = []
            bed_regions[scaffold].append((start, end))
    return bed_regions


def extract_subsequences(fasta_file, bed_regions):
    """Extract sequences based on BED regions and remove them from the primary FASTA."""
    extracted_sequences = {}
    primary_sequences = {}

    for record in SeqIO.parse(fasta_file, 'fasta'):
        scaffold = record.id
        sequence = str(record.seq)

        if scaffold in bed_regions:
            primary_sequences[scaffold] = list(sequence)  # Make modifiable list

            for start, end in bed_regions[scaffold]:
                subseq = sequence[start:end]
                extracted_sequences[(scaffold, start, end)] = subseq

                # Remove the subsequence from the primary sequence
                for i in range(start, end):
                    primary_sequences[scaffold][i] = '-'

            # Convert list back to a string and remove placeholders
            primary_sequences[scaffold] = "".join(primary_sequences[scaffold]).replace('-', '')
        else:
            # Keep the original sequence if no BED regions for this scaffold
            primary_sequences[scaffold] = sequence

    return extracted_sequences, primary_sequences


def reverse_complement_sequences(extracted_sequences):
    """Get reverse complement of each extracted sequence."""
    rev_comp_sequences = {}
    for key, seq in extracted_sequences.items():
        rev_comp_sequences[key] = str(Seq(seq).reverse_complement())
    return rev_comp_sequences


def introduce_mutation(sequences):
    """Introduce a random single-nucleotide mutation in each sequence."""
    mutated_sequences = {}
    for key, seq in sequences.items():
        seq_list = list(seq)
        pos = random.randint(0, len(seq_list) - 1)
        original_base = seq_list[pos]
        bases = ['A', 'T', 'C', 'G']
        bases.remove(original_base)
        seq_list[pos] = random.choice(bases)
        mutated_sequences[key] = "".join(seq_list)
    return mutated_sequences


def insert_sequences_back(primary_sequences, processed_sequences):
    """Insert processed sequences back into their original locations in the primary sequence."""
    updated_sequences = primary_sequences.copy()

    for (scaffold, start, end), seq in processed_sequences.items():
        if scaffold in updated_sequences:
            # Insert sequence back by reconstructing full string with processed subsequence
            updated_sequences[scaffold] = (
                    updated_sequences[scaffold][:start] + seq + updated_sequences[scaffold][start:]
            )
    return updated_sequences


def write_to_fasta(output_file, sequences):
    """Write sequences to a FASTA file."""
    with open(output_file, 'w') as f:
        for scaffold, seq in sequences.items():
            f.write(f'>{scaffold}\n')
            f.write(f'{seq}\n')


def process_fasta_with_bed(fasta_file, bed_file, output_file):
    """Main function to process FASTA using BED regions and perform modifications."""
    # Read and parse BED file
    bed_regions = read_bed_file(bed_file)

    # Extract sequences and modify the primary sequences accordingly
    extracted_sequences, primary_sequences = extract_subsequences(fasta_file, bed_regions)

    # Reverse complement the extracted sequences
    rev_comp_sequences = reverse_complement_sequences(extracted_sequences)

    # Introduce mutations into the reverse-complemented sequences
    processed_sequences = introduce_mutation(rev_comp_sequences)

    # Insert the processed sequences back into the primary sequence
    updated_sequences = insert_sequences_back(primary_sequences, processed_sequences)

    # Write the updated sequences to an output FASTA file
    write_to_fasta(output_file, updated_sequences)


process_fasta_with_bed('genome.fasta', 'regions.bed', 'updated_genome.fasta')