import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def generate_random_alignment(num_taxa, seq_length):
    """Generates a random multiple sequence alignment."""
    sequences = []
    for i in range(num_taxa):
        sequence = "".join(random.choice("ACGT") for _ in range(seq_length))
        seq_record = SeqRecord(Seq(sequence), id=f"taxon_{i+1}")
        sequences.append(seq_record)

    alignment = MultipleSeqAlignment(sequences)
    return alignment

if __name__ == "__main__":
    num_taxa = 12
    seq_length = 2000

    alignment = generate_random_alignment(num_taxa, seq_length)

    with open("test_files/generated_test.fasta", "w") as outfile:
        for record in alignment:
            outfile.write(f">{record.id}\n{record.seq}\n")

    print(f"Successfully generated alignment with {num_taxa} taxa and {seq_length} bp length to test_files/generated_test.fasta")
