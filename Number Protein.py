import re
from tkinter import filedialog
def count_proteins_in_fasta(fasta_file):
    count = 0
    with open(fasta_file) as f:
        a = re.compile(">\w*")
        all = f.readlines()
        for aa in all:
            if a.match(aa):
                count += 1
    return count

fasta_file = filedialog.askopenfilename()
protein_count = count_proteins_in_fasta(fasta_file)

print(f"Number of proteins in the FASTA file: {protein_count}")

