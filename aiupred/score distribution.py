from numpy import mean, array, concatenate
import matplotlib.pyplot as plt
import pickle as pickle

import aiupred_lib
from Bio import SeqIO

# Load the models and let AIUPred find if a GPU is available.
embedding_model, regression_model, device = aiupred_lib.init_models()
# Predict disorder of a sequence
sequence = ''

proteins = dict()

fasta_file = r"C:\Users\dell\Desktop\R with Yeast\Second Semester\uniprotkb_proteome_UP000002311_2024_11_29.fasta"
fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')

# i = 0
# for fasta in fasta_sequences:
#     name, sequence = fasta.id.split('|')[-1], str(fasta.seq)
#     proteins[name] = aiupred_lib.predict_disorder(sequence, embedding_model, regression_model, device)
#     i += 1
#     print(i)

# with open('aiupred_scores.p', 'wb') as pickle_file:
#     pickle.dump(proteins, pickle_3file, protocol=pickle.HIGHEST_PROTOCOL)

with open('aiupred_scores.p', 'rb') as fp:
    proteins = pickle.load(fp)

vals = concatenate(list(proteins.values()))
plt.hist(vals, bins=100)
plt.xlabel("Disorder scores")
plt.title("Distribution of all disorder scores from AIUPred")
plt.savefig("hist (bins_100) of all disorder scores.png")
