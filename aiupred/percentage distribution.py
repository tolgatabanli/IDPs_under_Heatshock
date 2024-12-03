import pickle
import matplotlib.pyplot as plt

with open('aiupred_scores.p', 'rb') as fp:
    proteins = pickle.load(fp)

percentages = dict()

for name, disorder_array in proteins.items():
    percentages[name] = (disorder_array > 0.5).sum() / len(disorder_array)

plt.hist(percentages.values(), bins=25)
plt.savefig("hist (bins_50) of percentage of proteins of residues alpha_0.5.png")
