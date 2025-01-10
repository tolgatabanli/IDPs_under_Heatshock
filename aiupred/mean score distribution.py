import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

with open('aiupred_scores.p', 'rb') as fp:
    proteins = pickle.load(fp)

# find the mode of scores in each protein with 10 bins
for prot, scores in proteins.items():
    # breakpoint()
    stat, bin_edges, binnumber = stats.binned_statistic(scores, scores, statistic='count', bins=20)

    # plt.hist(scores)
    # plt.show()

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    proteins[prot] = bin_centers[np.argmax(stat)]

    # print(bin_centers)
    # print(stat)
    # print(bin_centers[np.argmax(stat)])
    # plt.bar(bin_centers, stat, width=np.diff(bin_edges), align='center', edgecolor='black')
    # plt.show()

plt.hist(proteins.values(), bins=50)
plt.xlabel("Disorder scores")
plt.title("Distribution of mode disorder scores")
plt.savefig("hist (bins_50) of binned(20) mode disorder scores.png")
plt.show()

# save the dictionary in a tsv (protein_name, mode)
csv = "\n".join([prot + '\t' + str(mode) for prot, mode in proteins.items()])
with open("aiupred_modes_of_proteins.tsv", 'w') as output:
    output.write("uniprot\tmode\n")
    output.write(csv)

# Means
with open('aiupred_scores.p', 'rb') as fp:
    proteins = pickle.load(fp)
for prot, scores in proteins.items():
    proteins[prot] = np.mean(scores)
csv = "\n".join([prot + '\t' + str(mean) for prot, mean in proteins.items()])
with open("aiupred_means_of_proteins.tsv", 'w') as output:
    output.write("uniprot\tmean\n")
    output.write(csv)

plt.hist(proteins.values(), bins=50)
plt.xlabel("Distribution of mean disorder scores")
plt.title("Distribution of mean disorder scores")
plt.savefig("hist (bins_50) of mean disorder scores.png")
plt.show()


# By percent disorder with threshold 0.5
with open('aiupred_scores.p', 'rb') as fp:
    proteins = pickle.load(fp)
for prot, scores in proteins.items():
    proteins[prot] = np.sum(scores > 0.5) / scores.size
csv = "\n".join([prot + '\t' + str(ratio) for prot, ratio in proteins.items()])
with open("aiupred_perc_of_proteins_alpha50.tsv", 'w') as output:
    output.write("uniprot\tperc\n")
    output.write(csv)

plt.hist(proteins.values(), bins=50)
plt.xlabel("Percent disorder")
plt.title("Distribution of percent$_{\\alpha=50}$ disorders of proteins")
plt.savefig("hist (bins_50) of percent_alpha50 disorder scores.png")
plt.show()

# By percent disorder with threshold 0.8
with open('aiupred_scores.p', 'rb') as fp:
    proteins = pickle.load(fp)
for prot, scores in proteins.items():
    proteins[prot] = np.sum(scores > 0.8) / scores.size
csv = "\n".join([prot + '\t' + str(ratio) for prot, ratio in proteins.items()])
with open("aiupred_perc_of_proteins_alpha80.tsv", 'w') as output:
    output.write("uniprot\tperc\n")
    output.write(csv)

plt.hist(proteins.values(), bins=50)
plt.xlabel("Percent disorder")
plt.title("Distribution of percent$_{\\alpha=80}$ disorders of proteins")
plt.savefig("hist (bins_50) of percent_alpha80 disorder scores.png")
plt.show()

# By percent disorder with threshold 0.6
with open('aiupred_scores.p', 'rb') as fp:
    proteins = pickle.load(fp)
for prot, scores in proteins.items():
    proteins[prot] = np.sum(scores > 0.6) / scores.size
csv = "\n".join([prot + '\t' + str(ratio) for prot, ratio in proteins.items()])
with open("aiupred_perc_of_proteins_alpha60.tsv", 'w') as output:
    output.write("uniprot\tperc\n")
    output.write(csv)

plt.hist(proteins.values(), bins=50)
plt.xlabel("Percent disorder")
plt.title("Distribution of percent$_{\\alpha=60}$ disorders with AIUPred")
plt.savefig("hist (bins_50) of percent_alpha60 disorder aiupred.png")
plt.show()
