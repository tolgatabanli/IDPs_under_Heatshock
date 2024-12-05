import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import csv

with open('aiupred_scores.p', 'rb') as fp:
    proteins = pickle.load(fp)

# find the mode of scores in each protein with 10 bins
for prot, scores in proteins.items():
    # breakpoint()
    stat, bin_edges, binnumber = stats.binned_statistic(scores, scores, statistic='count', bins=30)

    # plt.hist(scores)
    # plt.show()

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    proteins[prot] = bin_centers[np.argmax(stat)]

    # print(bin_centers)
    # print(stat)
    # print(bin_centers[np.argmax(stat)])
    # plt.bar(bin_centers, stat, width=np.diff(bin_edges), align='center', edgecolor='black')
    # plt.show()

# for prot, scores in proteins.items():
#     proteins[prot] = np.median(scores)

plt.hist(proteins.values(), bins=50)
plt.xlabel("Disorder scores")
plt.title("Distribution of binned(10) mode disorder scores from AIUPred")
plt.show()
# plt.savefig("hist (bins_50) of binned(10) mode disorder scores.png")

# save the dictionary in a tsv (protein_name, mode)
csv = "\n".join([prot + '\t' + str(mode) for prot, mode in proteins.items()])
with open("aiupred_modes_of_proteins.tsv", 'w') as output:
    output.write("name_uniprot\tscore_mode\n")
    output.write(csv)


# By majority vote with separation into two categories around a threshold
with open('aiupred_scores.p', 'rb') as fp:
    proteins = pickle.load(fp)

for prot, scores in proteins.items():
    proteins[prot] = np.sum(scores > 0.8) / scores.size

csv = "\n".join([prot + '\t' + str(ratio) for prot, ratio in proteins.items()])
with open("aiupred_majority_vote.tsv", 'w') as output:
    output.write("name_uniprot\tratio\n")
    output.write(csv)

plt.hist(proteins.values(), bins=50)
plt.xlabel("Disorder ratio according to disorder vote")
plt.title("Distribution of disorder ratios")
plt.show()

vals = np.array([*proteins.values()])
plt.hist(vals[vals < 0.1], bins=50)
plt.show()

