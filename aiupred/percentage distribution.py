import pickle
import matplotlib.pyplot as plt
import numpy as np

with open('aiupred_scores.p', 'rb') as fp:
    proteins = pickle.load(fp)

percentages = dict()
alpha = 0.5
for name, disorder_array in proteins.items():
    percentages[name] = (disorder_array > alpha).sum() / len(disorder_array)

# Logit Transformation
# arr = [np.atleast_1d(perc) for perc in percentages.values()]
# arr = np.concatenate(arr)
# clipped = np.clip(arr, a_min=1e-6, a_max=1 - 1e-6)
# logit_values = np.log(clipped / (1 - clipped))

arr = [np.atleast_1d(perc) for perc in percentages.values()]
arr = np.concatenate(arr)
plt.hist(arr, bins=50)
plt.xlabel("Percentage of residues > " + str(alpha))
plt.title("Distribution of percentages of proteins")
plt.show()
#plt.savefig("hist (bins_50) of percentage of proteins of residues alpha_0.5.png")
