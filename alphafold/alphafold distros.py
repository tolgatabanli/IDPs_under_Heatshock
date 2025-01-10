import pandas as pd
import matplotlib.pyplot as plt

with open("../alphafold/local_plddts.csv") as f:
    prots = pd.read_csv(f)

prots = prots["x"].to_numpy()

plt.hist(prots, bins=50)
plt.xlabel("pLDDT scores")
plt.title("Distribution of all local pLDDT scores from AlphaFold")
plt.savefig("../alphafold/hist (bins_50) of all plddt scores.png")

