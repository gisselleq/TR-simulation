import pandas as pd
import msprime
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
import io

num_people = 1000
population_size = 10000
gen_list = [5000]
mutation_rate=  1e-4
# length of the genome sequence to simulate
sequence_length = 1
num_loci = 1000


model_specs = {
  "SMM": (msprime.SMM, {}),
    "TPM lo": (msprime.TPM, {"p": 0.9, "m": 0.5}),
    "TPM hi": (msprime.TPM, {"p": 0.1, "m": 0.1}),
  #  "EL2 base": (msprime.EL2,{"m": 0.5, "u": 0.5, "v": 0}),
   # "EL2 cont": (msprime.EL2, {"m":0.5 ,"u": 0.4, "v": 0.05}),
    #"EL2 bigjump": (msprime.EL2, {"m": 0.1, "u": 0.5, "v": 0}),
}

results = {}
models = {}
labels = {}

for name, (clas, params) in model_specs.items():
    models[name] = clas(**params)
    param_str = " ".join(f"{k}={v}" for k, v in params.items())
    labels[name] = f"{name} {param_str}".strip()


# add mutation models array here

for model_name, model in models.items():

    for num_gen in gen_list:

        tslist = msprime.sim_ancestry(
            samples=num_people,
            population_size=population_size,
            sequence_length=sequence_length,
            end_time=num_gen,
            num_replicates=num_loci,
            #random_seed=2,
        )
        het = []
        allele_freq = []
        num_alleles = []
        genotypes = []

        #vcf_path = f"{model_name.lower()}_{num_gen}.vcf"
        #wrote_header = False
        #skipped = []

        for i, ts in enumerate(tslist):
            mts = msprime.sim_mutations(
                ts, rate=mutation_rate, model=model, #random_seed=i + 1
            )
            # var alleles list the diff repeat counts seen at site
            # var geno stores which allele each chromosome has
            if mts.num_sites == 0:
                het.append(0)
                allele_freq.append(1)
                num_alleles.append(1)

                continue

            variant = next(mts.variants())
            g = variant.genotypes

            het_count = 0

            for person in range(num_people):
                copy1 = g[2 * person]
                copy2 = g[2 * person + 1]
                if copy1 != copy2:
                    het_count += 1

            het.append(het_count / num_people)
	        
            # counts how many time each index appears
            counts = np.bincount(g)
            num_alleles.append(int((counts > 0).sum()))

            # adds each freq to the list
            allele_freq.extend(counts[counts > 0] / len(g))

            #if len(variant.alleles) > 9:
            #    skipped.append(i)
            #else:
            #    buf = io.StringIO()
            #    mts.write_vcf(buf, contig_id=f"locus{i}", allow_position_zero=True)
             #   lines = buf.getvalue().splitlines(keepends=True)
             #   with open(vcf_path, "a" if wrote_header else "w") as f:
             #       if wrote_header:
             #           f.writelines(l for l in lines if not l.startswith("#"))
             #       else:
             #           f.writelines(lines)
              #  wrote_header = True
        print(model_name, num_gen, np.mean(het))
       # print(f"saved {vcf_path}, skipped {len(skipped)} loci with >9 alleles")
        results[(labels)[model_name], num_gen] = het

        # plots graphs inside loop
        #fig, (ax1) = plt.subplots(1, 1, figsize=(10, 4))
       # ax1.hist(het, bins=20, range=(0, 1))
        #ax1.set_title("heterozygosity per locus")
        # ax2.hist(allele_freq, bins = 20, range =(0,1))
        # ax2.set_title("allele freq")
        # lo, hi = min(num_alleles), max(num_alleles)
        # ax3.hist(num_alleles, bins = np.arange(lo - 0.5, hi + 1.5, 1))
        #  ax3.set_title("alleles per locus")
        # ax3.set_xlabel("number of distinct alleles")
        # ax3.set_ylabel("number of loci")

       # plt.savefig(f"{model_name.lower()}_{num_gen}.png", dpi=500)
        #print(f"saved {model_name.lower()}{num_gen}gen.png")



rows = []
for (label, num_gen), het_list in results.items():
    for h in het_list:
        rows.append({"mutation model": f"{label} {num_gen}", "het": h})

df = pd.DataFrame(rows)


fig, ax = plt.subplots(figsize=(10, 5))
sns.histplot(
    data=df,
    x="het",
    hue="mutation model",
    bins=20,
    binrange=(0, 1),
    ax=ax,
)


ax.set_title("het per locus")
ax.set_xlabel("heterozygosity")
plt.savefig("het_allmm.png", dpi=300)
plt.close()
