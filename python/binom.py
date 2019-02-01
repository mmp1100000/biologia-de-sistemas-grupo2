from scipy import stats
from utils import printerr, suppress_stderr
import pandas as pd
import sys

if (len(sys.argv) < 4):
    exit("You must provide the flagged interactions file, the confidence interval index and the number of total genes associated to the disease as argument.")
fichero = sys.argv[1]
conf_level = float(sys.argv[2])
number_assoc_genes = int(sys.argv[3])


df = pd.read_csv(fichero, header=None, sep="\t")
res_pval = {}
resul = []
probability = number_assoc_genes / 20000 # all human proteins
unique_drugs = set(df[0])

breast_cancer = df[df[2] == 1]
unique_genes = set(breast_cancer[1])
printerr("Disease-assocciated genes interacting: ")
printerr(unique_genes)

unique_drugs_interacting = set(breast_cancer[0])
printerr("Compounds interacting with disease-assocciated genes: " + str(len(unique_drugs_interacting)))

for i in list(unique_drugs):
    df_moment = df[df[0] == i]
    trials = len(list(df_moment[1]))
    success_rows = df[(df[0] == i) & (df[2] == 1)]
    success = len(success_rows[1])
    p_value = stats.binom_test(success, n=trials, p=probability, alternative='greater')
    res_pval.update({str(i): p_value})
for p in df[0]:
    resul.append(res_pval.get(p))

df = df.assign(p_value=pd.Series(resul).values)
df = df.sort_values(by=['p_value'])
df.to_csv("outputs/pvalued_interactions_" + sys.argv[2] + ".tsv", sep="\t", encoding='utf-8')
printerr(str(len(set(df[df["p_value"].values < (1 - conf_level)/len(unique_drugs)][0]))) + " out of " + str(
len(list(set(df[0])))) + " compounds allow us to reject the null hypothesis with a Bonferroni corrected p-value of " + str((1 - conf_level)/len(unique_drugs)) + ":")

for compound in set(df[df["p_value"].values < (1 - conf_level)/len(unique_drugs)][0]):
    printerr(compound)
