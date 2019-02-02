from scipy import stats
from utils import printerr, suppress_stderr
import pandas as pd
import sys


"""
    Este script aplica el test binomial sobre las interacciones, en las que debemos tener marcadas cuales
    pertenecen o no a la enfermedad. Devolviendonos una lista de los p-valores asociados a cada interacción.
    
"""

# Comprobamos que el archivo tenga las dimensiones correctas
if (len(sys.argv) < 4):
    exit("You must provide the flagged interactions file, the confidence interval index and the number of total genes associated to the disease as argument.")
# Parseamos las variables de la linea de comandos
fichero = sys.argv[1]
conf_level = float(sys.argv[2])
number_assoc_genes = int(sys.argv[3])

# Leemos el archivo
df = pd.read_csv(fichero, header=None, sep="\t")
res_pval = {}
resul = []

# Preparamos los datos para el test binomial
# Establecemos la probabilidad siguiendo la corrección de Bonferroni
probability = number_assoc_genes / 20000 # all human proteins
# Obtenemos la lista de fármacos
unique_drugs = set(df[0])
# Obtenemos la lista de targets
breast_cancer = df[df[2] == 1]
unique_genes = set(breast_cancer[1])

printerr("Disease-assocciated genes interacting: ")
printerr(unique_genes)

unique_drugs_interacting = set(breast_cancer[0])
printerr("Compounds interacting with disease-assocciated genes: " + str(len(unique_drugs_interacting)))

# Pasamos a calcular los p-valores con el test binomial
for i in list(unique_drugs):
    df_moment = df[df[0] == i]
    # Vemos el número total de interacciones para un fármaco
    trials = len(list(df_moment[1]))
    success_rows = df[(df[0] == i) & (df[2] == 1)]
    # Vemos el número total de interacciones con cáncer de mama
    success = len(success_rows[1])
    # Evaluamos el test binomial
    p_value = stats.binom_test(success, n=trials, p=probability, alternative='greater')
    res_pval.update({str(i): p_value})
for p in df[0]:
    # Asignamos a los diferentes fármacos
    resul.append(res_pval.get(p))

# Guardamos los resultados obtenidos
df = df.assign(p_value=pd.Series(resul).values)
df = df.sort_values(by=['p_value'])
df.to_csv("outputs/pvalued_interactions_" + sys.argv[2] + ".tsv", sep="\t", encoding='utf-8')
printerr(str(len(set(df[df["p_value"].values < (1 - conf_level)/len(unique_drugs)][0]))) + " out of " + str(
len(list(set(df[0])))) + " compounds allow us to reject the null hypothesis with a Bonferroni corrected p-value of " + str((1 - conf_level)/len(unique_drugs)) + ":")

for compound in set(df[df["p_value"].values < (1 - conf_level)/len(unique_drugs)][0]):
    printerr(compound)
