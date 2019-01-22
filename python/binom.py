from scipy import stats
from utils import printerr, suppress_stderr
import pandas as pd
import sys


if (len(sys.argv) < 3):
        exit("You must provide the flagged interactions file and the confidence interval percentage as argument.")
inputfile = sys.argv[1]
conf_level = float(sys.argv[2])/100


@suppress_stderr
def binomial(fichero, conf_level):
    df = pd.read_csv(fichero, header=None, sep="\t")
    res_pval = {}
    resul = []
    for i in list(set(df[0])):
        df_moment = df[df[0] == i]
        total_nodes = len(list(df_moment[1]))
        mama = df_moment[df[2] == 1]
        mama_nodes = len(mama[1])
        p_value = stats.binom_test(mama_nodes, n=total_nodes, p=1-conf_level, alternative='two-sided')
        res_pval.update({str(i): p_value})
    for p in df[0]:
        resul.append(res_pval.get(p))

    df = df.assign(p_value = pd.Series(resul).values)
<<<<<<< HEAD
    #df = df[df["p_value"].values < 0.05]
    df.to_csv("../outputs/pvalued_interactions_"+str(conf_level)+".tsv", sep="\t",  encoding='utf-8')
    return(df)


conf_level = 0.85
res = binomial(fichero="../outputs/flagged_interactions.tsv", conf_level = conf_level)
print("Tenemos un total de " + str(len(list(set(res[res["p_value"].values < 0.05][0])))) + " de " + str(len(list(set(res[0])))) + " fármacos, que actúan con un " + str(conf_level) + " de confianza sobre la familia funcional de cancer de mama.")
=======
    df.to_csv("outputs/pvalued_interactions_"+str(conf_level)+".tsv", sep="\t",  encoding='utf-8')
    return(df)


res = binomial(fichero=inputfile, conf_level = conf_level)
printerr(str(len(list(set(res[res["p_value"].values < 0.05][0])))) + " out of " + str(len(list(set(res[0])))) + " drugs allow us to reject the null hypothesis with a confidence interval of " + str(conf_level) + ".")
>>>>>>> f5ee13dc467357e69488d90ceaa2e643db02bb3a
