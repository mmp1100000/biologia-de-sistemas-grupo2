"""

    Test binomial

    Lo vamos a utilizar para sacar el p-valor asociado a cada nodo en base a sus interacciones.

    Para ello hace falta:

        Saber que interacciones son relacionadas con cancer de mama y cuales no
        Saber de cada nodo las interacciones que tiene

        Elegir el porcentage a partir del cual aceptaremos la hipótesis nula.

    :param
        n  : Numero de relaciones que tiene un nodo
        x  : list(Relaciones con cancer de mama, Relaciones que no sean con cancer de mama )
        alternative :  two-sided, greater, less
        p :


The null hypothesis is that the proteins in  are sampled from the same general population as the proteins in
, and thus the probability of observing a target of d as a relative of FF is the same as observing any protein as a
relative of FF i.e. PFF. Since we operate at a confidence level of 0.95, if p-value < 0.05 we reject the null hypothesis
and we consider that the probability of observing the targets of d among the relatives of FF is different from the
probability of observing any set of proteins among the relatives of FF. Therefore, the p-value reported indicates if
observing the targets of d among the relative of FF is likely to happen by chance. For p-values < 0.05 we consider
the corresponding drug-CATH functional family association statistically significant and not likely to happen by
chance.

"""

from scipy import stats
import pandas as pd

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
    df.to_csv("../outputs/pvalued_interactions_"+str(conf_level)+".tsv", sep="\t",  encoding='utf-8')
    return(df)


conf_level = 0.95
res = binomial(fichero="../outputs/flagged_interactions.tsv", conf_level = conf_level)
print("Tenemos un total de " + str(len(list(set(res[res["p_value"].values < 0.05][0])))) + " de " + str(len(list(set(res[0])))) + " fármacos, que actúan con un " + str(conf_level) + " de confianza sobre la familia funcional de cancer de mama.")