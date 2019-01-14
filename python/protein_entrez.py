from pprint import pprint

from Bio import Entrez
from tqdm import tqdm


def get_gene_id(id_string):
    """
    Given a string of accession ids separated by commas, returns its gene IDs as a list.
    :param id_string:
    :return:
    """
    gene_id_list = list()
    handle = Entrez.efetch(db="protein", id=id_string, retmode="xml")
    for record in tqdm(Entrez.read(handle)):
        try:
            entr = record['GBSeq_source-db']
            gene_id_raw = entr[entr.find('GeneID'):]
            gene_id = int(gene_id_raw[gene_id_raw.find(':') + 1:gene_id_raw.find(',')])
            gene_id_list.append(gene_id)
        except Exception:
            gene_id_list.append(None)
    handle.close()
    return gene_id_list



