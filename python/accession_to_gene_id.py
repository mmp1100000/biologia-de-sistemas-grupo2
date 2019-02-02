from Bio import Entrez
import fileinput

from tqdm import tqdm


Entrez.email = 'sosa@uma.es'

def get_gene_id(id_string):
    """
    Given a string of accession ids separated by commas, returns its gene IDs as a list.
    :param id_string:
    :return:
    """

    gene_id_list = list()
    handle = Entrez.efetch(db="protein", id=id_string, retmode="xml")
    for record in tqdm(Entrez.read(handle)):
        entr = record['GBSeq_source-db']
        gene_id_raw = entr[entr.find('GeneID'):]
        gene_id = int(gene_id_raw[gene_id_raw.find(':') + 1:gene_id_raw.find(',')])
        gene_id_list.append(gene_id)
    handle.close()
    return gene_id_list


for gene_id in get_gene_id(",".join([line.strip() for line in fileinput.input()])):
    print(gene_id)
