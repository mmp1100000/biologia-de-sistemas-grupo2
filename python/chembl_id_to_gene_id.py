# from pprint import pprint
import pickle
# from tqdm import tqdm
#
# from chembl_webresource_client.new_client import new_client
# from drug_gene import get_tsv_gene_ids
from protein_entrez import get_gene_id

#
# chembl_file = '../outputs/interactions.tsv'
#
#
# def get_accession_from_target(target):
#     return target[0]['target_components'][0]['accession']
#
#
# def chembl_id_to_accession_id(chembl_file):
#     """
#     Returns a dictionary of chembl_id:accession_id pairs, given the name of the interactions file
#     @TODO FIND OUT WHY SOME ACCESSION IDS CAN'T BE FOUND
#     :param chembl_file: interactions file
#     :return: dict
#     """
#     target = new_client.target
#     accession_ids = dict()
#     with open(chembl_file, 'r') as file:
#         for line in tqdm(file):
#             chembl_id = line.split('\t')[0]
#             if chembl_id not in accession_ids.keys():
#                 try:
#                     res2 = target.filter(target_chembl_id=chembl_id)
#                     accession_ids[chembl_id] = get_accession_from_target(res2)
#                 except Exception:
#                     accession_ids[chembl_id] = None
#     return accession_ids
#
#
# accession_dict = chembl_id_to_accession_id(chembl_file)
# accession_list = [x for x in list(accession_dict.values()) if x is not None]  # Remove None values
#
# pickle_out2 = open("dict", "wb")
# pickle.dump(accession_dict, pickle_out2)
# pickle_out2.close()
#
# pickle_out = open("list", "wb")
# pickle.dump(accession_list, pickle_out)
# pickle_out.close()
from python.protein_entrez import get_gene_id

file_list = open("list", "rb")
accession_list = pickle.load(file_list)
file_list.close()

gene_list = get_gene_id(','.join(accession_list))  # Join the accession list into a single string

pickle_out = open("genes", "wb")
pickle.dump(gene_list, pickle_out)
pickle_out.close()

# tsv_genes = get_tsv_gene_ids()  # Get gene ids from tsv
#
# for found_gene in gene_list:  # Compare found genes with tsv genes
#     if found_gene in tsv_genes:
#         print(found_gene)