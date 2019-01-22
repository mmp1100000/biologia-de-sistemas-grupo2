import fileinput, sys
from chembl_webresource_client.new_client import new_client


def printerr(line):
    sys.stderr.write(u'%s\n' % line)

def get_accession_from_target(target):
    return target[0]['target_components'][0]['accession']


def chembl_id_to_accession_id():
    """
    Returns a dictionary of chembl_id:accession_id pairs, given the name of the interactions file
    @TODO FIND OUT WHY SOME ACCESSION IDS CAN'T BE FOUND
    :param interactions_file: interactions file
    :return: dict
    """
    target = new_client.target
    accession_ids = list()
    for chembl_id in fileinput.input():
        chembl_id = chembl_id.strip()
        printerr(chembl_id)
        try:
            res2 = target.filter(target_chembl_id=chembl_id)
            accession_ids.append(get_accession_from_target(res2))
        except Exception:
            accession_ids.append(" ")
        # res2 = target.filter(target_chembl_id=chembl_id)
        # accession_ids.append(get_accession_from_target(res2))
    return accession_ids


accession_list = chembl_id_to_accession_id()

for gene_id in accession_list:
    print(gene_id)