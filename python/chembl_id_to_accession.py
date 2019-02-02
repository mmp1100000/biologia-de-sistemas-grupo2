from chembl_webresource_client.new_client import new_client
from tqdm import tqdm
import fileinput


def get_accession_from_target(target):
    # Sacamos el accession de un taget
    return target[0]['target_components'][0]['accession']


def chembl_id_to_accession_id():
    """"
        Cogemos un fichero con los ChemblID y los traducimos a accesions
    """
    target = new_client.target
    accession_ids = list()
    for chembl_id in tqdm(fileinput.input()):
        chembl_id = chembl_id.strip()
        res = target.filter(target_chembl_id=chembl_id)
        accession_ids.append(get_accession_from_target(res))
    return accession_ids


for accession in chembl_id_to_accession_id():
    print(accession)
