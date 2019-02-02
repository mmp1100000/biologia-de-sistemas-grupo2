from chembl_webresource_client.new_client import new_client
from utils import printerr
from tqdm import tqdm
import sys

if (len(sys.argv) < 2):
    exit("You must provide the minimum inhibition percentage at 1 micromol as argument.")
min_inhibition_percentage = int(sys.argv[1])

assay = new_client.assay

# Filtramos los assay para quedarnos solo con los GSK
assay_results = assay.filter(
    description__contains='GSK_PKIS',
    assay_type='B',
    relationship_type='D').only(['document_chembl_id'])

printerr('Number of Assay Search Results: ' + str(len(assay_results)))

# Vemos cuantos fármacos diferentes tenemos
doc_list = list()
for e in assay_results:
    doc_list.append(e['document_chembl_id'])

doc_list = list(set(doc_list))

printerr('Documents found in search: ' + str(doc_list))

printerr('Filtering activities with doc id in previous query...')

activity = new_client.activity

# Hacemos otra búsqueda para obtener las interacciones
interactions = activity.filter(
    document_chembl_id__in=doc_list,
    assay_type='B',
    relationship_type='D',
    assay_description__contains='at 1 uM', # pchembl = -log(IC50) = -log(10^-6)=6
    value__gte=min_inhibition_percentage).only(['target_chembl_id', 'molecule_chembl_id'])

printerr('Number of Bioactivity Search Results: ' + str(len(interactions)))

for interaction in tqdm(interactions):
    print(u"%s\t%s" % (interaction['target_chembl_id'], interaction['molecule_chembl_id']))
