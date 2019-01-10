import json

from tqdm import tqdm
from chembl_webresource_client.new_client import new_client

assay = new_client.assay

assay_results = assay.filter(
        description__contains='GSK_PKIS',
        assay_type='B',
        relationship_type='D').only(['document_chembl_id'])

print('Number of Assay Search Results: ' + str(len(assay_results)))

doc_list = list()
for e in assay_results:
        doc_list.append(e['document_chembl_id'])

doc_list = list(set(doc_list))

print('Documents found in search: ' + str(doc_list))

print('Filtering activities with doc id in previous query...')

activity = new_client.activity

res = activity.filter(
        document_chembl_id__in=doc_list, 
        assay_type='B', 
        pchembl_value__gte=6).only(['target_chembl_id', 'molecule_chembl_id'])

# print('Number of Bioactivity Search Results: ' + str(len(res)))
# input()
# pprint(res[0])

for interaction in res:
        print(interaction['target_chembl_id'] + '\t' + interaction['molecule_chembl_id'])

# activities = new_client.target
# res2 = activities.filter(target_chembl_id="CHEMBL4367",
#                          pchembl_value__isnull=False)
