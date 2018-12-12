import json

from conda._vendor.tqdm import tqdm

from chembl_webresource_client.new_client import new_client

assay = new_client.assay
res = assay.search('GSK_PKIS')

print('Number of Assay Search Results: ' + str(len(res)))

activity = new_client.activity
doc_list = list()
for e in res:
    doc_list.append(e['document_chembl_id'])
doc_list = list(set(doc_list))
print('Documents found in search: ' + str(doc_list))

print('Filtering activities with doc id in previous query...')

res = activity.filter(document_chembl_id__in=doc_list)
print('Number of Bioactivity Search Results: ' + str(len(res)))
input()
print(res[0])

with open('interactions.txt', 'w') as file:
    for interaction in tqdm(res):
        file.write(interaction['target_chembl_id']+'\n')

# chembl_target_ids = set()

# for interaction in tqdm(res):
#    chembl_target_ids.add(interaction['target_chembl_id'])

# print(chembl_target_ids)

activities = new_client.target
res2 = activities.filter(target_chembl_id="CHEMBL4367",
                         pchembl_value__isnull=False)

for e in res2:
    print(e)