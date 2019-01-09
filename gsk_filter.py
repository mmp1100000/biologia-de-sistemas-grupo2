import json

from conda._vendor.tqdm import tqdm

from chembl_webresource_client.new_client import new_client

assay = new_client.assay
res = assay.search('GSK_PKIS')

print('Number of Assay Search Results: ' + str(len(res)))

activity = new_client.activity
doc_list = list()
for e in res:
    if e['assay_type'] == 'B' and e['relationship_type'] == 'D':
        doc_list.append(e['document_chembl_id'])
doc_list = list(set(doc_list))
print('Documents found in search: ' + str(doc_list))

print('Filtering activities with doc id in previous query...')

res = activity.filter(document_chembl_id__in=doc_list)
print('Number of Bioactivity Search Results: ' + str(len(res)))
input()
# pprint(res[0])
file1 = open('targets.txt', 'w')
file2 = open('compounds.txt', 'w')
for interaction in tqdm(res):
    if interaction['assay_type'] == 'B' and interaction['pchembl_value'] is not None and interaction[
        'pchembl_value'] >= 6:
        file1.write(interaction['target_chembl_id'] + ',')
        file2.write(interaction['molecule_chembl_id'] + ',')
    # file.write(interaction['molecule_chembl_id'] + '-' + interaction['target_chembl_id'] + '\n')
    # Esto devuelve todos los tagets > Tenemos que quedarnos con el GSK y el target

activities = new_client.target
res2 = activities.filter(target_chembl_id="CHEMBL4367",
                         pchembl_value__isnull=False)
