import json
import sys

from tqdm import tqdm
from chembl_webresource_client.new_client import new_client


def printerr(line):
        sys.stderr.write(u'%s\n' % line)

assay = new_client.assay

assay_results = assay.filter(
        description__contains='GSK_PKIS',
        assay_type='B',
        relationship_type='D').only(['document_chembl_id'])

printerr('Number of Assay Search Results: ' + str(len(assay_results)))

doc_list = list()
for e in assay_results:
        doc_list.append(e['document_chembl_id'])

doc_list = list(set(doc_list))

printerr('Documents found in search: ' + str(doc_list))

printerr('Filtering activities with doc id in previous query...')

activity = new_client.activity

res = activity.filter(
        document_chembl_id__in=doc_list, 
        assay_type='B',
        relationship_type='D',
        assay_description__contains='at 1 uM',
        value__gte=50).only(['target_chembl_id', 'molecule_chembl_id'])

printerr('Number of Bioactivity Search Results: ' + str(len(res)))

for interaction in tqdm(res):
        print(u"%s\t%s" % (interaction['target_chembl_id'], interaction['molecule_chembl_id']))
