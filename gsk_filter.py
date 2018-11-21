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
print('Number of Bioactivity Search Results: '+str(len(res)))
