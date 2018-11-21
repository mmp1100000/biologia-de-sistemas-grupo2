import re

from chembl_webresource_client import CompoundResource
from chembl_webresource_client.new_client import new_client

target = new_client.activity
tsvname = 'journal.pcbi.1004120.s004.TSV'
tsvdict = dict()

with open(tsvname, mode='r') as tsvfile:
    is_key = False
    for line in tsvfile:
        tab_pos = re.search("\d", line)
        if tab_pos:
            genes = line[tab_pos.start()+1:].strip('\n').split('/')
            tsvdict[line[0:tab_pos.start()].strip()] = genes

breast_genes = [x for d in tsvdict for x in tsvdict[d] if d.startswith("breast")]
print(breast_genes[0])
print(target.search(str(breast_genes[0])))
compounds = CompoundResource()
