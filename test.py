import json
from pprint import pprint

from conda._vendor.tqdm import tqdm

from chembl_webresource_client.new_client import new_client

assay = new_client.assay
res = assay.search('GSK_PKIS')

for r in res:
    pprint(r)