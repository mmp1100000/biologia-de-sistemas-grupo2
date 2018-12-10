from pprint import pprint

from Bio import Entrez

handle = Entrez.efetch(db="protein", id="P42681", retmode="xml")
#records = Entrez.parse(handle)
print(str(Entrez.read(handle)[0]['GBSeq_source-db']))['GeneID':]
#for record in records:
#    pprint(record['GBSeq_source-db'])
handle.close()

