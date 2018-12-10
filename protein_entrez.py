from Bio import Entrez

handle = Entrez.efetch(db="protein", id="P42681", retmode="xml")
records = Entrez.parse(handle)
for record in records:
    print(record)
handle.close()