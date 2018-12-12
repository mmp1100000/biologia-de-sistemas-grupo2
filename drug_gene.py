import re

tsvname = 'journal.pcbi.1004120.s004.TSV'
tsvdict = dict()


def get_tsv_gene_ids():
    """
    Obtain breast neoplasm genes from journal tsv
    :return: gene list
    """
    with open(tsvname, mode='r') as tsvfile:
        is_key = False
        for line in tsvfile:
            tab_pos = re.search("\d", line)
            if tab_pos:
                genes = line[tab_pos.start() + 1:].strip('\n').split('/')
                tsvdict[line[0:tab_pos.start()].strip()] = genes

    breast_genes = [int(x) for d in tsvdict for x in tsvdict[d] if d.startswith("breast")]
    return breast_genes
