import fileinput
import re


for line in fileinput.input():
    if fileinput.isfirstline():
        continue
    first_digit = re.search("\d[^\s-]", line).start()
    disease_name = line[:first_digit].strip()
    genes = line[first_digit:].strip()
    print(u"%s\t%s" % (disease_name, genes))
