import sys
all_exons = [] 
geneline = ''

for line in sys.stdin: 
    if not line.startswith('#'):
        s = line.strip().split('\t')
        if s[2] == 'gene': 
            if len(all_exons) == 1: 
                sys.stdout.write(geneline)
            geneline = line
            all_exons = []
            
        if s[2] == "exon":
            all_exons.append(line)
        

if len(all_exons) == 1: 
        sys.stdout.write(geneline)
        geneline = line
               
    