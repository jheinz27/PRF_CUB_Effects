import sys
cor = 0
tot = 0 

found = 0
tot = 0

found_prfs = 0
all_prfs = 0
tot_na = 0 

for line in sys.stdin: 
    if not line.startswith('#'):
        regs = line.strip().split('\t')[-2].split(',')

        starts = []
        ends = []

        if regs[0] == 'n/a':
            regs = ['0-0']
            tot_na += 1
        for r in regs:
            sp = r.split('-')
            starts.append(int(sp[0]))
            ends.append(int(sp[1]))
        prfs = line.strip().split('\t')[-1].split(',')


        #delete out prf as soon as found to avoid overlapping regions conflict

        for i in range(len(starts)):
        
            tot += 1 
            toremove = []
            for p in prfs: 
                if starts[i] <= int(p) <= ends[i]: 
                    found += 1
                    toremove.append(p)
            for temp in toremove:
                prfs.remove(temp)
                
        prfs = line.strip().split('\t')[-1].split(',')
   
        for prf in prfs:
            all_prfs += 1 
            for i in range(len(starts)):
                
                if starts[i] <= int(prf) <= ends[i]: 
                    found_prfs += 1
                 
                    break 
                    
        
                

print('Prob of regions: ' + str(found) + '/' + str(tot - tot_na) + ' prop= ' + str(round(found/(tot - tot_na), 6)))

print('Filtered Prfs: ' + str(found_prfs) + '/' + str(all_prfs) + ' prop= ' + str(round(found_prfs/all_prfs, 6)))                      
            
                    
           

        