import sys

#function to write the double line of the rel adapt from online onto single lines
def to_write(lis): 
    out = lis[0]
    for i in range(1, len(lis)): 
        out +=  '\t' + lis[i]
    return out + '\n' 
    
for line in sys.stdin: 
    s = line.strip().split() 
    sys.stdout.write(to_write(s[:7]))
    sys.stdout.write(to_write(s[7:]))
    
    

    