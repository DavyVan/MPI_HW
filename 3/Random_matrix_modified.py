import random
import sys

## Generates a random dense matrix
## Use: python randomMatrix.py M N

if (len(sys.argv) > 1):

    N = int(sys.argv[1])

    f = open('input.dat','w');


    table= [ [ 0 for i in range(N) ] for j in range(N) ]
    b = [0 for i in range(N)]

    for i in range(N):
        for j in range(N):
            table[i][j] = random.random()
            f.write(str(table[i][j]))
            f.write('\t')
        f.write('\n')

    for i in range(N):
        b[i] = random.random()
        f.write(str(b[i]))
        f.write('\t')
    #f.write('\n')



    f.close()
