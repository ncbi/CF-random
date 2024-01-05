#! /usr/local/bin/python3

import sys
import numpy as np

if __name__ == '__main__':

    f = open(sys.argv[1]).read().splitlines()

    plddt = [float(x.split()[-2]) for x in f if len(x.split()) > 3 and \
             x.split()[2] == 'CA']

    plddt = np.array(plddt)

    print(sys.argv[1], np.average(plddt))
    #print(plddt)
