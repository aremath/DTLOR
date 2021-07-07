import sys,copy
import trees,familiesDTLORstuff
from Tree import *
from multiprocessing import Pool

import tblib.pickling_support
tblib.pickling_support.install()
class ExceptionWrapper(object):

    def __init__(self, ee):
        self.ee = ee
        _, _, self.tb = sys.exc_info()

    def re_raise(self):
        raise self.ee.with_traceback(self.tb)

# costs
D = 1 # duplication
T = 1 # transfer
L = 1 # loss
O = 1 # origin
R = 1 # rearrangment

## funcs

def loadD(fn):
    ''''''
    D={}
    with open(fn,"r") as tempf:
        while True:
            s=tempf.readline()
            if s=='':
                break
            L=s.rstrip().split("\t")
            gene=L[0]
            if L[1].isdigit():
                value=int(L[1])
            else:
                value = L[1]
            D[gene]=value
    return D

## main
#TODO: proper CLI

def recon_wrap(args):
    try:
        return familiesDTLORstuff.reconcile(args)
    except Exception as e:
        return ExceptionWrapper(e)

if __name__ == "__main__":

    speciesTreeFN = sys.argv[1]
    geneTreeFN = sys.argv[2]

    # load stuff
    speciesRtreeO = Rtree()
    speciesRtreeO.fromNewickFileLoadSpeciesTree(speciesTreeFN)
    geneUtreeO = Utree()
    geneUtreeO.fromNewickFile(geneTreeFN)

    bigTipMapD = loadD("tipMap.tsv")
    tipMapD = {} # cut down to those in this gene tree
    for leaf in geneUtreeO.leaves():
        tipMapD[leaf]=bigTipMapD[leaf]
    
    locusMapD = loadD("locusMap.tsv")
    gtLocusMapD = familiesDTLORstuff.reduceLocusMap(geneUtreeO,locusMapD)
    
    argT = (speciesRtreeO,geneUtreeO,tipMapD,gtLocusMapD,D,T,L,O,R)

    argumentL = [argT]

    with Pool(processes=2) as p:
        for result in p.imap_unordered(recon_wrap, argumentL):
            if isinstance(result, ExceptionWrapper):
                result.re_raise()
            optGeneRTreeO,optG = result
    #optRootedGeneTree,optMPR = familiesDTLORstuff.reconcile(argT)

            print("Rooted tree:")
            print(optGeneRTreeO)
            print()
            print("G:")
            print(optG)
