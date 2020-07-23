import sys,copy
import trees,familiesDTLORstuff
from Tree import *

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

if __name__ == "__main__":

    speciesTreeFN = sys.argv[1]
    geneTreeFN = sys.argv[2]

    # load stuff
    speciesTree = Rtree()
    speciesTree.fromNewickFileLoadSpeciesTree(speciesTreeFN)
    geneTree = Utree()
    geneTree.fromNewickFile(geneTreeFN)

    bigTipMapD = loadD("tipMap.tsv")
    tipMapD = {} # cut down to those in this gene tree
    for leaf in geneTree.leaves():
        tipMapD[leaf]=bigTipMapD[leaf]
    
    locusMapD = loadD("locusMap.tsv")
    gtLocusMapD = familiesDTLORstuff.reduceLocusMap(geneTree,locusMapD)
    
    argT = (speciesTree,geneTree,tipMapD,gtLocusMapD,D,T,L,O,R)

    optRootedGeneTree,optMPR = familiesDTLORstuff.reconcile(argT)

    print("Rooted tree:")
    print(optRootedGeneTree)
    print()
    print("MPR:")
    print(optMPR)
