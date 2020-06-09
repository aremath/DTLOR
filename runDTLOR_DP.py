import sys,copy
import trees,familiesDTLORstuff

# costs
D = 0.3 # duplication
T = 0.4 # transfer
L = 0.4 # loss
O = 0.1 # origin
R = 0.2 # rearrangment

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
            gene=int(L[0])
            if L[1].isdigit():
                value=int(L[1])
            else:
                value = L[1]
            D[gene]=value
    return D

## main

if __name__ == "__main__":

    speciesTreeFN = sys.argv[1]
    geneTreeFN = sys.argv[2]

    # load stuff
    speciesTree = trees.readTree(speciesTreeFN)
    geneTree = trees.loadOneGeneTree(geneTreeFN)

    bigTipMapD = loadD("tipMap.tsv")
    tipMapD = {} # cut down to those in this gene tree
    for leaf in trees.leafList(geneTree):
        tipMapD[leaf]=bigTipMapD[leaf]
    
    locusMapD = loadD("locusMap.tsv")
    gtLocusMapD = familiesDTLORstuff.reduceLocusMap(geneTree,locusMapD)
    locusMapForRootingD = trees.createLocusMapForRootingD(geneTree,copy.deepcopy(gtLocusMapD))
    
    argT = (speciesTree,geneTree,tipMapD,gtLocusMapD,locusMapForRootingD,D,T,L,O,R)

    optRootedGeneTree,optMPR = familiesDTLORstuff.reconcile(argT)

    print("Rooted tree:")
    print(optRootedGeneTree)
    print()
    print("MPR:")
    print(optMPR)
