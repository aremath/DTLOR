import sys,copy
import argparse
import trees,familiesDTLORstuff
from Tree import *

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

# Main
if __name__ == "__main__":
    parser = argparse.ArgumentParser("Maximum Parsimony DTLOR reconciliation")
    parser.add_argument("species", metavar="<species_tree>", help="file path for species tree")
    parser.add_argument("gene", metavar="<gene_tree>", help="file path for gene tree")
    parser.add_argument("--mapping", metavar="<mapping_file>", default="tipMap.tsv",
                        required=False, help="file path for tip mapping")
    parser.add_argument("--locus", metavar="<locus_file>", default="locusMap.tsv",
                        required=False, help="file path for locus mapping")
    parser.add_argument("-d", type=int, metavar="<duplication_cost>", default=1,
                        help="cost incurred by duplication")
    parser.add_argument("-t", type=int, metavar="<duplication_cost>", default=1,
                        help="cost incurred by transfer")
    parser.add_argument("-l", type=int, metavar="<duplication_cost>", default=1,
                        help="cost incurred by loss")
    parser.add_argument("-o", type=int, metavar="<duplication_cost>", default=1,
                        help="cost incurred by origin")
    parser.add_argument("-r", type=int, metavar="<duplication_cost>", default=1,
                        help="cost incurred by rearrangement")
    args = parser.parse_args()

    speciesTreeFN = args.species
    geneTreeFN = args.gene

    # load stuff
    speciesTree = Rtree()
    speciesTree.fromNewickFileLoadSpeciesTree(speciesTreeFN)
    geneTree = Utree()
    geneTree.fromNewickFile(geneTreeFN)

    bigTipMapD = loadD(args.mapping)
    tipMapD = {} # cut down to those in this gene tree
    for leaf in geneTree.leaves():
        tipMapD[leaf]=bigTipMapD[leaf]
    
    locusMapD = loadD(args.locus)
    gtLocusMapD = familiesDTLORstuff.reduceLocusMap(geneTree,locusMapD)
    
    argT = (speciesTree,geneTree,tipMapD,gtLocusMapD,args.d,args.t,args.l,args.o,args.r)

    optRootedGeneTree,optMPR = familiesDTLORstuff.reconcile(argT)

    print("Rooted tree:")
    print(optRootedGeneTree)
    print()
    print("MPR:")
    print(optMPR)
