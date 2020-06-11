import trees,DTLOR_DP,random

def reduceLocusMap(geneTree,locusMapD):
    '''Create a new locus map D with only entries for genes in geneTree.'''
    gtLocusMapD={}
    for leaf in trees.leafList(geneTree):
        # the leaf is a gene number
        gtLocusMapD[leaf] = locusMapD[leaf]
    return gtLocusMapD
        
def reconcile(argT):
    '''Reconcile a single gene tree.'''

    speciesTree,geneTree,tipMapD,gtLocusMapD,locusMapForRootingD,D,T,L,O,R = argT

    # species tree to right format
    speciesTree=trees.parseTreeForDP(speciesTree,parasite=False)

    # get all possible rootings
    allRootingsL=trees.get_all_rerootings(geneTree, locusMapForRootingD)
    if allRootingsL==[]:  #all rerooting not valid (all nodes have the same loc)
        allRootingsL=[geneTree]

    # try all rootings and record the ones and their
    # solutions with the best scores
    best_score=float('inf')
    bestMPRs=[]
    for geneTree in allRootingsL:
        geneTreeD=trees.parseTreeForDP(geneTree,parasite=True) # gene tree to right format
        MPR,cost=DTLOR_DP.DP(speciesTree, geneTreeD, tipMapD, gtLocusMapD, D, T, L, O, R)
        print("Min Cost: {}".format(cost))
        if cost<best_score: 
            #if the score is better than current best
            #update best score, clear record and add new record
            best_score=cost
            bestMPRs=[]
            bestMPRs.append((geneTree,MPR))
        elif cost==best_score:
            bestMPRs.append((geneTree,MPR))

    #sample one MPR from the MPRs for this specific unrooted tree
    optRootedGeneTree,optMPR=random.choice(bestMPRs) 

    return optRootedGeneTree,optMPR    
