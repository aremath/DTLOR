import trees, DTLOR_DP, random
import new_DTLOR_DP

def reduceLocusMap(geneUtreeO,locusMapD):
    '''Create a new locus map D with only entries for genes in geneUtreeO.'''
    gtLocusMapD={}
    for leaf in geneUtreeO.leaves():
        # the leaf is a gene number
        gtLocusMapD[leaf] = locusMapD[leaf]
    return gtLocusMapD

def reconcile(argT):
    '''Reconcile a single gene tree.'''

    speciesRtreeO,geneUtreeO,tipMapD,gtLocusMapD,D,T,L,O,R = argT
    
    # convert species tree to dp format
    speciesTreeD = speciesRtreeO.createDtlorD(True)

    # try all rootings and record the ones and their
    # solutions with the best scores
    best_score=float('inf')
    bestMPRs=[]

    for geneRtreeO in geneUtreeO.iterAllRootedTrees():

        geneTreeD = geneRtreeO.createDtlorD(False) # put in dp format
        #MPR,cost=DTLOR_DP.DP(speciesTreeD, geneTreeD, tipMapD, gtLocusMapD, D, T, L, O, R)
        cost, G = new_DTLOR_DP.compute_dtlor_graph(speciesTreeD, geneTreeD, tipMapD, gtLocusMapD, D, T, L, O, R)
        #old_cost, old_G, old_nmprs = DTLOR_DP.DP(speciesTreeD, geneTreeD, tipMapD, gtLocusMapD, D, T, L, O, R)
        nmprs = new_DTLOR_DP.count_MPRs(G)[(new_DTLOR_DP.NodeType.ROOT,)]
        MPR = new_DTLOR_DP.find_MPR(G)
        node_median_graph = new_DTLOR_DP.build_node_median_graph(G)
        event_median_graph = new_DTLOR_DP.build_event_median_graph(G)
        print("Min Cost: {}, {}".format(cost, nmprs))
        #print("Old Min Cost: {}, {}".format(old_cost, nmprs))
        if cost < best_score: 
            # If the score is better than current best
            # Update best score, clear record and add new record
            best_score=cost
            bestMPRs=[(geneRtreeO, MPR)]
        elif cost==best_score:
            bestMPRs.append((geneRtreeO,MPR))

    #sample one MPR from the MPRs for this specific unrooted tree
    optRootedGeneTree, optMPR = random.choice(bestMPRs) 
    return optRootedGeneTree,optMPR
