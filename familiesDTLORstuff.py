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
        #event_mpr = new_DTLOR_DP.build_event_median_graph(MPR)
        # Get the median graph
        #node_median_graph = new_DTLOR_DP.build_node_median_graph(G)
        # Here is how to get the events for one MPR from the graph
        event_median_graph = new_DTLOR_DP.build_event_median_graph(G)
        event_mpr = new_DTLOR_DP.find_MPR(event_median_graph) # Note: use rand=True to get a random median MPR
        #print(new_DTLOR_DP.get_events(event_mpr))
        # Ensure that they each have all location mappings
        #mpr_maps = [node for node in MPR if node[0] is new_DTLOR_DP.NodeType.LOCATION_MAPPING]
        #event_maps = [node for node in event_mpr if node[0] is new_DTLOR_DP.NodeType.LOCATION_MAPPING]
        #assert len(mpr_maps) == len(event_maps)
        events = new_DTLOR_DP.get_events(event_mpr)
        event_cost = new_DTLOR_DP.score_events(events, D, T, L, O, R)
        assert cost == event_cost, "Alg: {}, Events: {}".format(cost, event_cost)
        print("Min Cost: {}".format(cost))
        print("MPRS: {}".format(nmprs))
        print()
        if cost < best_score: 
            # If the score is better than current best
            # Update best score, clear record and add new record
            best_score=cost
            bestMPRs=[(geneRtreeO, MPR)]
        elif cost==best_score:
            bestMPRs.append((geneRtreeO, MPR))

    #sample one MPR from the MPRs for this specific unrooted tree
    optRootedGeneTree, optMPR = random.choice(bestMPRs) 
    return optRootedGeneTree,optMPR
