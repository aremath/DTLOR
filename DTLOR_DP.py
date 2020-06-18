# DTLOR_DP.py
# Nuo Liu, HMC 2019-2020 Thesis project
# Modified from Ran Libeskind-Hadas(June 2015)

# The basic DP algorithm for reconciling pairs of trees
# Altered and expanded by Carter Slocum and Annalise Schweickart
# A tree is represented as a dictionary of key-value pairs where a key is an
# edge name and the value is a tuple of the form
# (start vertex, end vertex, left child edge name, right child edge name)
# An edge name may be None.  The "dummy" edge leading to the root of the
# parasite tree, denoted e^P in the technical report, must be named "pTop".
# Edited by Annalise Schweickart and Carter Slocum, July 2015 to return
# the DTL reconciliation graph that uses frequency scoring, as well as the
# number of reconciliations of the host and parasite trees
# Edited by Ross Mawhorter in 2020 to reduce the running time per Rose Liu's formulation.

# import newickFormatReader
from Greedy import *
import copy
import sys, glob, os
import time
from random import choice

def valid_star(is_top_star, allsynteny):
    """
    returns the set of valid syntenies for descendants
    """
    if is_top_star:
        new_locs = copy.deepcopy(allsynteny)
        new_locs.append("*")
        return new_locs
    else:
        return allsynteny

def valid(synteny, allsynteny):
    """
    returns the set of valid syntenies for the descendents
    """
    is_top_star = (synteny == "*")
    return valid_star(is_top_star, allsynteny)

def delta_cost(synteny1, synteny2, O, R):
    if synteny1==synteny2:
        return 0
    elif synteny1=="*":
        return O
    else:
        return R

Infinity = float('inf')

def preorder(tree, rootEdgeName):
    """ Takes a tree as input (see format description above) and returns a 
    list of the edges in that tree in preorder (high edges to low edges)"""

    value = tree[rootEdgeName]
    _,_,leftChildEdgeName,rightChildEdgeName = value

    # base case
    if leftChildEdgeName == None: # then rightChildEdgeName == None also
        return [rootEdgeName]
    else:
        return [rootEdgeName] + \
                preorder(tree, leftChildEdgeName) + \
                preorder(tree, rightChildEdgeName)

def postorder(tree, rootEdgeName):
    """ Takes a tree as input (see format description above) and returns a 
    list of the edges in that tree in postorder (low edges to high edges)"""

    value = tree[rootEdgeName]
    _,_,leftChildEdgeName,rightChildEdgeName = value
    # base case
    if leftChildEdgeName == None: # then rightChildEdgeName == None also
        return [rootEdgeName]
    else:
        return postorder(tree, leftChildEdgeName) + \
               postorder(tree, rightChildEdgeName) + \
               [rootEdgeName]

# returns ce = oplus_i ce_i
#TODO: rename
def combine_costs(events_list):
    """
    Combines the events list into a single events list
    Each element of events_list is a tuple of (cost, [events])
    Which holds the DP entry for some part of the DP and the associated events
    This combines them into the entry of lowest cost, combining
    all events that have lowest cost
    """
    cost = Infinity
    events = []
    for c,e in events_list:
        print("\tCost: {}".format(c))
        print("\tEvents: {}".format(events))
        assert type(e) is list, e
        assert (type(c) is float) or (type(c) is int), (c, e)
        #print(c, e)
        if cost > c:
            cost = c
            events = e
        elif cost == c:
            events.extend(e)
    print("Events: {}".format(events))
    return (cost, events)

#TODO: optimization where valid returns only the synteny locations for the clade UNDER a given gene node
#TODO: delta no longer needs to use * since null counts origins separately...?
#TODO: iterable rather than actually creating lists
#TODO: Refactor other code to reflect the tree representation change:
# ("C", None, None) versus ("C", (None, None, None, None), (None, None, None, None))
# and ("L", (...), None) versus ("L", (...), (None, None, None, None)

def DP(hostTree, parasiteTree, phi, locus_map, D, T, L, Origin, R):
    """ Takes a hostTree, parasiteTree, tip mapping function phi, a locus_map, 
        and duplication cost (D), transfer cost (T), loss cost (L), 
        origin cost (O) and rearrange cost(R) and returns the DTLOR graph in the form of a dictionary, 
        as well as a the number of maximum parsimony reconciliations. The notation and 
        dynamic programming algorithm are explained in the tech report.
        Cospeciation is assumed to cost 0. """

    A = {}  # A, C, O, and bestSwitch are all defined in tech report
    C = {}
    O = {}
    eventsDict = {} # Dictionary to keep track of events that correspond to the min cost reconciliation 
    bestSwitch = {} 
    allsynteny=list(locus_map.values())
    # Capture O, R for ease of use
    delta = lambda s1, s2: delta_cost(s1, s2, Origin, R)
    #print("The dimensions is %d by %d by %d by %d"%(len(postorder(parasiteTree, "pTop")),len(Allsynteny), len(Allsynteny),len(postorder(hostTree, "hTop"))))
    for ep in postorder(parasiteTree, "pTop"):
        _,vp,ep1,ep2 = parasiteTree[ep]
        # is vp a tip?
        if ep1 == None: # then ep2 == None too and vp is a tip!
            vpIsATip = True
        else:
            vpIsATip = False
        #TODO: rename l_bottom
        for l_bottom in allsynteny:  #for start and end vertex of gene edge
            for eh in postorder(hostTree, "hTop"):
                _,vh,eh1,eh2 = hostTree[eh]
                eventsDict[(vp, vh, l_bottom)] = []
                # is vh a tip?
                if eh1 == None: # then eh2 == None too and vh is a tip!
                    vhIsATip = True
                else:
                    vhIsATip = False
                # Compute A[(ep, eh, l_bottom)]
                if vhIsATip:
                    if vpIsATip and phi[vp] == vh and locus_map[vp]==l_bottom:
                        A[(ep, eh, l_bottom)] = (0, [("C", None, None)])
                    else: 
                        A[(ep, eh, l_bottom)] = (Infinity, [])
                else: # vh is not a tip
                    # Compute cospeciation events
                    if not vpIsATip:
                        def get_cospeciations(l1, l2):
                            synteny_cost = delta(l_bottom, l1) + delta(l_bottom, l2)
                            co1 = (synteny_cost + C[(ep1, eh1, l1)][0] + C[(ep2, eh2, l2)][0], \
                                    [("S", (ep1, eh1, l1), (ep2, eh2, l2))])
                            co2 = (synteny_cost + C[(ep1, eh2, l1)][0] + C[(ep2, eh1, l2)][0], \
                                    [("S", (ep1, eh2, l1), (ep2, eh1, l2))])
                            return combine_costs([co1, co2])
                        cospeciation_list = [get_cospeciations(l1, l2) for l1 in allsynteny for l2 in allsynteny]
                        cospeciations = combine_costs(cospeciation_list)
                    else:
                        cospeciations = (Infinity, [])
                    # Compute loss events
                    # eh1 is the branch where ep is lost
                    loss_eh1 = (C[(ep, eh2, l_bottom)][0] + L, [("L", (vp, eh2, l_bottom), None)])
                    # eh2 is the branch where ep is lost
                    loss_eh2 = (C[(ep, eh1, l_bottom)][0] + L, [("L", (vp, eh1, l_bottom), None)])
                    losses = combine_costs([loss_eh1, loss_eh2])

                    # Determine which event occurs for A[(ep, eh, l_bottom)] 
                    A[(ep, eh, l_bottom)] = combine_costs([cospeciations, losses])

                # Compute C[(ep, eh,l_top, l_bottom)]
                # First, compute duplications
                if not vpIsATip:
                    DUPepeh=Infinity
                    # List to keep track of lowest cost duplication event
                    dupList=[Infinity]
                    def get_duplication(l1, l2):
                        synteny_cost = delta(l_bottom, l1) + delta(l_bottom, l2)
                        dup_cost = C[(ep1, eh, l1)][0] + C[(ep2, eh, l2)][0] + D
                        dup_event = ("D", (ep1, vh, l1), (ep2, vh, l2))
                        return (dup_cost, [dup_event])
                    dup_list = [get_duplication(l1, l2) for l1 in allsynteny for l2 in allsynteny]
                    duplications = combine_costs(dup_list)
                else:
                    duplications = (Infinity, [])
               
                # Compute transfer table
                if not vpIsATip:
                    switchList = [] # List to keep track of lowest cost switch
                    SWITCHepeh=Infinity
                    #need to find all possible children syntenies
                    def get_transfer(l1, l2):
                        synteny_cost = delta(l_bottom, l1) + delta(l_bottom, l2)
                        # Cost to transfer ep2
                        ep2_cost, ep2_locations = bestSwitch[(ep2, eh, l2)]
                        ep2_switch_cost = T + synteny_cost + C[(ep1, eh, l1)][0] + ep2_cost
                        ep2_switch_events = [("T", (ep1, vh, l1), (ep2, location[1], l2)) for location in \
                                ep2_locations]
                        ep2_switch = (ep2_switch_cost, ep2_switch_events)
                        # Cost to transfer ep1
                        ep1_cost, ep1_locations = bestSwitch[(ep1, eh, l1)]
                        ep1_switch_cost = T + synteny_cost + C[(ep2, eh, l2)][0] + ep1_cost
                        ep1_switch_events = [("T", (ep2, vh, l2), (ep1, location[1], l1)) for location in \
                                ep1_locations]
                        ep1_switch = (ep1_switch_cost, ep_1_switch_events)
                        return combine_costs([ep2_switch, ep1_switch])
                    transfer_list = [get_transfer(l1,l2) for l1 in all_synteny for l2 in all_synteny]
                    transfers = combine_costs(transfer_list)
                else:
                    transfers = (Infinity, [])

                # Compute C[(ep, eh, l_top, l_bottom)] and add the event or events with that cost
                # to the dictionary eventsDict
                C[(ep, eh, l_bottom)] = \
                        combine_costs([A[(ep, eh, l_bottom)], duplications, transfers])

                # Compute O[(ep, eh, l_bottom)]
                if vhIsATip: 
                    O[(ep, eh, l_bottom)] = (C[(ep, eh, l_bottom)][0], [(vp, vh, l_bottom)])
                else: 
                    o_c = (C[(ep, eh, l_bottom)][0], [(vp, vh, l_bottom)])
                    o_l = (O[(ep, eh1, l_bottom)])
                    o_r = (O[(ep, eh2, l_bottom)])
                    O[(ep, eh, l_bottom)] = combine_costs([o_c, o_l, o_r])

            # Compute bestSwitch values
            bestSwitch[(ep, "hTop", l_bottom)] = (Infinity, [])
            for eh in preorder(hostTree, "hTop"):
                _, vh, eh1, eh2 = hostTree[eh]

                # is vh a tip?
                if eh1 == None: # then eh2 == None too and vh is a tip!
                    vhIsATip = True
                else:
                    vhIsATip = False
                # Find the best switches and switch locations
                if eh1 != None and eh2 != None: # not a tip
                    ep_bestSwitch = bestSwitch[(ep, eh, l_bottom)]
                    O_eh2=O[(ep, eh2, l_bottom)]
                    O_eh1=O[(ep, eh1, l_bottom)]
                    bestSwitch[(ep, eh1, l_bottom)] = combine_costs([ep_bestSwitch, O_eh2, O_eh1])

        #TODO: what if ep1 or ep2 is None? i.e. a tip
        def get_single_null(eh, l):
            # Left child stays null
            l_map = (ep1, "hTop", "*")
            r_map = (ep2, eh, l)
            left_null_cost = C[l_map][0] + C[r_map][0] + Origin
            left_null_event = ("N", l_map, r_map)
            left_null = (left_null_cost, [left_null_event])
            # Right child stays null
            l_map = (ep1, eh, l)
            r_map = (ep2, "hTop", "*")
            right_null_cost = C[l_map][0] + C[r_map][0] + Origin
            right_null_event = ("N", l_map, r_map)
            right_null = (right_null_cost, [right_null_event])
            return combine_costs([left_null, right_null])
        single_null_list = [get_single_null(eh, l) \
                for eh in postorder(hostTree, "hTop") for l in allsynteny]
        single_null = combine_costs(single_child_list)

        # Neither child gets a synteny
        l_map = (ep1, "hTop", "*")
        r_map = (ep2, "hTop", "*")
        both_null_cost = C[l][0] + C[r][0]
        both_null_event = ("N", l_map, r_map)
        both_null = (both_null_cost, [both_null_event])

        def get_neither_child_null(eh1, l1, eh2, l2):
            l_map = (ep1, eh1, l1)
            r_map = (ep2, eh2, l2)
            neither_null_cost = C[l_map][0] + C[r_map][0] + 2 * Origin
            neither_null_event = ("N", l_map, r_map)
            return (neither_null_cost, neither_null_event)

        neither_null_list = [get_neither_child_null(eh1, l1, eh2, l2) \
                for eh1 in postorder(hostTree, "hTop") for eh2 in postorder(hostTree, "hTop") \
                for l1 in allsynteny for l2 in allsynteny]
        neither_null = combine_costs(neither_null_list)

        C[(g, "hTop", "*")] = combine_costs([single_null, both_null, neither_null])

    # Cost for assigning the root
    def get_root(eh, l):
        m = C[("pTop", eh, l)][0]
        return (m[0] + Origin, m[1])
    root_list = [get_root_not_null(eh, l) for eh in postorder(hostTree, "hTop") for l in allsynteny + ["*"]]
    min_cost, _ = combine_costs(root_not_null_list)

    # Find the mapping nodes involving pTop of minimum cost
    best_roots = [m for m,c in C.items if m[0] == "pTop" and c[0] == min_cost]

    # This picks a random MPR from the optimal ones
    MPR = find_MPR(best_roots, C)
    return MPR, min_cost

def find_MPR(best_roots, C):
    """
    Find a single MPR for the given C dict and best_roots.
    """
    return find_MPR_helper(best_roots, C, {})

def find_MPR_helper(nodes, C, MPR):
    """
    Recursively find a single MPR using C. Does the work for find_MPR.
    """
    mapping = choice(nodes)
    event = choice(C[mapping][1])
    MPR[mapping] = event
    e_type, e_left, e_right = event
    if e_left is not None:
        _ = findOneMPR([e_left], C, MPR)
    if e_right is not None:
        _ = findOneMPR([e_right], C, MPR)
    return MPR

def MPR_graph(best_roots, C):
    """
    Find the MPR graph for a given C dict and best_roots.
    """
    return MPR_graph_helper(best_roots, C, {})

def MPR_graph(nodes, C, G):
    """
    Recursively create the entire MPR graph. Does the work for MPR_graph.
    """
    for mapping in nodes:
        events = c[mapping][1]
        G[mapping] = events
        for e_type, e_left, e_right in events:
            MPR_graph([e_left], C, G)
            MPR_graph([e_right], C, G)
    return G

def preorderDTLORsort(DTLOR, ParasiteRoot):
    """This takes in a DTL reconciliation graph and parasite root and returns 
    a sorted list, orderedKeysL, that is ordered by level from largest to 
    smallest, where level 0 is the root and the highest level has tips."""

    keysL = orderDTLOR(DTLOR, ParasiteRoot)
    uniqueKeysL = sortHelper(DTLOR, keysL)
    orderedKeysL = []
    levelCounter = 0
    while len(orderedKeysL) < len(keysL):
        for mapping in keysL:
            if mapping[-1] == levelCounter:
                orderedKeysL = orderedKeysL + [mapping]
        levelCounter += 1
    
    # lastLevel = orderedKeysL[-1][1]
    return orderedKeysL

def addScores(treeMin, DTLORDict, ScoreDict):
    """Takes the list of reconciliation roots, the DTLOR reconciliation graph, 
    a dictionary of parent nodes, and a dictionary of score values, and 
    returns the DTLOR with the normalized frequency scores calculated."""
    newDTLOR = copy.deepcopy(DTLORDict)
    parentsDict = {}
    ParasiteRoot=treeMin[0][0]
    preOrder = preorderDTLORsort(DTLORDict, ParasiteRoot)
    for root in preOrder:
        #TODO: when would this be true?
        if root != (None, None): #format (key, level)
            vertices = root[0] #format (vp, vh, l, l)
            if root[1] == 0:
                parentsDict[vertices] = ScoreDict[vertices]
            for n in range(len(DTLORDict[vertices])-1): #the last item is always the cost
                _,child1,child2,oldScore = DTLORDict[vertices][n]
                newDTLOR[vertices][n][3] = parentsDict[vertices] * \
                (1.0 * oldScore / ScoreDict[vertices])
                if child1 is not None:
                    if child1 in parentsDict:
                        parentsDict[child1] += newDTLOR[vertices][n][3]
                    else: 
                        parentsDict[child1] = newDTLOR[vertices][n][3] 
                if child2 is not None:
                    if child2 in parentsDict:
                        parentsDict[child2] += newDTLOR[vertices][n][3]
                    else: 
                        parentsDict[child2] = newDTLOR[vertices][n][3]
    normalize = newDTLOR[preOrder[-1][0]][0][-1]  #score of the top most event
    for key in newDTLOR:
        for event in newDTLOR[key][:-1]:
            event[-1] = event[-1]/normalize
    
    return newDTLOR, normalize

