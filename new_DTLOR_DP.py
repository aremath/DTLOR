from collections import defaultdict
from enum import Enum, auto
from itertools import product
from DTLOR_DP import check_tip, preorder, postorder, delta_r, find_min_events

Infinity = float("inf")

class GraphType(Enum):
    """
    Defines how each NodeType interacts with the graph when
    choosing a reconciliation
    """
    CHOOSE = auto()
    ALL = auto()

class NodeType(Enum):
    """
    Encodes the type of a node in a reconciliation graph
    """
    # DTL events
    COSPECIATION = (auto(), GraphType.ALL)
    DUPLICATION = (auto(), GraphType.ALL)
    LOSS = (auto(), GraphType.ALL)
    TRANSFER = (auto(), GraphType.ALL)
    # Assigns locations to the left and right children of a gene node
    LOCATION_ASSIGNMENT = (auto(), GraphType.ALL)
    # Maps a gene tree node to a list of locations
    LOCATION_LIST = (auto(), GraphType.CHOOSE)
    # Maps a gene tree node to a syntenic location
    LOCATION_MAPPING = (auto(), GraphType.CHOOSE)
    # Maps a gene tree node to a list of species tree nodes
    # Used for choosing the optimal root of a reconciliation
    SPECIES_LIST = (auto(), GraphType.CHOOSE)
    # Maps gene tree node to species tree node
    SPECIES_MAPPING = (auto(), GraphType.CHOOSE)
    # When a gene gets a true location
    ORIGIN = (auto(), GraphType.ALL)
    # Should not have children
    CONTEMPORANEOUS = (auto(), None)
    # "Root" event for choosing the syntenic location of the root
    ROOT = (auto(), GraphType.CHOOSE)
    def __init__(self, e_id, graph_type):
        self._e_id = e_id
        self.graph_type = graph_type
    def __repr__(self):
        return "{}".format(self.name)

def get_synteny_clades(ep, parasite_tree, locus_map, clades):
    """
    Get the set of syntenic locations associated with each subtree of the parasite tree.
    Returns a dictionary where the key is a parasite edge and the value is a set of syntenic
    locations in the subtree below that edge.
    """
    _, vp, ep1, ep2 = parasite_tree[ep]
    if check_tip(vp, ep1, ep2):
        clades[ep] = set([locus_map[ep]])
    else:
        get_synteny_clades(ep1, parasite_tree, locus_map, clades)
        get_synteny_clades(ep2, parasite_tree, locus_map, clades)
        clades[ep] = clades[ep1] | clades[ep2]
    return clades

def DP(host_tree, parasite_tree, phi, locus_map, D, T, L, O, R):
    parasite_root = next(iter(parasite_tree))
    host_root = next(iter(host_tree))
    # First, compute C
    C, C_star, C_graph = DTL_reconcile(host_tree, parasite_tree, phi, D, T, L)
    S, S_star, S_graph = synteny_reconcile(host_tree, parasite_tree, locus_map, R)
    # Union the two graphs before adding Null events
    G = {**C_graph, **S_graph}
    Origin = {}
    Null = {}
    for ep in postorder(parasite_tree):
        _, vp, ep1, ep2 = parasite_tree[ep]
        o_cost = C_star[ep][0] + S_star[ep][0] + O
        c_choice = (NodeType.SPECIES_LIST, ep)
        s_choice = (NodeType.LOCATION_LIST, ep)
        o_event = (NodeType.ORIGIN, ep)
        G[o_event] = [c_choice, s_choice]
        Origin[ep] = o_cost, [o_event]
        if check_tip(vp, ep1, ep2):
            Null[ep] = Infinity
        else:
            # Null or Origin for left child
            left_null = Null[ep1], [(NodeType.LOCATION_MAPPING, ep1, "*")]
            left_origin = Origin[ep1]
            left_cost, left_nodes = find_min_events([left_null, left_origin])
            # Null or Origin for right child
            right_null = Null[ep1], [(NodeType.LOCATION_MAPPING, ep2, "*")]
            right_origin = Origin[ep2]
            right_cost, right_nodes = find_min_events([left_null, left_origin])
            a_nodes = []
            for l_node, r_node in product(left_nodes, right_nodes):
                node = (NodeType.LOCATION_ASSIGNMENT, l_node, r_node)
                S_graph[node] = [l_node, r_node]
                a_nodes.append(node)
            cost = left_cost + right_cost
            Null[ep] = cost
            G[(NodeType.LOCATION_MAPPING, ep, "*")] = a_nodes
    # Compute the final choice node and cost: does the root get a location or "*"?
    root_null = (Null[parasite_root], [(NodeType.LOCATION_MAPPING, parasite_root, "*")])
    root_origin = Origin[parasite_root]
    root_cost, root_events = find_min_events([root_null, root_origin])
    G[(NodeType.ROOT,)] = root_events
    G = prune_graph(G)
    return root_cost, G

def DTL_reconcile(host_tree, parasite_tree, phi, D, T, L):
    A = {}
    C = {}
    # node -> [node]
    C_graph = {}
    O = {}
    # gene_node -> (cost, [species_node])
    C_star = {}
    best_switch = {}
    host_root = next(iter(host_tree))
    for ep in postorder(parasite_tree):
        C_star[ep] = (Infinity, [])
        _, vp, ep1, ep2 = parasite_tree[ep]
        vp_is_a_tip = check_tip(vp, ep1, ep2)
        for eh in postorder(host_tree):
            _, vh, eh1, eh2 = host_tree[eh]
            vh_is_a_tip = check_tip(vh, eh1, eh2)
            # Relevant mapping nodes, for convenience
            ep_eh_m = (NodeType.SPECIES_MAPPING, ep, eh)
            ep_eh1_m = (NodeType.SPECIES_MAPPING, ep, eh1)
            ep_eh2_m = (NodeType.SPECIES_MAPPING, ep, eh2)
            ep1_eh_m = (NodeType.SPECIES_MAPPING, ep1, eh)
            ep2_eh_m = (NodeType.SPECIES_MAPPING, ep2, eh)
            # Compute A
            if vh_is_a_tip:
                if vp_is_a_tip and phi[ep] == eh:
                    e = (NodeType.CONTEMPORANEOUS, None, None)
                    A[ep_eh_m] = (0, [e])
                else:
                    A[ep_eh_m] = (Infinity, [])
            else:
                # Cospeciation
                if not vp_is_a_tip:
                    co1_event = (NodeType.COSPECIATION,
                            (NodeType.SPECIES_MAPPING, ep1, eh1),
                            (NodeType.SPECIES_MAPPING, ep2, eh2))
                    co1_cost = C[(NodeType.SPECIES_MAPPING, ep1, eh1)] + \
                            C[(NodeType.SPECIES_MAPPING, ep2, eh2)]
                    co1 = (co1_cost, [co1_event])
                    co2_event = (NodeType.COSPECIATION,
                            (NodeType.SPECIES_MAPPING, ep1, eh2),
                            (NodeType.SPECIES_MAPPING, ep2, eh1))
                    co2_cost = C[(NodeType.SPECIES_MAPPING, ep1, eh2)] + \
                            C[(NodeType.SPECIES_MAPPING, ep2, eh1)]
                    co2 = (co2_cost, [co2_event])
                    cospeciation = find_min_events([co1, co2])
                else:
                    cospeciation = (Infinity, [])
                # Loss
                loss_eh1 = (C[ep_eh2_m] + L, [(NodeType.LOSS,
                    (NodeType.SPECIES_MAPPING, ep, eh2), None)])
                loss_eh2 = (C[ep_eh1_m] + L, [(NodeType.LOSS,
                    (NodeType.SPECIES_MAPPING, ep, eh1), None)])
                loss = find_min_events([loss_eh1, loss_eh2])
                A[ep_eh_m] = find_min_events([cospeciation, loss])
            # Compute C
            # Duplication
            if not vp_is_a_tip:
                dup_event = (NodeType.DUPLICATION,
                        (NodeType.SPECIES_MAPPING, ep1, eh), 
                        (NodeType.SPECIES_MAPPING, ep2, eh))
                duplication = (D + C[ep1_eh_m] + C[ep2_eh_m], [dup_event])
            else:
                duplication = (Infinity, [])
            # Transfer
            if not vp_is_a_tip:
                ep2_cost, ep2_locations = best_switch[ep2_eh_m]
                ep2_switch_cost = T + C[ep1_eh_m] + ep2_cost
                ep2_switch_events = [(NodeType.TRANSFER,
                    (NodeType.SPECIES_MAPPING, ep1, eh),
                    (NodeType.SPECIES_MAPPING, ep2, location)) \
                        for location in ep2_locations]
                ep2_switch = (ep2_switch_cost, ep2_switch_events)
                ep1_cost, ep1_locations = best_switch[ep1_eh_m]
                ep1_switch_cost = T + C[ep2_eh_m] + ep1_cost
                ep1_switch_events = [(NodeType.TRANSFER, 
                    (NodeType.SPECIES_MAPPING, ep2, eh),
                    (NodeType.SPECIES_MAPPING, ep1, location)) \
                        for location in ep1_locations]
                ep1_switch = (ep1_switch_cost, ep1_switch_events)
                transfer = find_min_events([ep2_switch, ep1_switch])
            else:
                transfer = (Infinity, [])
            cost, events = find_min_events([A[ep_eh_m], duplication, transfer])
            C[ep_eh_m] = cost
            C_graph[ep_eh_m] = events
            for event in C_graph[ep_eh_m]:
                t,l,r = event
                event_nodes = []
                if l is not None:
                    event_nodes.append(l)
                if r is not None:
                    event_nodes.append(r)
                C_graph[event] = event_nodes
            C_star[ep] = find_min_events([C_star[ep], (C[ep_eh_m], [eh])])
            # Compute O: (cost, [mapping_node])
            O_c = (C[ep_eh_m], [ep_eh_m])
            if vh_is_a_tip:
                O[ep_eh_m] = O_c
            else:
                O[ep_eh_m] = find_min_events([O_c, O[ep_eh1_m], O[ep_eh2_m]])
        # Compute best_switch
        best_switch[(NodeType.SPECIES_MAPPING, ep, host_root)] = (Infinity, [])
        for eh in preorder(host_tree):
            _, vh, eh1, eh2 = host_tree[eh]
            vh_is_a_tip = check_tip(vh, eh1, eh2)
            ep_eh_m = (NodeType.SPECIES_MAPPING, ep, eh)
            ep_eh1_m = (NodeType.SPECIES_MAPPING, ep, eh1)
            ep_eh2_m = (NodeType.SPECIES_MAPPING, ep, eh2)
            if not vh_is_a_tip:
                best_switch[ep_eh1_m] = find_min_events([best_switch[ep_eh_m], O[ep_eh2_m]])
                best_switch[ep_eh2_m] = find_min_events([best_switch[ep_eh_m], O[ep_eh1_m]])
        # Now that we're done computing C_star[ep], add appropriate choice events
        species_choice = (NodeType.SPECIES_LIST, ep)
        C_graph[species_choice] = [(NodeType.SPECIES_MAPPING, ep, eh) for eh in C_star[ep][1]]
    return C, C_star, C_graph

#TODO gene/species vs. parasite/host
def synteny_reconcile(host_tree, parasite_tree, locus_map, R):
    # gene_node -> [location]
    S_star = {}
    # (g,s) -> cost
    S = {}
    # node -> [[node]]
    S_graph = {}
    allsynteny = set(locus_map.values())
    # Capture R for convenience
    delta = lambda l1, l2: delta_r(l1, l2, R)
    for ep in postorder(parasite_tree):
        S_star[ep] = (Infinity, [])
        _, vp, ep1, ep2 = parasite_tree[ep]
        for lp in allsynteny:
            ep_lp_m = (NodeType.LOCATION_MAPPING, ep, lp)
            if check_tip(vp, ep1, ep2):
                if locus_map[ep] == lp:
                    S[ep_lp_m] = 0
                    S_graph[ep_lp_m] = []
                else:
                    S[(NodeType.LOCATION_MAPPING, ep, lp)] = Infinity
            else:
                # Synteny cost for the left child
                l_keep_m = (NodeType.LOCATION_MAPPING, ep1, lp)
                l_keep = (S[l_keep_m], [l_keep_m])
                l_rearrange = (S_star[ep1][0] + R, [(NodeType.LOCATION_LIST, ep1)])
                l_cost, l_nodes = find_min_events([l_keep, l_rearrange])
                # Synteny cost for the right child
                r_keep_m = (NodeType.LOCATION_MAPPING, ep2, lp)
                r_keep = (S[r_keep_m], [r_keep_m])
                r_rearrange = (S_star[ep2][0] + R, [(NodeType.LOCATION_LIST, ep2)])
                r_cost, r_nodes = find_min_events([r_keep, r_rearrange])
                a_nodes = []
                # Create the appropriate assignment nodes
                for l_node, r_node in product(l_nodes, r_nodes):
                    node = (NodeType.LOCATION_ASSIGNMENT, l_node, r_node)
                    S_graph[node] = [l_node, r_node]
                    a_nodes.append(node)
                cost = l_cost + r_cost
                S[(NodeType.LOCATION_MAPPING, ep, lp)] = cost
                S_graph[ep_lp_m] = a_nodes
                S_star[ep] = find_min_events([S_star[ep], (cost, [lp])])
        # Create the appropriate choice node
        location_choice = (NodeType.LOCATION_LIST, ep)
        S_graph[location_choice] = [(NodeType.LOCATION_MAPPING, ep, lp) for lp in S_star[ep][1]]
    return S, S_star, S_graph

#TODO: methods of Graph class?
def prune_graph(G):
    """
    Removes the unnecessary nodes from G.
    """
    new_G = {}
    extant_nodes = [(NodeType.ROOT,)]
    while len(extant_nodes) != 0:
        node = extant_nodes.pop()
        if node[0].graph_type is not None:
            children = G[node]
            new_G[node] = children
            extant_nodes.extend(children)
    return new_G

def find_MPR(G, MPR, rand=False):
    MPR = {}
    extant_nodes = [(NodeType.ROOT,)]
    while len(extant_nodes) != 0:
        node = extant_nodes.pop()
        node_children = G[node]
        if node[0].graph_type is GraphType.CHOOSE:
            if rand:
                choice = random.choice(node_children)
            else:
                choice = node_children[0]
            MPR[node] = choice
            extant_nodes.append(choice)
        elif node[0].graph_type is GraphType.ALL:
            MPR[node] = node_children
            extant_nodes.extend(node_children)
        elif node[0].graph_type is None:
            pass
        else:
            assert False, "Bad GraphType"
    return MPR
