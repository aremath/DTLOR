from DTLOR_DP import check_tip, preorder, postorder, delta_r

Infinity = float("inf")

# key: parasite node, value: set of syntenic locations in the clade below that node
def get_synteny_clades(ep, parasite_tree, locus_map, clades):
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
    synteny_clades = get_synteny_clades(parasite_root, parasite_tree, locus_map, {})
    # First, compute C
    C = DTL_reconcile(host_tree, parasite_tree, phi, D, T, L)
    C_g = {g: min([C[(g, s)] for s in preorder(host_tree)]) for g in preorder(parasite_tree)}
    S = synteny_reconcile(host_tree, parasite_tree, locus_map, R, synteny_clades)
    S_g = {g: min([S[(g, l)] for l in synteny_clades[g]]) for g in preorder(parasite_tree)}
    Origin = {}
    Null = {}
    for ep in postorder(parasite_tree):
        _, vp, ep1, ep2 = parasite_tree[ep]
        Origin[ep] = C_g[ep] + S_g[ep] + O
        if check_tip(vp, ep1, ep2):
            Null[ep] = Infinity
        else:
            none_null = Null[ep1] + Null[ep2]
            left_null = Null[ep1] + Origin[ep2]
            right_null = Origin[ep1] + Null[ep2]
            neither_null = Origin[ep1] + Origin[ep2]
            Null[ep] = min([none_null, left_null, right_null, neither_null])
    return min(Null[parasite_root], Origin[parasite_root])

def DTL_reconcile(host_tree, parasite_tree, phi, D, T, L):
    A = {}
    C = {}
    O = {}
    best_switch = {}
    host_root = next(iter(host_tree))
    for ep in postorder(parasite_tree):
        _, vp, ep1, ep2 = parasite_tree[ep]
        vp_is_a_tip = check_tip(vp, ep1, ep2)
        for eh in postorder(host_tree):
            _, vh, eh1, eh2 = host_tree[eh]
            vh_is_a_tip = check_tip(vh, eh1, eh2)
            # Compute A
            if vh_is_a_tip:
                if vp_is_a_tip and phi[ep] == eh:
                    A[(ep, eh)] = 0
                else:
                    A[(ep, eh)] = Infinity
            else:
                # Cospeciation
                if not vp_is_a_tip:
                    cospeciation = min(C[(ep1, eh1)] + C[(ep2, eh2)],
                             C[(ep1, eh2)] + C[(ep2, eh1)])
                else:
                    cospeciation = Infinity
                # Loss
                loss = L + min(C[(ep, eh1)], C[(ep, eh2)])
                A[(ep, eh)] = min(cospeciation, loss)
            # Compute C
            # Duplication
            if not vp_is_a_tip:
                duplication = D + C[(ep1, eh)] + C[(ep2, eh)]
            else:
                duplication = Infinity
            # Transfer
            if not vp_is_a_tip:
                transfer = T + min(C[(ep1, eh)] + best_switch[(ep2, eh)],
                                   C[(ep2, eh)] + best_switch[(ep1, eh)])
            else:
                transfer = Infinity
            C[(ep, eh)] = min(A[(ep, eh)], duplication, transfer)
            # Compute O
            if vh_is_a_tip:
                O[(ep, eh)] = C[(ep, eh)]
            else:
                O[(ep, eh)] = min(C[(ep, eh)], O[(ep, eh1)], O[(ep, eh2)])
        # Compute best_switch
        best_switch[(ep, host_root)] = Infinity
        for eh in preorder(host_tree):
            _, vh, eh1, eh2 = host_tree[eh]
            vh_is_a_tip = check_tip(vh, eh1, eh2)
            if not vh_is_a_tip:
                best_switch[(ep, eh1)] = min(best_switch[(ep, eh)], O[(ep, eh2)])
                best_switch[(ep, eh2)] = min(best_switch[(ep, eh)], O[(ep, eh1)])
    return C

def synteny_reconcile(host_tree, parasite_tree, locus_map, R, synteny_clades):
    S = {}
    allsynteny = set(locus_map.values())
    # Capture R for convenience
    delta = lambda l1, l2: delta_r(l1, l2, R)
    for ep in postorder(parasite_tree):
        _, vp, ep1, ep2 = parasite_tree[ep]
        for lp in allsynteny:
            if check_tip(vp, ep1, ep2):
                if locus_map[ep] == lp:
                    S[(ep, lp)] = 0
                else:
                    S[(ep, lp)] = Infinity
            else:
                # Optimal synteny for left child
                l = min([delta(lp, l1) + S[(ep1, l1)] for l1 in synteny_clades[ep]])
                # Optimal synteny for right child
                r = min([delta(lp, l2) + S[(ep2, l2)] for l2 in synteny_clades[ep]])
                S[(ep, lp)] = l + r
    return S
