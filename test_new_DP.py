import trees, familiesDTLORstuff, new_DTLOR_DP

## Main

argT = (933, ('s0', ('s1', ('s2', ('E_coli_ATCC11775', (), (), None), ('E_coli_K12', (), (), None), None), ('E_fergusonii', (), (), None), None), ('S_bongori', (), (), None), None), ('root', (10542, (), (), None), ('g0', ('g1', (710, (), (), None), (14240, (), (), None), None), (6319, (), (), None), None), None), {10542: 'S_bongori', 710: 'E_coli_ATCC11775', 14240: 'E_coli_K12', 6319: 'E_fergusonii'}, {10542: 3140, 710: 3141, 14240: 3141, 6319: 3140}, {10542: 3140, 710: 3141, 14240: 3141, 6319: 3140, 'g1': 3141, 'g0': '*', 'root': '*'}, 0.3, 0.4, 0.4, 0.1, 0.2)

initFamNum,speciesTree,geneTree,tipMapD,gtLocusMapD,locusMapForRootingD,D,T,L,O,R = argT

args = (speciesTree, geneTree, tipMapD, gtLocusMapD, locusMapForRootingD, D, T, L, O, R)

speciesTreeD = trees.parseTreeForDP(speciesTree, parasite=False)
for genetree in trees.get_all_rerootings(geneTree, locusMapForRootingD):
    geneTreeD = trees.parseTreeForDP(geneTree, parasite=True)
    print(new_DTLOR_DP.DP(speciesTreeD, geneTreeD, tipMapD, gtLocusMapD, D, T, L, O, R))

