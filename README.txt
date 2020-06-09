
This directory is for running DTLOR detached from the rest of xenoGI.

Run like so:

python3 runDTLOR_DP.py speciesTree.tre initFam001601.tre # small tree, 1.6 sec on purves
python3 runDTLOR_DP.py speciesTree.tre initFam000220.tre # larger tree, 11.6 sec on purves
python3 runDTLOR_DP.py speciesTree.tre initFam000060.tre # still larger, 51 min on purves
python3 runDTLOR_DP.py speciesTree.tre initFam000001.tre # big tree, ~2 days on purves

The files trees.py and familiesDTLORstuff.py contain some necessary functions from xenoGI. You shouldn't need to modify these.

FYI, reconcile within familiesDTLORstuff.py is taking an unrooted gene tree, rooting it in every possible way and reconciling that with the species tree. (we do eliminate some possible rootings if there are subtrees where every tip has the same locus.)

DTLOR_DP.py has the dtlor code.
