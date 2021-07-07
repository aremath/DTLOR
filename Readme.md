This directory is for running DTLOR detached from the rest of xenoGI.

Run like so:

```bash
$ python3 runDTLOR_DP.py speciesTree.tre initFam001601.tre
$ python3 runDTLOR_DP.py speciesTree.tre initFam001601.tre -d 2 -t 3 -l 1 -o 2 -r 1
```

Can also specify the locus mapping or the tip mapping with `--mapping` and `--locus` respectively.

Reconciliation is done using an unrooted gene tree, rooting it in every possible way and reconciling that with the species tree. (we eliminate some possible rootings if there are subtrees where every tip has the same locus.)
