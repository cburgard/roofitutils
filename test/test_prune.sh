fit.py -i test/test_workspace.root -o prune.txt -d asimovData --hesse --write-result
order_NPs.py --pois mu --fitResult prune.root -o order.txt
prune_NPs.py --input test/test_workspace.root --fitResult prune.root --pois mu --snapshot-name prune_1pct --order order.txt --output pruned.root --percentages 1 --write-list pruned.txt
