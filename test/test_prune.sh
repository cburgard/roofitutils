python scripts/fit.py -i test/test_workspace.root -o prune.txt -d asimovData --hesse --write-result
python scripts/order_NPs.py --pois mu --fitResult prune_fitresult.root --output order.txt
python scripts/prune_NPs.py --input test/test_workspace.root --fitResult prune_fitresult.root --pois mu --snapshot-name prune_1pct --order order.txt --output pruned.root --percentages 1
