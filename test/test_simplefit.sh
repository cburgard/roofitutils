scripts/fit.py --input test/test_workspace.root --workspace Test --data obsData --hesse --correlation-matrix --print -o test.json
scripts/plotpois.py -i color=red test.json -o pois.tex --values
scripts/plotcorrelation.py -i ../test.json -o corrmat.tex
