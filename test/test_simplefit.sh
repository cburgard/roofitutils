fit.py --input test/test_workspace.root --workspace Test --data obsData --hesse --correlation-matrix --print -o test.json
plotpois.py -i color=red test.json -o pois.tex --values
plotcorrelation.py -i test.json -o corrmat.tex
