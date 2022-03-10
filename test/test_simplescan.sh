fit.py --input test/test_workspace.root --workspace Test --data asimovData --hesse --scan mu 10 0 2 --print --output scan1d.json
plotscan.py --input color=red scan1d.json --output scan1d.tex --label "\mu"
# pdflatex scan1d.tex
