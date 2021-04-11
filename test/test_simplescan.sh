scripts/fit.py --input test/test_workspace.root --workspace Test --data asimovData --hesse --scan mu 10 0 2 --print --output scan1d.json
scripts/plotscan.py --input color=red scan.json --output scan1d.tex
# pdflatex scan1d.tex
