scripts/fit.py --input test/test_workspace.root --workspace Test --data asimovData --hesse --scan mu 10 0 2 --print --output scan.json
scripts/plotscan.py --input color=red scan.json --output scan.tex
pdflatex scan.tex
