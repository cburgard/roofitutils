fit.py --input test/test_workspace.root --workspace Test --data asimovData --hesse --scan mu 20 0 2 --scan norm_bkg 20 0.9 1.1 --print --output scan2d.json
plotscan.py --input color=red scan2d.json --output scan2d.tex
# pdflatex scan2d.tex
