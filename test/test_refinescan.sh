###This is a test code to help using the --refine-scan flag in fit.py
###It can be used to scan a workspace using as a starting point the contour of a limit 
###as a first step creates a nice round contour:
fit.py --input test/test_workspace.root --workspace Test --data asimovData --hesse --scan mu 20 -2 4 --scan norm_bkg 20 0.9 1.1 --print --output scan2d.json

###creates the pointsto scan
fit.py --input test/test_workspace.root --workspace Test --data asimovData --print --refine-scan scan2d.json --output scan2d_ref.json --write-submit out_test_refine_scan --job-size 1000 --refine-scan-threshold 95 --refine-scan-spread 0.001

###plot some points and the contour
plotscan.py --points color=black  out_test_refine_scan/coords_0.txt -i color=red scan2d.json -o scan2d_points.tex

###or launch the scans with:
#source out_test_refine_scan/jobs.txt

# pdflatex scan2d.tex
