python scripts/fit.py --input test/workspace.root --workspace Test --data obsData --output ../fitresults/20200123_test_obsData.txt --scan "mu" "50" "-0.3" "3" --no-findSigma
#python scripts/fit.py --input test/workspace.root --workspace Test --data asimovData --output ../fitresults/20200123_test_asimovData.txt
python scripts/plotscan.py --input "color=black" ../fitresults/20200123_test_obsData.txt --labels "\\mu" --poi "mu" --output ../plots/20200123_test_1d_obsData.tex 
