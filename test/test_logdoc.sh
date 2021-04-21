scripts/fit.py --input test/test_workspace.root --workspace Test --data obsData --hesse --logsave log.txt
scripts/logdoc.py log.txt --fail-on-print --print unknown
