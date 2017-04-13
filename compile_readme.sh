pandoc README.md -t json | python2.7 pandocfilter-include.py | pandoc -f json -o README.md

