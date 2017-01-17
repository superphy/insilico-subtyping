pandoc README.md -t json | python pandocfilter-include.py | pandoc -f json -o README.md

