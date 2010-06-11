#	$Id: cs2ls.py 291 2008-11-18 20:22:32Z rumble $

import sys

def print_ls(line):
	ls = ""

	lookup = {}
	lookup['A'] = {'0' : 'A', '1' : 'C', '2' : 'G', '3' : 'T' }
	lookup['C'] = {'0' : 'C', '1' : 'A', '2' : 'T', '3' : 'G' }
	lookup['G'] = {'0' : 'G', '1' : 'T', '2' : 'A', '3' : 'C' }
	lookup['T'] = {'0' : 'T', '1' : 'G', '2' : 'C', '3' : 'A' }

	if line[0] != 'A' and line[0] != 'C' and \
	   line[0] != 'G' and line[0] != 'T':
		print "Error: invalid cs read [%s]" % (line)
		return;

	last = line[0]
	first = True
	for c in line:
		if first:
			ls = ls + c
			first = False
			continue

		if c != '0' and c != '1' and \
		   c != '2' and c != '3':
			print "Error: invalid cs read [%s]" % (line)
			return
		
		last = lookup[last][c]
		ls = ls + last

	print "Colourspace: %s" % (line)
	print "Letterspace: %s" % (ls)

for line in sys.stdin:
	line = line.strip()
	if len(line) > 0:
		print_ls(line.upper())
