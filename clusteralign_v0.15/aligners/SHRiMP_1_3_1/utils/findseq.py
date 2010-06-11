#	$Id: findseq.py 245 2008-06-06 18:24:28Z rumble $

import sys
from utils import *

if len(sys.argv) != 3:
	print >> sys.stderr, "usage: %s [contig_file] [sequence]" % (sys.argv[0])
	sys.exit(1)

seen_fasta_line = False
lineset = []
fd = open_gz_or_ascii(sys.argv[1])
for line in fd:
	if line.startswith(">"):
		if seen_fasta_line:
			print >> sys.stderr, "ERROR: this only supports one contig per file"
			sys.exit(1)
		seen_fasta_line = True
	else:
		''' Simple concatenation is wayyyyyy too slow (many orders of magnitude slower)'''
		lineset.append(line.strip())
fd.close()

contig = "".join(lineset).upper()
lineset = None

token = sys.argv[2].upper()

positions = []
start = 0
while True:
	pos = contig[start:].find(token)
	if pos == -1:
		break
	else:
		positions.append(start + pos + 1)
	start = start + pos + 1
	if start >= len(contig):
		break

if len(positions) == 0:
	print "Sequence Not Found."
else:
	print "Sequence Found at Positions:"
	for pos in positions:
		print "    " + str(pos)
