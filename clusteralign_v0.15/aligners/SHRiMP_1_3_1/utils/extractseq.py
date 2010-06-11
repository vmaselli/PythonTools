#	$Id: extractseq.py 245 2008-06-06 18:24:28Z rumble $

import sys
from utils import *

if len(sys.argv) != 4:
	print >> sys.stderr, "usage: %s [contig_file] [start_index] [end_index]" % (sys.argv[0])
	sys.exit(1)

start = int(sys.argv[2])
if start < 1:
	print >> sys.stderr, "error: start index must be >= 1"
	sys.exit(1)

end = int(sys.argv[3])
if end < start:
	print >> sys.stderr, "error: end index less than start index"
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

contig = "".join(lineset)
lineset = None

''' Input is 1 to n, indexing is 0 to n-1 '''
start = start - 1

if start >= len(contig):
	print >> sys.stderr, "error: start index is greater than contig length"
	sys.exit(1)

if end >= len(contig):
	print >> sys.stderr, "error: end index is greater than contig length"
	sys.exit(1)

print >> sys.stderr, "Outputting %d bases (%d to %d, inclusive) of %d contig bases" % (end - start, start + 1, end, len(contig))
print contig[start : end]
