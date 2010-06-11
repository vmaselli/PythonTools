#	$Id: splittigs.py 245 2008-06-06 18:24:28Z rumble $

import sys
from utils import *

if len(sys.argv) != 2:
	print >> sys.stderr, "usage: %s [fasta_file]" % (sys.argv[0])
	sys.exit(1)

outfd = None
nfiles = 0
skipping = False

fd = open_gz_or_ascii(sys.argv[1])
for line in fd:
	if line.startswith(">"):
		fname = line[1:].strip() + ".fa"
		if outfd != None:
			outfd.close()
		outfd = open(fname, "w")
		print >> sys.stderr, "splitting into file [%s]" % (fname)
		nfiles = nfiles + 1

	if outfd == None:
		if not skipping:
			print >> sys.stderr, "warning: no contig label yet; skipping line"
			skipping = True
		continue
	skipping = False

	outfd.write(line)

if outfd != None:
	outfd.close()

if nfiles != 0:
	print >> sys.stderr, "------------------------------------------"
	print >> sys.stderr, "created %d individual contig file(s)" % (nfiles)
	print >> sys.stderr, "------------------------------------------"
