#	$Id: extractunmapped.py 335 2009-05-15 18:51:03Z rumble $

import sys
from utils import *

if len(sys.argv) != 2:
	print >> sys.stderr, "usage: %s [rmapper_output_file]" % (sys.argv[0])
	sys.exit(1)

seen_unmapped_comment = False
fd = open_gz_or_ascii(sys.argv[1])
for line in fd:
	if line.startswith("#UNMAPPED READS:"):
		if seen_unmapped_comment:
			print >> sys.stderr, "ERROR: multiple unmapped comments in one file?!?"
			sys.exit(1)
		seen_unmapped_comment = True
	if seen_unmapped_comment and not line.startswith("#"):
		print line.strip()
