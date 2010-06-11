#	$Id: colourise.py 245 2008-06-06 18:24:28Z rumble $

import sys
from utils import *

MAXCOLS		= 71
BASE_N_COLOUR	= 4
BASE_T		= 3
ALT_BASE_NUM	= 4

colourmat = (	(0, 1, 2, 3, BASE_N_COLOUR),	\
		(1, 0, 3, 2, BASE_N_COLOUR),	\
		(2, 3, 0, 1, BASE_N_COLOUR),	\
		(3, 2, 1, 0, BASE_N_COLOUR),	\
		(BASE_N_COLOUR, BASE_N_COLOUR,	\
		 BASE_N_COLOUR, BASE_N_COLOUR,	\
		 BASE_N_COLOUR))

lettermap = { "A" : 0, "a" : 0, "C" : 1, "c" : 1, "G" : 2, "g" : 2, \
	      "T" : 3, "t" : 3 }

if len(sys.argv) != 2:
	print >> sys.stderr, "usage: %s input_file" % (sys.argv[0])
	sys.exit(1)

col = 0
last = BASE_T
str = ""

fd = open_gz_or_ascii(sys.argv[1])
for line in fd:
	if line.startswith(">"):
		if len(str) > 0:
			print str
		col = 1
		print "%s" % (line.strip())
		str = "T"
		last = BASE_T
		continue

	for char in line.strip():
		if col == MAXCOLS:
			col = 0
			print str
			str = ""

		letter = -1
		if char in lettermap:
			letter = lettermap[char]
		else:
			letter = ALT_BASE_NUM

		col = col + 1

		str = str + "%d" % (colourmat[last][letter])
		last = letter

if len(str) > 0 and str != "T":
	print str
