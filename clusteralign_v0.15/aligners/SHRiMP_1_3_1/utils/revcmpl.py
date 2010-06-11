#	$Id: revcmpl.py 245 2008-06-06 18:24:28Z rumble $

import sys
from utils import *

WIDTH = 50

if len(sys.argv) != 2:
	print >> sys.stderr, "usage: %s [file_name]" % (sys.argv[0])
	sys.exit(1)

def revcmpl(sequence):
	lookup = []

	i=0
	while i < 256:
		if i == ord('A'):
			lookup.append('T')
		elif i == ord('a'):
			lookup.append('t')
		elif i == ord('C'):
			lookup.append('G')
		elif i == ord('c'):
			lookup.append('g')
		elif i == ord('G'):
			lookup.append('C')
		elif i == ord('g'):
			lookup.append('c')
		elif i == ord('T'):
			lookup.append('A')
		elif i == ord('t'):
			lookup.append('a')
		elif i == ord('N'):
			lookup.append('N')
		elif i == ord('n'):
			lookup.append('n')
		else:
			lookup.append('!')
		i = i + 1

	rc_list = ['!'] * len(sequence)
	i = len(sequence) - 1
	j = 0
	while i >= 0:
		rc_list[j] = lookup[ord(sequence[i])]
		i = i - 1
		j = j + 1

	i = 0
	while i < len(sequence):
		if rc_list[i] == '!':
			print >> sys.stderr, "error: saw unhandled base [%c]" % (sequence[i])
			sys.exit(1)
		i = i + 1

	revcmpl_sequence = "".join(rc_list)
	if len(revcmpl_sequence) != len(sequence):
		print >> sys.stderr, "REVCMPL ERROR!!!"
		sys.exit(1)	

	return (revcmpl_sequence)

lines = []
seen_contig = False

fd = open_gz_or_ascii(sys.argv[1])
for line in fd:
	if line.startswith(">"):
		if seen_contig:
			print >> sys.stderr, "error: only one fasta line per file is supported"
			sys.exit(1)
		seen_contig = True
		contig = line.strip() + "_revcmpl"
		print contig
		continue
	elif line.startswith("#"):
		print line
		continue
	elif line.strip() == "":
		continue

	lines.append(line.strip())

fd.close()
lines.reverse()

cols = 0
str = ""

for line in lines:
	for char in revcmpl(line):
		str = str + char
		cols = cols + 1
		if cols == WIDTH:
			if str != "":
				print str
			str = ""
			cols = 0

if str != "":
	print str
