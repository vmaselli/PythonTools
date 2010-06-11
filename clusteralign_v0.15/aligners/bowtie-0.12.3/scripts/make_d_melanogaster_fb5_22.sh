#!/bin/sh

#
# Downloads sequence for a D. melanogaster from flybase.  Currently set
# to download 5.22, but F, REL, and IDX_NAME can be edited to reflect a
# different version number.  (But note that you will usually also have
# to change the date in REL.)
#

GENOMES_MIRROR=ftp://ftp.flybase.net/genomes/Drosophila_melanogaster
F=dmel-all-chromosome-r5.22.fasta
REL=dmel_r5.22_FB2009_09
IDX_NAME=d_melanogaster_fb5_22

get() {
	file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o `basename $1` $1
		return $?
	else
		wget -O `basename $1` $1
		return $?
	fi
}

BOWTIE_BUILD_EXE=./bowtie-build
if [ ! -x "$BOWTIE_BUILD_EXE" ] ; then
	if ! which bowtie-build ; then
		echo "Could not find bowtie-build in current directory or in PATH"
		exit 1
	else
		BOWTIE_BUILD_EXE=`which bowtie-build`
	fi
fi

if [ ! -f $F ] ; then
	FGZ=$F.gz
	get ${GENOMES_MIRROR}/$REL/fasta/$FGZ || (echo "Error getting $FGZ" && exit 1)
	gunzip $FGZ || (echo "Error unzipping $FGZ" && exit 1)
fi

CMD="${BOWTIE_BUILD_EXE} $* $F $IDX_NAME"
echo "Running $CMD"
if $CMD ; then
	echo "$IDX_NAME index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi
