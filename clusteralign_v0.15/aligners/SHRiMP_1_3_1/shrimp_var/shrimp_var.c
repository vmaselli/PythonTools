#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>

// file handling
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>



#include "shrimp_var.h"

#define LINE_MAX 256

static int Rflag = false;
static int nlines = 0;
static FILE *outfile;

static int inputtype = NONE;

int main(int argc, char **argv) {

	char *progname = argv[0];
	outfile = stdout;
	
	
	char ch;
	
	while ((ch = getopt(argc, argv, "Ro:rpv")) != -1) {
		switch (ch) {
			case 'R':
				Rflag = true;
				break;
			case 'r': // rmapper input
				inputtype = RMAPPER_CUR;
				break;
			case 'p':
				inputtype = PROBCALC_CUR;
				break;
			case 'v':
				inputtype = RMAPPER_V09;
				break;
			case 'o':
				outfile = fopen(optarg, "w");
				break;
			default:
				usage(progname);
		}
	}
	
	if (inputtype == NONE) {
		usage(progname);
	}
	
	argc -= optind;
	argv += optind;
	
	fprintf(stderr, "#assuming format:\n"
			">readname contigname strand contigstart contigend readstart readend "
		    "readlength score editstring %snormodds pgenome pchance\n",
		    (Rflag) ? "readsequence " : "");
	
	file_iterator_n(argv, argc);
	
	fclose(outfile);
	
	return 0;
}

/*
 * look at all the passed paths, call file_iterator(path)
 */
int
file_iterator_n(char **paths, int npaths)
{
	int files;
	int i;

	for (i = files = 0; i < npaths; i++)
		files += file_iterator(paths[i]);

	return (files);
}



/*
 * Given a path, if it's a regular file (or symlink to one), call fh on it.
 * If it's a directory, call fh on all regular files within it (or symlinks to
 * regular files).
 *
 * Returns the number of files fh was called on.
 */
int file_iterator(char *path) {
	char fpath[2048];
	struct stat sb;
	DIR *dp;
	struct dirent *de;
	int files;

	/* is a regular file... */
	if (stat(path, &sb) != 0) {
		fprintf(stderr, "error: failed to stat [%s]: %s\n", path,
		    strerror(4));
		exit(1);
	}
	
	if (S_ISREG(sb.st_mode)) {
		variant_transform(path);
		return (1);
	}

	/* is (hopefully) a directory... */
	dp = opendir(path);
	if (dp == NULL) {
		fprintf(stderr, "error: failed to open directory [%s]: %s\n", path, 
				strerror(3));
		exit(1);
	}

	files = 0;
	while (1) {
		de = readdir(dp);
		
		if (de == NULL)
			break;

		fprintf(stderr, "processing file: %s\n", de->d_name);
		
		strcpy(fpath, path);
		
		if (fpath[strlen(path) - 1] != '/')
			strcat(fpath, "/");
		
		strcat(fpath, de->d_name);
		
		if (stat(fpath, &sb) != 0) {
			fprintf(stderr, "error: failed to stat [%s]: %s\n", path,
			    strerror(4));
			exit(1);
		}

		if (!S_ISREG(sb.st_mode) && !S_ISLNK(sb.st_mode))
			continue;

		/* ensure it's a regular file or link to one */
		if (S_ISREG(sb.st_mode)) {
			variant_transform(fpath);
			files++;
		} else {
			fprintf(stderr, "warning: [%s] is neither a regular "
			    "file, nor a link to one. mode:%i, ; skipping...\n", fpath, 
			    sb.st_mode);
			continue;
		}
		
	}

	closedir(dp);

	return (files);
}


/* 
 * read in the file probcalc output file in path, do output with a
 * more detailed variant info
 */
int variant_transform(char *fpath) {
	
	FILE *fp;
	
	if (strcmp("-", fpath) == 0) {
		fp = stdin;
	} else {
		fp = fopen(fpath, "r");
	}
	if (fp == NULL) {
		fprintf(stderr, "error: could not open file [%s]: %s\n",
		    fpath, strerror(3));
		exit(1);
	}
	
	int nlines_local = 0;
	
	// reading in variables
	char readname[256], contigname[256], strand[256], editstring[100];
	char readsequence[256];
	long contigstart = -1; 
	long contigend = -1; 
	long readstart  = -1; 
	long readend = -1;
	
	int readlength = -1;
	
	float score = -1; 
	float normodds = -1; 
	float pgenome = -1; 
	float pchance  = -1;
	
	char line[LINE_MAX];
	int linelen = 0;
	int i;
	
	char ech;
	char field[256];
	int fieldnr = 1;
	int fieldindex = 0;
	
	
	while (fgets(line, LINE_MAX, fp) != NULL) {
		
		if (line[0] == '#')
			continue;
		
		linelen = strlen(line);
		for (i = 0; i < linelen; i++) {
		
			ech = line[i];
				
			if (ech == '\t') {
				field[fieldindex] = '\0';
				fieldindex = 0;
				
				switch (fieldnr) {
					case 1:
						strcpy(readname, field);
						break;
					case 2:
						strcpy(contigname, field);
						break;
					case 3:
						strcpy(strand, field);
						break;
					case 4:
						contigstart = atoll(field);
						break;
					case 5:
						contigend = atoll(field);
						break;
					case 6:
						readstart = atoll(field);
						break;
					case 7:
						readend = atoll(field);
						break;
						
					case 8:
						readlength = atoi(field);
						break;
					case 9:
						score = atof(field);
						break;
					case 10:
						strcpy(editstring, field);
						break;
						
					case 11:
						if (Rflag)
							strcpy(readsequence, field);
						else
							normodds = atof(field);
						break;
						
					case 12:
						if (Rflag)
							normodds = atof(field);
						else
							pgenome = atof(field);
						break;
						
					case 13:
						if (Rflag)
							pgenome = atof(field);
						else
							pchance = atof(field);
						break;
						
					case 14:
						if (Rflag)
							pchance = atof(field);
						else {
							fprintf(stderr, "no Flag, and too many fields");
							exit(1);
						}
						break;
						
					default:
						fprintf(stderr, "error: failed to parse file line [%s] fieldnr:%i\n", line, fieldnr);
						break;
						
				}	
				fieldnr++;
				fieldindex = 0;
				
			} else {
				field[fieldindex++] = ech;
				
			}
		}
		
		fieldnr = 1;
		fieldindex = 0;
		
		nlines_local++;
		fprintf(outfile, "%s\t%s\t%li", readname, editstring, contigstart);
		editstr_to_stats(editstring, contigstart, (strcmp(strand, "+") == 0));
		fprintf(outfile, "\n");
			
	}
	
	nlines += nlines_local;
	return nlines_local;
}



/*
 * usage information
 */
static void usage(char *progname) {
	char *slash;

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage: %s (-v|-p|-r) [-R] -o outfile results_dir1|results_file1 "
	    "results_dir2|results_file2 ...\n", progname);
	exit(1);
}



/* 
 * shrimp edit string to stats 
 */  
void editstr_to_stats(char * str, long readloc, int is_forward) {

	// control variables
	int inins = 0; 
	int indel = 0;
	int innum = 0;

	// counts
	int nr_snps = 0;
	int nr_ins = 0;
	int nr_dels = 0;
	
	// temporary sections
	int inssize = -1; 
	int delsize = -1;
	char num[5];
	char rev_num[5];
	char ins[100];
	char rev_ins[100];
	
	// current character
	char ech;
	
	// current outstring
	char outstring[256];
	outstring[0] = '\0';
	
	int isnuc = -1;

	int i, j;
	int slen = strlen(str);
	for (i = 0; i < slen; i++) {
		
		if (is_forward) {
			ech = str[i];
		} else {
			ech = str[slen - i - 1];
		}
			
		
		// make sure we have a valid character
		if(!assert_editstring_char(ech)) {
			fprintf(stderr, "Unrecognized character: %c\n", ech);
			exit(1);
		}
		
		// if its a number
		if (ech == '1' || ech == '2' || ech == '3' || ech == '4' ||
				ech == '5' || ech == '6' || ech == '7' || ech == '8' || 
				ech == '9' || ech == '0') {
			num[innum] = ech;
			innum++;
			num[innum] = '\0';
			
		} else if (innum > 0) { // not in a number, but was
			if (is_forward) {
				readloc += atoi(num);
			} else {
				for (j = 0; j < (int) strlen(num); j++) {
					rev_num[j] = num[strlen(num) - j - 1];
				}
				rev_num[strlen(num)] = '\0';
				readloc += atoi(rev_num);
			} 
			innum = 0;
			num[0] = '\0';
		}
		
		// check if it is a nucleotide
		isnuc = (ech == 'A' || ech == 'C' || ech == 'T' || ech == 'G' || ech == 'N');
		
		// SNP
		if (!inins & isnuc) {
			nr_snps++;
			if (is_forward) {
				sprintf(&outstring[strlen(outstring)], "s-%c-%li\t", ech, readloc);
			} else { 
				sprintf(&outstring[strlen(outstring)], "s-%c-%li\t", complement(ech), readloc);
			}
			readloc++;
			continue;
		}
		
		// DELETIONS
		if (ech == '-' && indel) {
			delsize++;
			continue;
		} else if (ech == '-') { // just starting a deletion
			indel = 1;
			delsize = 1;
			
		} else if (indel) { // was in deletin, but not a deletion anymore
			indel = 0;
			//write deletion down
			sprintf(&outstring[strlen(outstring)], "d-%i-%li\t", delsize, readloc);
			nr_dels++;
			readloc += delsize;
			delsize = 0;
		} // ow not in deletion and not a delection char '-'  

		// INSERTIONS
		// opening bracket
		if ((is_forward && ech == '(') || (!is_forward && ech == ')')) {
			assert(!indel); // should not be in a deletion
			inins = 1;
			inssize = 0;
			
		} else if (isnuc && inins) {
			ins[inssize] = ech;
			inssize++;
			ins[inssize] = '\0';
			continue;
			
		// closing bracket
		} else if ((is_forward && ech == ')') || (!is_forward && ech == '(')) {
			nr_ins++;
			if (is_forward) {
				sprintf(&outstring[strlen(outstring)], "i-%s-%li\t", ins, readloc - 1);
			} else {
				for (j = 0; j < (int) strlen(ins); j++) {
					rev_ins[j] = complement(ins[j]); // no need for ins[strlen(ins) - j - 1], since already reading backwards
				}
				rev_ins[strlen(ins)] = '\0';
				sprintf(&outstring[strlen(outstring)], "i-%s-%li\t", rev_ins, readloc - 1);
			} 
			inins = 0;
			inssize = 0;
			continue;
		}
	}	
	
	fprintf(outfile, "\t%i %i %i\t", nr_snps, nr_ins, nr_dels);
	fprintf(outfile, outstring);
}

/*
 * Assert that it is the type of letter we expect
 */
int assert_editstring_char(char echar) {
	return
		echar == 'A' || echar == 'C' || echar == 'G' || echar == 'T' ||
	        echar == 'N' || 
		echar == '1' || echar == '2' || echar == '3' || echar == '4' || 
		echar == '5' || echar == '6' || echar == '7' || echar == '8' || 
		echar == '9' || echar == '(' || echar == ')' || echar == '-' ||
		echar == 'x' || echar == '0';
		
}

/* 
 * return the DNA complement residue
 */
int complement(char ech) {
	if (ech == 'A') return 'T';
	if (ech == 'T') return 'A';
	if (ech == 'C') return 'G';
	if (ech == 'G') return 'C';
	if (ech == 'N') return 'N';
	assert(0);
	return -1;
}


