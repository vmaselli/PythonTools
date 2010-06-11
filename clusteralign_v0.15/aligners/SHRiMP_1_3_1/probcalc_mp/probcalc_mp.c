// $Id: probcalc_mp.c 383 2009-09-30 22:47:56Z matei $

/****************************************************************************
 * The probcalc_mp project takes in a binary file filled with mapping_t 
 * structures as defined in dbtypes.h. This is a file with read mappings as
 * outputted by SHRiMP's probcalc tool. Assuming the binary file is sorted by 
 * read id, probcalc_mp will find all the mate pair mappings, assign them 
 * probability-based confidence values, and output them. How and which 
 * computations are done, and which matepair are outputted depends on the 
 * given parameters. Due to the fact that this software is aimed at being used 
 * in current NGS research, many parameters are provided. See accompanying 
 * readme.
 ***************************************************************************/

#include <alloca.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <inttypes.h>

// file handling
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <dirent.h>
#include <sys/mman.h>

#include "dbtypes.h"
#include "probcalc_mp.h"

// global statistics
static double gl_mean  = 0;				// mean of distances
static double gl_stdev = 0;				// standard deviaiton
static uint64_t gl_hist[HIST_BINS]; 	// histogram of distances
static double gl_hist_cumsum[HIST_BINS];// cumulative histogram of distances from mean
static int do_gl_hist = 1;				// boolean whether to do histogram
static uint64_t gl_good_mps  = 0;		// number of good mappings (good: d < M, ++, --) 
static int gl_done_mean = 0;			// bool: done the mean calculation
static uint64_t gl_uniq_reads = 0;		// bool: done the mean calculation

// settings
static int do_unique = 0;				// calc mean only on unique mappings 
static int print_max = PRINT_MAX;		// the number of matepair scores to print
static char * fwd_suffix = (char *)"";			// forward suffix
static int fwdsuflen = 0;				// forward suffix length
static char * rev_suffix = (char *)"";			// reverse (backward) suffix
static int revsuflen = 0;				// reverse suffix length
static int Rflag = false;				// Rflag included in probcalc output
static uint64_t distcutoff = 0;			// hard distance cutoff (M).
static uint64_t hist_distcutoff = 0;	// hard distance cutoff used in histograms (M).
static int max_reads = MAX_READS;		// maximum number of a read's mappings
static uint64_t genome_length = 0;		// genome length
static int discordant = 0;				// output only discordant or only concordant mate pairs
static int sort_field = SORT_PGENOME;	// the field to sort the results by
static uint64_t gl_mean_nr = MEAN_MPS;	// number of USEFUL mp to include in mean comp
static uint64_t gl_printed_mp = 0; 		// the number of printed mps
static double pgenome_cutoff = PGENOME_CUTOFF; // pchance cutoff
static double pchance_cutoff = PCHANCE_CUTOFF; // pgenome cutoff 
static int quickmode = 0;				// TODO boolean quickmode option
static double nr_stdev = 2;				// number of standard deviations for second M assignment
static int allow_diff_chr = 1;			// allow the same chromosome or not
static int input_file_type = ASCII;	// input file type, ascii or binary

// debug variables
static uint64_t gl_debug_call = 0;		// the number of times mp_analysis is called
static FILE *outdebugfp;				// debugging output file
static int debuglimit = DEBUGLIMIT;		// if DEBUG is on, the number of debugs before quitting 
static int debugmode = 0;				// debug mode...
static int debugseek = 0;				// seek into the read file
static int gl_debug_lines = 0;			// the number of binary "lines" read

// time variables
clock_t start_clock[TOTAL_TIME + 1];
clock_t end_clock[TOTAL_TIME + 1];
double elapsed_clock[TOTAL_TIME + 1];

//other
static int calledMP = 0;


static uint64_t
abs64(int64_t x)
{
	if (x < 0)
		return (-x);
	else
		return (x);
}

int main(int argc, char **argv) {

	start_clock[TOTAL_TIME] = clock();

	char *progname = argv[0];
	int ch;
	
	char *mappingfilename = (char *)"";
	char * inputfiletype = (char *)"ascii";
	
	// check for given input
	int givenr = 1;
	int givenm = 0;
	int givenf = 0;
	int givenb = 0;
	int giveng = 0;
	int givenM = 0;
	
	while ((ch = getopt(argc, argv, "m:x:R:f:b:M:g:duL:T:D:C:G:qs:cei:")) != -1) {
		switch (ch) {
		
			case 'R': // Rflag was included in probcalc runs
				Rflag = true;
				break;
			case 'x': // max number of reads to expect, per direction
				max_reads = atoi(optarg);
			case 'd': // output only discordant mate pairs
				discordant = 1;
				break;		
			case 'u': // compute mean with only unique matepair mappings
				do_unique = 1;
				break;		
			case 'L': // number of USEFUL reads to include in mean comp.
				gl_mean_nr = atoll(optarg); // gl_mean_nr = 0 implies look at all
				break;		
			case 'T': // maximum matepairs to print
				print_max = atoi(optarg);
				break;
			case 'G': // PGENOME CUTOFF
				pgenome_cutoff = atof(optarg);
				break;
			case 'C': // PCHANCE CUTOFF
				pchance_cutoff = atof(optarg);
				break;
			case 'c': //don't allow same chromos
				allow_diff_chr = 0;
				break;
			case 'm': // binary mapping file
				mappingfilename = optarg;
				givenm = 1;
				break;
			case 'f': // forward suffix
				fwd_suffix = strdup(optarg);
				fwdsuflen = strlen(fwd_suffix);
				givenf = 1;
				break;
			case 'b': // reverse (backward) suffix
				rev_suffix = strdup(optarg);
				revsuflen = strlen(rev_suffix);
				givenb = 1;
				break;
			case 'g': // genome length
				genome_length = atoll(optarg);
				giveng = 1;
				break;
			case 'D': // debug mode
				debugmode = 1;
				debuglimit = atoll(optarg);
				break;
			case 'q': // quick mode
				quickmode = 1;
				break;
			case 's': // the number of standard deviations for new M comp.
				nr_stdev = atof(optarg);
				break;
			case 'i': // input file type. Default: ASCII
				inputfiletype = optarg;
				if (strcmp(inputfiletype, "ascii") == 0) 
					input_file_type = ASCII;
				else 
					input_file_type = BINARY;
				
			case 'M': // hard distance cut off. 
				distcutoff = atoll(optarg);
				hist_distcutoff = distcutoff;
				givenM = 1;
				break;			
				
			default:
				usage(progname);
		}
	}

	// checking of needed parameters
	if (!givenr) { fprintf(stderr, "Read filename not given\n"); usage(progname); exit(1); }
	if (!givenm) { fprintf(stderr, "Mapping filename not given\n"); usage(progname); exit(1); }
	if (!givenf) { fprintf(stderr, "Fwd suffix not given\n"); usage(progname); exit(1); }
	if (!givenb) { fprintf(stderr, "Rev suffix not given\n"); usage(progname); exit(1); }
	if (!giveng) { fprintf(stderr, "Genome Length not given\n"); usage(progname); exit(1); }
	if (!givenM) { fprintf(stderr, "Distance cutoff not given\n"); usage(progname); exit(1); }
	
	// print settings
	print_dashed_line();
	fprintf(stderr, "Some settings:\n");
	fprintf(stderr, "pgenome cutoff:  %f\n", pgenome_cutoff);
	fprintf(stderr, "pchance cutoff:  %f\n", pchance_cutoff);
	
	if (input_file_type == ASCII)
		fprintf(stderr, "input file type: %s\n", "ASCII");
	else
		fprintf(stderr, "input file type: %s\n", "Binary");
	
	print_dashed_line();
	
	
	// open debug file
	outdebugfp = fopen("debug.log", "w");
	
	// iterator variables
	int i;
	
	// time init
	for (i = 0; i < TOTAL_TIME; i++) 
		elapsed_clock[i] = 0;
	
	
	
	/*******************************************************************
	 * Do The Work...                                                  *
	 *******************************************************************/
	uint64_t nr_reads;
	
	// compute mean, stdev, histogram 
	if (do_gl_hist) {
		fprintf(stderr, "\n");
		print_dashed_line();
		fprintf(stderr, "Computing mean, standard deviation and histogram ...\n");
		
		start_clock[STATS_TIME] = clock();
		for (i = 0; i < HIST_BINS; i++) 
			gl_hist[i] = 0;
		
		nr_reads = filepass (mappingfilename, MEAN_PASS);
		end_clock[STATS_TIME] = clock();
		elapsed_clock[STATS_TIME] = 
			((double)end_clock[STATS_TIME] - start_clock[STATS_TIME]) / CLOCKS_PER_SEC;
		
		fprintf(stderr, "Done general statistics in %.2f seconds\n"
						"      Mean:  %.2f. STDev: %.2f \n"
						"      Read:  %s lines (read mappings)\n"
						"      Found: %s total read pairs\n"
						"      Used:  %s good read pairs\n"
						"      Unique: %s reads\n"
						"      Mean histogram in: \"histogram.log\" \n"
						"      Percentile histogram in: \"chistogram.log\" \n",
						elapsed_clock[STATS_TIME], gl_mean, sqrt(gl_stdev/gl_good_mps),
						comma_integer(gl_debug_lines), comma_integer(nr_reads), 
						comma_integer(gl_good_mps), comma_integer(gl_uniq_reads));
		print_dashed_line();
	}

		
	//compute gl_hist_cumsum
	if (do_gl_hist) {
		compute_cumsum();
		print_hists();
	}		
	
	distcutoff = (uint64_t) ceil(gl_mean + (double)nr_stdev * sqrt(gl_stdev/ (double)gl_good_mps));
	fprintf(stderr, "new M cutoff: %" PRIu64 " = %.2f + %.2f * %.2f\n", 
			distcutoff, gl_mean, nr_stdev, sqrt(gl_stdev/gl_good_mps));
	
	
	// calculate statistics and produce output. This is the biggie.
	fprintf(stderr, "\n");
	print_dashed_line();
	nr_reads = filepass (mappingfilename, OUTPUT_PASS);
	print_dashed_line();
	
	/*******************************************************************
	 * Done The Work... print various "informative" output             *
	 *******************************************************************/
	
	
	
	end_clock[TOTAL_TIME] = clock();
	elapsed_clock[TOTAL_TIME] = ((double) 
			(end_clock[TOTAL_TIME] - start_clock[TOTAL_TIME])) / CLOCKS_PER_SEC;
	
	char rev[30] = "$Revision: 383 $";  	
	rev[16] = '\0';
	char *prev = &rev[11];
	
	print_dashed_line();
	fprintf(stderr, "Done probcalc_mp %s in %.2f seconds\n"
					"      Looked at: %s total Unique read pairs from %s read mappings\n"
					"      Printed a total of: %s mat pairs\n",
					prev, elapsed_clock[TOTAL_TIME], 
					comma_integer(nr_reads), comma_integer(gl_debug_lines), 
					comma_integer(gl_printed_mp));
	print_dashed_line();
	
	fprintf(stderr, "Time Spent:\n"
			"  total   time : %4.2f, \n"
			"  mapping calls: %4.2f, \n"
			"  mp_analysis  : %4.2f, \n"
			"  statistics   : %4.2f, \n", 
			elapsed_clock[TOTAL_TIME], elapsed_clock[SYS_TIME], 
			elapsed_clock[ANALYSIS_TIME], elapsed_clock[STATS_TIME]);
	print_dashed_line();
	
	return 0;
}



/*
 * Pass through the readfile and mappingfile once through the method
 * specified by pass_type (MEAN_PASS or OUTPUT_PASS)
 */
uint64_t filepass(char * mappingfilename, int pass_type) {
	
	// open files
	FILE *mappingfile = fopen(mappingfilename, "r");
	if (mappingfile == NULL) {
		fprintf(stderr, "Error: could not open readfile: %s\n", mappingfilename);
		exit(1);
	}
	struct stat fs;
	stat(mappingfilename, &fs);
	uint64_t nr_mappings = (uint64_t) (fs.st_size * 1.0 / sizeof(mapping_t)); 
	// TODO: could check if it actually rounds to an integer.
	
	if (pass_type != MEAN_PASS || gl_mean_nr == 0)
		fprintf(stderr, "Detected number of mappings: %s. Progress:\n",
				comma_integer(nr_mappings));
	
	int perc = 0;
	if (!debugmode) {
		fprintf(stderr, "|0        10        20        30        40        50"
				"        60        70        80        90       100\n");
	}
	
	
	// start later into the readfile. 
	if (debugseek > 0) {
		fseek(mappingfile, debugseek * sizeof(mapping_t), SEEK_SET);
	}
	
	// (re) setup the some counts
	gl_good_mps = 0;
	gl_debug_lines = 0;
	
	// prepare structure arrays
	// note: max_reads + 1 because we read dirrectly to these arrays, but may need to
	// switch
	//mapping_t *fwd_maps = (mapping_t *) malloc(sizeof(mapping_t) * (max_reads + 1));

	// Sun Studio does not allow variable-sized array declarations
	//mapping_t fwd_maps[max_reads + 1];
	mapping_t *fwd_maps = (mapping_t *)alloca(sizeof(mapping_t) * (max_reads + 1));
	assert(fwd_maps);
	//mapping_t *rev_maps = (mapping_t *) malloc(sizeof(mapping_t) * (max_reads + 1));
	//mapping_t rev_maps[max_reads + 1];
	mapping_t *rev_maps = (mapping_t *)alloca(sizeof(mapping_t) * (max_reads + 1));
	assert(rev_maps);
	mapping_t *mapping = fwd_maps;
	
	// indeces into the arrays. These need to be reset after each mp group 
	int fwd_index = 0; 
	int rev_index = 0;
	
	uint64_t debugwithzero = 0;
	uint64_t debugwithoutzero = 0;
	
	uint64_t nr_reads = 0;
	
	int is_forward = -1; // is forward read
	int do_analysis = 1; // boolean on whether to do theanalysis with the current read 
	
	char cur_name[READNAME_LEN]; 	// the read root we are currently adding to  
	cur_name[0] = '\0';
	char test_name[READNAME_LEN]; 	// the read just taken from fread. 
		
	// loop over all reads...
	///while(fread (mapping, sizeof(mapping_t), 1, mappingfile) == 1) {
	while(read_probcalc_line(mappingfile, mapping) == 1) {
		
		// process bar
		gl_debug_lines++;
		
		int condition = 0;
		if ( !debugmode && pass_type == MEAN_PASS && gl_mean_nr != 0) {
			if (gl_good_mps * 100.0 / gl_mean_nr > perc)
				condition = 1;
		} else if (! debugmode && (pass_type == OUTPUT_PASS || gl_mean_nr == 0)) {
			if (gl_debug_lines * 100.0 / nr_mappings > perc) {
				condition = 1;
			}
		} 
		
		if (condition) {
			if (perc % 10 == 0) 
				fprintf(stderr, "|");
			else
				fprintf(stderr, "-");
			perc++;
		}
		
		// check if the read is forward, 
		strcpy(test_name, mapping->readname);
		int readlen = strlen(test_name);
		is_forward = is_forward_test(test_name, readlen);
		
		if (is_forward) 
			test_name[readlen - fwdsuflen] = '\0';
		else 
			test_name[readlen - revsuflen] = '\0';
		
		// see if this is a new read: only compare reads up to suffixes.
		if (strcmp(cur_name, test_name) != 0) { 
			
			if (debugmode && pass_type == OUTPUT_PASS)
				fprintf(outdebugfp, "fwd_nr:%i, rev_nr:%i\n", fwd_index, rev_index);
			
			if (nr_reads > 0 && fwd_index > 0 && rev_index > 0 && do_analysis) {
				start_clock[ANALYSIS_TIME] = clock();
				mp_analysis(fwd_maps, rev_maps, fwd_index, rev_index, pass_type);
				end_clock[ANALYSIS_TIME] = clock();
				elapsed_clock[ANALYSIS_TIME] += ((double) 
						(end_clock[ANALYSIS_TIME] - start_clock[ANALYSIS_TIME]))
						/ CLOCKS_PER_SEC;
			} 
			
			if (fwd_index > 0) gl_uniq_reads++;
			if (rev_index > 0) gl_uniq_reads++;
			
			do_analysis = 1;
			
			// if doing the mean, and we're done the mean
			if (pass_type == MEAN_PASS && gl_done_mean)  {
				break;
			}
			
			// copy the new readname to cur_name (up to suffix)
			strcpy(cur_name, test_name);
			
			// if you are in debugmode, lower the debuglimit
			if (debugmode) {
				if (fwd_index * rev_index == 0) debugwithzero += fwd_index + rev_index;
				else debugwithoutzero += fwd_index + rev_index;
				
				if (nr_reads % 100000 == 0)
					fprintf(stderr, "DEBUGLINE: #reads:%"
					    PRIu64 ", #mappings_now:%i "
					    "| w/0:%" PRIu64 ", wo/0:%"
					    PRIu64 "\n", 
					    nr_reads, fwd_index * rev_index, 
					    debugwithzero, debugwithoutzero);
				
				if (pass_type != MEAN_PASS) debuglimit--;
				if (pass_type != MEAN_PASS && debuglimit < 0) break;
			}
			
			// reset variables as needed
			fwd_index = 0;
			fwd_maps[fwd_index] = *mapping;
			rev_index = 0;
	
			// the number of times mp_analysis is called
			gl_debug_call++;

			nr_reads++;
		} 
		
		// make sure the max_reads variable holds
		assert(fwd_index < max_reads);
		assert(rev_index < max_reads);
		
		// don't do the analysis if you are calculating the mean with unique 
		// mappings only, and you already have more than one mapping in the 
		// fwd (if its a fwd read) or rev (if its a reverse read), since its 
		// clearly not a unique mapping. The check is done here to avoid 
		// seeking (costly) if its not necessary.
		if (pass_type == MEAN_PASS && do_unique &&
				((is_forward && fwd_index >= 1) ||
						(!is_forward && rev_index >= 1))) {
			do_analysis = 0;
		}
		
		
		start_clock[SYS_TIME] = clock();
		
		if (is_forward && do_analysis) {
			fwd_index++;
			mapping = &fwd_maps[fwd_index]; // the new mapping space
		} else if (do_analysis) { 
			rev_maps[rev_index++] = *mapping;
		}
		
		end_clock[SYS_TIME] = clock();
		elapsed_clock[SYS_TIME] += ((double) 
				(end_clock[SYS_TIME] - start_clock[SYS_TIME])) / CLOCKS_PER_SEC;
	}
	
	if (!debugmode)
		fprintf(stderr, "|\n");
	
	// close the files
	fclose(mappingfile);
	
	return nr_reads;
}



/*
 * The Mate Pair analysis. This does one or two passes over the combinations
 * of mate pairs...
 */
void mp_analysis(mapping_t *fwd_maps, mapping_t *rev_maps, 
		int fwd_nr, int rev_nr, int pass_type) {

	// iteration variables
	int i,j;
		
	// the number of "good" mp mappings.
	// good: d < M and R+F+, F-R-
	int good_mps = 0;  
	uint64_t good_mps_dist = 0; // the distance of a good mp
	uint64_t dist = 0; 			// local distance
	
	// looking in every combination
	if (pass_type == MEAN_PASS || discordant) {
		for(i = 0; i < fwd_nr; i++) {
			for (j = 0; j < rev_nr; j++) {
				
				// compute the distance between these two IF the mate pair is good
				// good: (d < M) and (R+F+ or F-R-)
				// if the mate pair is not good, good_mp_dst returns 0
				dist = good_mp_dst(&fwd_maps[i], &rev_maps[j]);
				
				if (dist > 0) {
					good_mps_dist = dist;
					good_mps++;
				}
				
				if (pass_type == MEAN_PASS && good_mps > 1)
					break; 
			}
			
			if (pass_type == MEAN_PASS && good_mps > 1)
				break;
		}
	}
	
	// add to m, stdev, etc
	if (pass_type == MEAN_PASS && good_mps == 1) {
		increments_stats(good_mps_dist);
	}
	
	// if the discordant flag is NOT set, 
	// OR discordant is set and the mp is NOT conc2, 
	if (pass_type == OUTPUT_PASS && (!discordant || good_mps == 0)) {   
		
		double totnormodds = 0;
		int mp_set_index = 0;
		
		mate_pair_val *mp_set = (mate_pair_val *) 
				malloc(sizeof(mate_pair_val) * fwd_nr * rev_nr);
		assert(mp_set);
		
		// look through all the combinations
		for(i = 0; i < fwd_nr; i++) {
			for (j = 0; j < rev_nr; j++) {
				add_p_stats(&fwd_maps[i], &rev_maps[j], mp_set, 
						&totnormodds, &mp_set_index);
			}
		}
				
		// fix norm odds (normalize)
		for (i = 0; i < mp_set_index; i++) {
			mp_set[i].normodds = mp_set[i].normodds/totnormodds; 
		}
		
		//qsort
		qsort(mp_set, mp_set_index, sizeof(*mp_set), mate_pair_val_cmp);
		
		// output
		if (!calledMP) {
			printf("#FORMAT: fwd_name fwd_chr fwd_editstring fwd_strand fwd_start fwd_end fwd_pg"
					"rev_name rev_chr rev_editstring rev_strand rev_start rev_end rev_pg"
					"distance normodds pgenome pchance\n");
			calledMP = true;
		}
		
		// output
		for (i = 0; i < mp_set_index; i++) {
			
			if (i >= print_max) {
				if (mate_pair_val_cmp(&mp_set[i-1], &mp_set[i]) != 0)
					break;
			}
			
			
			printf("%lli\t", (long long int) gl_printed_mp);
			gl_printed_mp++;
			printf("%s\t%s\t%s\t%c\t%lli\t%lli\t%1.3f\t", 
					&(mp_set[i].fwd_rs->readname[1]), mp_set[i].fwd_rs->contigname,
					mp_set[i].fwd_rs->editstring, mp_set[i].fwd_rs->strand, 
					(long long int) mp_set[i].fwd_rs->contigstart, 
					(long long int) mp_set[i].fwd_rs->contigend, mp_set[i].fwd_rs->pgenome);
			printf("%s\t%s\t%s\t%c\t%lli\t%lli\t%1.3f\t", 
					&(mp_set[i].rev_rs->readname[1]), mp_set[i].rev_rs->contigname,
					mp_set[i].rev_rs->editstring, mp_set[i].rev_rs->strand, 
					(long long int) mp_set[i].rev_rs->contigstart, 
					(long long int) mp_set[i].rev_rs->contigend, mp_set[i].rev_rs->pgenome);
			printf("%lli\t%1.3f\t%1.3f\t%1.10f\n",
					(long long int) mp_set[i].dist, 
					mp_set[i].normodds, mp_set[i].pgenome, mp_set[i].pchance);
		}
		
		free(mp_set);
	}
	
}

/*
 * test if the string readname is a forward or a reverse read
 */ 
int is_forward_test(char * readname, int readlen) {
	char * fwd_inputsuff = &readname[readlen - fwdsuflen];
	char * rev_inputsuff = &readname[readlen - revsuflen];
	
	if (strcmp(fwd_suffix, fwd_inputsuff) == 0) {
		return true;
	} else if (strcmp (rev_suffix, rev_inputsuff) == 0) {
		return false;
	} else {
		fprintf(stderr, "error: read is neither forward nor reverse\n");
		fprintf(stderr, "read name: %s   fwd_suffix: %s   rev_suffix: %s\n",
				readname, fwd_suffix, rev_suffix);
		exit(1);
	}
	
}

/*
 * compare mate pairs
 */
int mate_pair_val_cmp(const void *a, const void *b) {
	
	mate_pair_val *mppva = (mate_pair_val *)a;
	mate_pair_val *mppvb = (mate_pair_val *)b;
	double cmp_a, cmp_b;

	// sort byt the correct field
	switch (sort_field) {
		case SORT_PGENOME:
			/* reversed - want big first */
		    cmp_a = mppvb->pgenome;
			cmp_b = mppva->pgenome;
			break;
		case SORT_PCHANCE:
			cmp_a = mppva->pchance;
			cmp_b = mppvb->pchance;
			break;
		case SORT_NORMODDS:
			/* reversed - want big first */
		    cmp_a = mppvb->normodds;
			cmp_b = mppva->normodds;
			break;
			
		default:
			/* shut up, gcc */
			cmp_a = cmp_b = 0;
			assert(0);
	}

	if (cmp_a == cmp_b)
		return (0);
	else if (cmp_a < cmp_b)
		return (-1);
	else
		return (1);
}

/*
 * usage of the program.
 */
static void usage(char *progname) {
	char *slash;

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage for Revision 1.166. Current $Revision: 383 $  \n"
			"%s -m binary_mapping_filename -f forward_suffix -b reverse_suffix"
			" -g genome_length -M hard_distance_limit [-L nr_mate_pairs] [-q] "
			"[-C PCHANCE_CUTOFF] [-G PGENOME_CUTOFF] [-R] "
			"[-x max_reads_to_expect] [-d] [-u] [-D] [-T max_reads_to_output] "
			"[-s nr_stdev] [-c]", progname);
	exit(1);
}

/*
 * Compute the cumsum histogram
 */
void compute_cumsum() {
	uint64_t subtract;
	gl_hist_cumsum[0] = 1;
	int mean_bin = (int) floor((gl_mean * 1.0 / hist_distcutoff) * HIST_BINS);
	
	int i;
	for (i = 1; i < HIST_BINS; i++) {
		subtract = 0;
		if (mean_bin + (i-1) < HIST_BINS)
			subtract += gl_hist[mean_bin + (i-1)];
		if (mean_bin - (i-1) > 0 && (i-1) != 0)
			subtract += gl_hist[mean_bin - (i-1)];
		
		gl_hist_cumsum[i] = gl_hist_cumsum[i-1] - (subtract*1.0 / gl_good_mps);
		gl_hist_cumsum[i] = MAX(gl_hist_cumsum[i], 0);
	}
}

/*
 * Print histograms to log files.
 */
void print_hists() {
	FILE *hfp = fopen("histogram.log", "w");
	if (hfp == NULL) {
		fprintf(stderr, "Error: could not open histfile: histogram.log\n");
		return;
	}
	
	FILE *cfp = fopen("chistogram.log", "w");
	if (cfp == NULL) {
		fprintf(stderr, "Error: could not open chistfile: chistogram.log\n");
		return;
	}
	
	int i;
	for (i = 0; i < HIST_BINS; i++) {
		fprintf(hfp, "[bin %i: %.2f-%.2f] %" PRIu64 "\n", 
				i, i*1.0*hist_distcutoff/HIST_BINS, 
				(i+1)*1.0*hist_distcutoff/HIST_BINS, gl_hist[i]);
		fprintf(cfp, "[bin %i: %.2f-%.2f] %f\n", 
				i, i*1.0*hist_distcutoff/HIST_BINS, 
				(i+1)*1.0*hist_distcutoff/HIST_BINS, gl_hist_cumsum[i]);
	}
	
	fclose(hfp);
	fclose(cfp);	
}

/*
 * if the mate pair is "good" (see definition of good above), 
 * the method returns the distance between the two, including the reads
 * themselves. 
 * 
 * otherwise, returns 0
 */
uint64_t good_mp_dst(mapping_t * fwd_map, mapping_t * rev_map) {

	// testing variables
	uint64_t dist; // distance
	uint64_t cs_fwd; // contigstart/end of fwd read
	uint64_t cs_rev; // contigstart/end of rev read
	char s_fwd; // strand of fwd read
	char s_rev; // strand of reverse read
	
	// booleans
	int are_plus_strands, are_minus_strands, is_small_dist, is_plus_order;
	int is_minus_order;
	
	// get the distance
	if (fwd_map->contigstart < rev_map->contigstart) {
		cs_fwd = fwd_map->contigstart; 
		cs_rev = rev_map->contigend;					
	} else {
		cs_fwd = fwd_map->contigend; 
		cs_rev = rev_map->contigstart;
	}
	dist = abs64(cs_fwd - cs_rev);
	is_small_dist = (dist < distcutoff);
	
	// get the strands
	s_fwd = fwd_map->strand;
	s_rev = rev_map->strand;
	are_plus_strands = ((s_fwd == s_rev) && (s_fwd == '+'));
	are_minus_strands = ((s_fwd == s_rev) && (s_fwd == '-'));
		
	// get the order
	is_plus_order = (are_plus_strands && (cs_rev < cs_fwd));
	is_minus_order = (are_minus_strands && (cs_fwd < cs_rev));
	
	// return distance if this is a good mapping
	if ((is_small_dist) &&  (is_plus_order || is_minus_order)) {
		return dist;
	} else {
		return 0;
	}
}

/* 
 * increment the mean, standard deviation, and histogram, using the
 * given distance
 */
void increments_stats(uint64_t good_mps_dist) {
	
	int binnr;
	
	// increase number of reads...
	gl_good_mps++;
	
	// mean
	double prev_mean = gl_mean;
	gl_mean = gl_mean + (good_mps_dist - gl_mean)/gl_good_mps;
	
	// stdev
	gl_stdev = gl_stdev + (good_mps_dist - prev_mean) * (good_mps_dist - gl_mean);
	
	// histogram
	if (do_gl_hist) {
		
		assert(good_mps_dist < distcutoff);
		assert(hist_distcutoff == distcutoff);
		binnr = (int) floor((good_mps_dist * 1.0 / hist_distcutoff) * HIST_BINS);
		assert(binnr < HIST_BINS);
		
		gl_hist[binnr] = gl_hist[binnr] + 1;
	
		if (debugmode) {
			if (gl_good_mps % 100000 == 0)
				fprintf(stderr, "DEBUGLINE: mean:%.2f, "
				    "stdev:%.2f, nr:%" PRIu64 ";"
				    "nr_called:%" PRIu64 "; binnr:%i; "
				    "glhist_binnr:%" PRIu64 " ||| dist:%"
				    PRIu64 "\n", 
				    gl_mean, sqrt(gl_stdev/gl_good_mps),
				    gl_good_mps, gl_debug_call, binnr,
				    gl_hist[binnr], good_mps_dist);
		}
	}
	
	if ((gl_mean_nr != 0) && (gl_good_mps >= gl_mean_nr) && fabs(prev_mean - gl_mean) < 1.0) {
		gl_done_mean = 1;
		
		
//		fprintf(stderr, "mean signal sent: gl_good_mps:%lli, gl_mean_nr:%lli, ABS:%i\n",
//				(long long int) gl_good_mps, (long long int) gl_mean_nr, 
//				ABS(prev_mean - gl_mean) < 0.2);
		
	}
}

/*
 * Compute the probabilities, given a mapping. Add it to the mate pair array
 * if the probabilities pass the thresholds.
 */
inline int add_p_stats(mapping_t * fwd_map, mapping_t * rev_map, 
		mate_pair_val *mp_set, double * totnormodds, int * mp_set_index ) {
	
	// if discordant are asked for and the contig names don't match
	if (!allow_diff_chr && 
			(strcmp(fwd_map->contigname, rev_map->contigname) != 0))
		return 0;
	
	
	
	// contig start and end
	uint64_t cs_fwd; // contigstart/end of fwd read
	uint64_t cs_rev; // contigstart/end of rev read
	
	// get the distance
	if (fwd_map->contigstart < rev_map->contigstart) {
		cs_fwd = fwd_map->contigstart; 
		cs_rev = rev_map->contigend;					
	} else {
		cs_fwd = fwd_map->contigend; 
		cs_rev = rev_map->contigstart;
	}
	uint64_t dist = abs64(cs_fwd - cs_rev);
	
	
	
	// pgenome
	double pgenome_fwd = fwd_map->pgenome;
	double pgenome_rev = rev_map->pgenome;
	int pgenome_bin = 0;
	double pgenome_cumsum = 0;
	double pgenome = 0;
		
	if (discordant) {
		pgenome = pgenome_fwd * pgenome_rev;
	} else {
		pgenome_bin = (int) floor((fabs((double)(dist) - gl_mean) 
				* 1.0 / hist_distcutoff) * HIST_BINS);
	
		if (pgenome_bin >= HIST_BINS)
			pgenome_cumsum = 0;
		else
			pgenome_cumsum = gl_hist_cumsum[pgenome_bin];
	
		pgenome = pgenome_fwd * pgenome_rev * pgenome_cumsum;
	}
	
	pgenome = MIN(ALMOST_ONE, pgenome);
	
	if (pgenome < pgenome_cutoff)
		return 0;

	//fprintf(stderr, "pgenome1: %f, pgenome2: %f, pgenome: %f, pgenome_bin: %i,,small_fabs:%f pgenome_cumsum:%f\n", 
		//		pgenome_fwd, pgenome_rev, pgenome, pgenome_bin, fabs((double)(dist) - gl_mean), pgenome_cumsum);

	
	// pchance
	double pchance;
		
	// pchance
	double pchance_fwd = fwd_map->pchance;
	double pchance_rev = rev_map->pchance;
	
	if (discordant || quickmode) {
		pchance = pchance_fwd * pchance_rev;
	} else { 
		double pchance_fwd_alt = 1 - pow(1 - pchance_fwd, 
				fabs((double)(dist) - gl_mean + 1) * 1.0 / (double)genome_length);
		double pchance_rev_alt = 1 - pow(1 - pchance_rev, 
				fabs((double)(dist) - gl_mean + 1) * 1.0 / (double)genome_length);
		pchance = (pchance_fwd * pchance_rev_alt + 
				pchance_rev * pchance_fwd_alt) / 2;
	}
	pchance = MAX(ALMOST_ZERO, pchance);
	
	if (pchance > pchance_cutoff)
		return 0;
	

	
	
	// assigning the values
	mp_set[*mp_set_index].fwd_rs = fwd_map;
	mp_set[*mp_set_index].rev_rs = rev_map;
	mp_set[*mp_set_index].pchance = pchance;
	mp_set[*mp_set_index].pgenome = pgenome;
	mp_set[*mp_set_index].normodds = pgenome/pchance;
	mp_set[*mp_set_index].dist = dist;
	mp_set_index[0]++;
	
	totnormodds[0] += pgenome/pchance;
	
	return 1;
}


/*
 * Source of fuction: SHRiMP (original code probably by Steve Rumble)
 * 
 * Return a string on the stack corresponding to an unsigned integer that also
 * features commas. E.g.: int 1000 yields "1,000".
 *
 * There is a pool of sequentially allocated buffers returned, so this should
 * be safe to use multiple times in function arguments.
 */
char * comma_integer(uint64_t val) {
	static char rets[50][32];	// no malloc, allow uses in fn args, etc
	static int col = 0;

	char *ret = rets[(col++ % (sizeof(rets) / sizeof(rets[0])))];
	char str[sizeof(rets[0])];
	int skip, i, j;

	memset(str, 0, sizeof(str));	// XXX - shut up, valgrind
	snprintf(str, sizeof(str), "%lli" , (long long int) val);

	skip = 3 - (strlen(str) % 3);
	for (i = j = 0; str[i] != '\0'; i++) {
		if ((i + skip) % 3 == 0 && i != 0)
			ret[j++] = ',';
		ret[j++] = str[i];
	}
	ret[j] = '\0';

	return (ret);
}

/*
 * Print a dashed line of 80 chars, including the newline char at the end
 */
void print_dashed_line() {
	fprintf(stderr, "----------------------------------------");
	fprintf(stderr, "----------------------------------------\n");
}


/*
 * Read a line from the mapping file, and assign the apropriate values to 
 * the mapping_t mapping.
 * The mappingfile can be binary or ascii, with input_file_type appropriately 
 * set.  
 */
int read_probcalc_line(FILE *mappingfile, mapping_t * mapping) {
	
	if (input_file_type == BINARY) {
		return fread (mapping, sizeof(mapping_t), 1, mappingfile);
		
	} else {
		
		char line[MAXLINELEN];
		int linelen = 0;
		int i;
		
		char ech; 	// temp character
		char field[256];
		int fieldnr = 1;
		int fieldindex = 0;
		
		// look though all of the tab-delimited fields
		while (fgets(line, MAXLINELEN, mappingfile) != NULL) {
			if (line[0] == '#')
				continue;
			
			linelen = strlen(line);
			for (i = 0; i < linelen; i++) {
				ech = line[i];
					
				if (ech == '\t') {
					field[fieldindex] = '\0';
					
					switch (fieldnr) {
						case 1:
							strcpy(mapping->readname, field);
							break;
							
						case 2:
							strcpy(mapping->contigname, field);
							break;
							
						case 3:
							if (strlen(field) != 1) {
								fprintf(stderr, "Error, strand is not single character: [len:%zu] %s\n", strlen(field), field);
								exit(1);
							}
							mapping->strand = field[0];
							break;
							
						case 4:
							mapping->contigstart = atoll(field);
							break;
							
						case 5:
							mapping->contigend = atoll(field);
							break;
							
						case 6:
							mapping->readstart = atoll(field);
							break;
							
						case 7:
							mapping->readend = atoll(field);
							break;
							
						case 8:
							mapping->readlength = atoi(field);
							break;
							
						case 9:
							mapping->score = atoi(field);
							break;
							
						case 10:
							strcpy(mapping->editstring, field);
							break;
							
						case 11:
							if (!Rflag)
								mapping->normodds = atof(field);
							break;
							
						case 12:
							if (Rflag)
								mapping->normodds = atof(field);
							else
								mapping->pgenome = atof(field);
							break;
							
						case 13:
							if (Rflag)
								mapping->pgenome = atof(field);
							else
								mapping->pchance = atof(field);
							break;
							
						case 14:
							if (Rflag)
								mapping->pchance = atof(field);
							else {
								fprintf(stderr, "no R Flag, and too many fields. line:\n");
								fprintf(stderr, "%s\n", line);
								exit(1);
							}
							break;
							
						default:
							fprintf(stderr, "error: failed to parse file line [%s]"
									" fieldnr:%i\n", line, fieldnr);
							break;
							
					}
					
					fieldnr++;
					fieldindex = 0;
					
				} else {
					field[fieldindex++] = ech;
				}
			} // for loop through the line
			
			return 1;
		} // while loop skipping over comment lines
		return 0; // since the result was therefore NULL
		
	}

	// avoid unreachable compiler warning; we'd get a no return warning if it were
	//assert(0);
	//return -1;
}

