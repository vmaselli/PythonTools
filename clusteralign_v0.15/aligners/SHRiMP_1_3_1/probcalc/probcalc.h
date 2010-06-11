/*	$Id: probcalc.h 383 2009-09-30 22:47:56Z matei $	*/

#define DEF_PCHANCE_CUTOFF	0.001
#define DEF_PGENOME_CUTOFF	0.0
#define DEF_NORMODDS_CUTOFF	0.0
#define DEF_TOP_MATCHES		10
#define DEF_NUMBER_MATCHES	10


/* Stats stuff */
double maxCount(int ins, int dels, int len, double delev, double deln, double insev, double insn); /* max indel Z */
double minCount(int ins, int dels, int len, double delev, double deln, double insev, double insn); /* min indel Z */
double subCount(int subs, int len); /* subs Z */
double fastchoose(int n, int m); /* fast choose function: uses fastlchoose */
double fastlchoose(int n, int m); /* fast log choose function: uses lgamma */
double fastfact(int n);
void initStats(int maxlen) ; /* initiate statistics - build the lookup tables */
void editstr_to_stats(char * str, int slen, int readlen, int * delFreq, int * insFreq);
void readIndelStats(int readlen, char * editstr, double * delev, double * insev, double * deln, double * insn);
