#ifndef __DBTYPES_H__
#define __DBTYPES_H__

/* Read mapping struct */
#define READNAME_LEN    32
#define CONTIGNAME_LEN  32
#define EDITSTRING_LEN  32

typedef struct {
    char        readname[READNAME_LEN];
    char        contigname[CONTIGNAME_LEN];
    uint64_t    contigstart;
    uint64_t    contigend;
    char        strand;
    uint8_t     readstart;
    uint8_t     readend;
    uint8_t     readlength;
    uint8_t     score;
    char        editstring[EDITSTRING_LEN];
    double      normodds;
    double      pgenome;
    double      pchance;
} mapping_t;

/* Read-hit index struct */
typedef struct {
    char        readname[READNAME_LEN];
    uint64_t    offset;
} readname_idx_t;
            
#endif /* !defined(__DBTYPES_H__) */
