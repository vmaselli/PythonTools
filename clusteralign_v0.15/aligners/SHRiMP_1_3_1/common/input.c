/*	$Id: input.c 383 2009-09-30 22:47:56Z matei $	*/

#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "../common/sw-full-common.h"
#include "../common/input.h"
#include "../common/util.h"

enum {
	F_SKIP,
	F_READNAME,
	F_CONTIGNAME,
	F_STRAND,
	F_CONTIGSTART,
	F_CONTIGEND,
	F_READSTART,
	F_READEND,
	F_READLENGTH,
	F_SCORE,
	F_EDITSTRING,
	F_NORMODDS,
	F_PGENOME,
	F_PCHANCE,
	F_READSEQUENCE,
	F_UNKNOWN
};

struct _field_table {
	char const *field;
	int   type;
} field_table[] = {
	{ "#FORMAT:",		F_SKIP		},
	{ "readname",		F_READNAME	},
	{ "contigname",		F_CONTIGNAME	},
	{ "strand",		F_STRAND	},
	{ "contigstart",	F_CONTIGSTART	},
	{ "contigend",		F_CONTIGEND	},
	{ "readstart",		F_READSTART	},
	{ "readend",		F_READEND	},
	{ "readlength",		F_READLENGTH	},
	{ "score",		F_SCORE		},
	{ "editstring",		F_EDITSTRING	},
	{ "normodds",		F_NORMODDS	},
	{ "pgenome",		F_PGENOME	},
	{ "pchance",		F_PCHANCE	},
	{ "readsequence",	F_READSEQUENCE	},
	{ NULL,			0		}
};

struct format_spec {
	int nfields;
	int *fields;
};

bool
editstr_to_sfr(const char *editstr, struct sw_full_results *sfrp)
{
	char *scratch;
	size_t len, i, j;
	bool inparen = false;

	len = strlen(editstr);
	scratch = (char *)xmalloc(len + 1);
	scratch[0] = '\0';

	memset(sfrp, 0, sizeof(*sfrp));

	for (i = j = 0; i <= len; i++) {
		assert(j <= len);

		if (inparen) {
			if (editstr[i] == ')') {
				scratch[j] = '\0';
				if (strcspn(scratch, "ACGTUMRWSYKVHDBXNacgtumrwsykvhdbxn") != 0) {
					free(scratch);
					return (false);
				}
				if (shrimp_mode == MODE_COLOUR_SPACE)
					assert(strspn(scratch, "acgtumrwsykvhdbxn") == 0);
				sfrp->deletions += strlen(scratch);
				inparen = false;
				j = 0;
				scratch[0] = '\0';
			} else {
				scratch[j++] = editstr[i];
			}
		} else {
			if (!isdigit((int)editstr[i]) && j != 0) {
				scratch[j] = '\0';
				sfrp->matches += strtoul(scratch, NULL, 0);
				j = 0;
				scratch[0] = '\0';
			}

			switch (editstr[i]) {
			case '-':
				sfrp->insertions++;
				break;
			case '(':
				inparen = true;
				break;
			case '\0':
				break;
			case 'X':
			case 'x':
				sfrp->crossovers++;
				break;
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'N':
			case 'a':
			case 'c':
			case 'g':
			case 't':
			case 'n':
				sfrp->mismatches++;
				break;
			default:
				if (isdigit((int)editstr[i]))
					scratch[j++] = editstr[i];
				else
					return (false);
			}
		}
	}

	free(scratch);

	return (!inparen);
}

static struct format_spec * 
format_get_default()
{
	struct format_spec *fsp;
	int i;

	fsp = (struct format_spec *)xmalloc(sizeof(*fsp));
	memset(fsp, 0, sizeof(*fsp));
	fsp->nfields = sizeof(field_table)/sizeof(field_table[0]) - 2;
	fsp->fields = (int *)xmalloc(sizeof(*fsp->fields) * fsp->nfields);

	for (i = 0; i < fsp->nfields; i++)
		fsp->fields[i] = field_table[i + 1].type;

	return (fsp);
}

static struct format_spec *
format_get_from_string(char *format)
{
	struct format_spec *fsp;
	char *field;
	int i, next;

	fsp = (struct format_spec *)xmalloc(sizeof(*fsp));
	memset(fsp, 0, sizeof(*fsp));

	field = strtok(format, " ");
	while (field != NULL) {
		field = strtrim(field);

		next = F_UNKNOWN;
		for (i = 0; field_table[i].field != NULL; i++) {
			if (strcmp(field_table[i].field, field) == 0) {
				next = field_table[i].type;
				break;
			}
		}

		if (next != F_SKIP) {
			fsp->nfields++;
			fsp->fields = (int *)xrealloc(fsp->fields, sizeof(*fsp->fields) * fsp->nfields);
			fsp->fields[fsp->nfields-1] = next;
		}

		if (next == F_UNKNOWN)
			fprintf(stderr, "warning: unknown format field [%s]\n", field);
		
		field = strtok(NULL, " ");
	}

	return (fsp);
}

static void
format_free(struct format_spec *fsp)
{

	free(fsp->fields);
	free(fsp);
}

static void
handle_field(struct input *inp, int type, char *val)
{
	struct sw_full_results sfr;

	assert(inp != NULL && val != NULL);

	if (is_whitespace(val))
		return;

	switch (type) {
	case F_READNAME:
		inp->read = xstrdup(val);
		break;
	case F_CONTIGNAME:
		inp->genome = xstrdup(val);
		break;
	case F_STRAND:
		if (*val == '-')
			inp->flags |= INPUT_FLAG_IS_REVCMPL;
		break;
	case F_CONTIGSTART:
		inp->genome_start = strtoul(val, NULL, 0) - 1;
		break;
	case F_CONTIGEND:
		inp->genome_end = strtoul(val, NULL, 0) - 1;
		break;
	case F_READSTART:
		inp->read_start = (uint16_t)(strtoul(val, NULL, 0) - 1);
		break;
	case F_READEND:
		inp->read_end = (uint16_t)(strtoul(val, NULL, 0) - 1);
		break;
	case F_SCORE:
		inp->score = strtoul(val, NULL, 0);
		break;
	case F_NORMODDS:
		inp->normodds = atof(val);
		inp->flags |= INPUT_FLAG_HAS_NORMODDS;
		break;
	case F_PGENOME:
		inp->pgenome = atof(val);
		inp->flags |= INPUT_FLAG_HAS_PGENOME;
		break;
	case F_PCHANCE:
		inp->pchance = atof(val);
		inp->flags |= INPUT_FLAG_HAS_PCHANCE;
		break;
	case F_READSEQUENCE:
		inp->read_seq = xstrdup(val);
		break;
	case F_READLENGTH:
		inp->read_length = (uint16_t)strtoul(val, NULL, 0);
		break;
	case F_EDITSTRING:
		inp->edit = xstrdup(val);
		editstr_to_sfr(val, &sfr);
		inp->matches = (uint16_t)sfr.matches;
		inp->mismatches = (uint16_t)sfr.mismatches;
		inp->deletions = (uint16_t)sfr.deletions;
		inp->insertions = (uint16_t)sfr.insertions;
		inp->crossovers = (uint16_t)sfr.crossovers;
		break;
	case F_UNKNOWN:
		break;
	default:
		fprintf(stderr, "warning: crazy field in input [%s]\n", val);
	}
}

void
input_free(struct input *inp)
{

	if (inp->read != NULL)
		free(inp->read);
	if (inp->read_seq != NULL)
		free(inp->read_seq);
	if (inp->genome != NULL)
		free(inp->genome);
	if (inp->edit != NULL)
		free(inp->edit);
}

/*
 * Parse the key-value paired output created by output_normal in output.c.
 * This will only parse lines beginning with '>', so it'll work just fine
 * with both rmapper's standard and pretty printed output formats.
 *
 * Returns false on EOF.
 */
bool
input_parseline(gzFile fp, struct input *inp)
{
	char buf[8192];
	struct format_spec *fsp;
	bool ret = true;

	assert(fp != NULL && inp != NULL);

	memset(inp, 0, sizeof(*inp));
	fsp = format_get_default();

	while (true) {
		/* XXX: Assume sizeof(buf) is always big enough. */
		if (fast_gzgets(fp, buf, sizeof(buf)) == NULL) {
			ret = false;
			break;
		}

		if (buf[0] == '#' && strncmp(buf, "#FORMAT:", 8) == 0) {
			format_free(fsp);
			fsp = format_get_from_string(buf);
		} else if (buf[0] == '>') {
			char *str = &buf[1];
			char *val;
			int i;

			val = strtok(str, "\t");
			for (i = 0; val != NULL; i++) {
				if (i < fsp->nfields)
					handle_field(inp, fsp->fields[i], val);

				val = strtok(NULL, "\t");
			}

			break;
		}
	}

	format_free(fsp);

	return (ret);
}
