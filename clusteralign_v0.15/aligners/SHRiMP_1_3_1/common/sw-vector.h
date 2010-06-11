/*	$Id: sw-vector.h 348 2009-06-16 23:26:27Z rumble $	*/

int	sw_vector_setup(int, int, int, int, int, int, int, int, int, bool);
void	sw_vector_stats(uint64_t *, uint64_t *, uint64_t *, double *);
int	sw_vector(uint32_t *, int, int, uint32_t *, int, uint32_t *, int, bool);
