#ifdef _OPENMP
#define using_OMP  .true.
#define _OMPTGT_(x) $OMP x
#else
#define using_OMP  .false.
#define _OMPTGT_(x) disabled
#endif
