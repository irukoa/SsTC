#ifdef _OPENMP
#define using_OMP  .true.
#define _OMPTGT_(x) $OMP x
#else
#define using_OMP  .false.
#define _OMPTGT_(x) disabled
#endif

#ifdef _OPENACC
#define using_ACC  .true.
#define _ACCTGT_(x) $ACC x
#else
#define using_ACC  .false.
#define _ACCTGT_(x) disabled
#endif