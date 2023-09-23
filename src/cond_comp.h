#ifdef _OPENMP
#define _OMPTGT_(x) $OMP x
#ifdef _OMPOFFLOADTGT
#define _OMPOFFLOADTGT_(x) $OMP x
#else
#define _OMPOFFLOADTGT_(x) disabled
#endif
#else
#define _OMPTGT_(x) disabled
#define _OMPOFFLOADTGT_(x) disabled
#endif
