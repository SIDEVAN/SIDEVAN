#ifndef __GLOTTAL_TYPES_H__
#define __GLOTTAL_TYPES_H__

#define	FALSE	0
#define	TRUE	1

#include <limits.h>

typedef unsigned char boolean_T;
typedef char char_T;
typedef signed char int8_T;
typedef unsigned char uint8_T;
typedef short int16_T;
typedef unsigned short uint16_T;
typedef int int32_T;
typedef unsigned int uint32_T;
typedef float real32_T;
typedef double real64_T;

typedef double real_T;


#ifndef typedef_creal_T
#define typedef_creal_T
typedef struct
{  
	real_T re;  
	real_T im;  
} creal_T;  
#endif

#ifndef typedef_TDimensiones
#define typedef_TDimensiones
typedef struct
{
    int32_T n_Param;
    int32_T n_OutFlags;
    int32_T n_ResultsPointRef;
    int32_T n_GloSourcePoints;
    int32_T n_GlottalSourceLevd;
    int32_T n_MinArg;
} TDimensiones;
#endif /*typedef_TDimensiones*/

#endif
 
