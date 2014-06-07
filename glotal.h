#ifndef __GLOTTAL_H__
#define __GLOTTAL_H__
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "rtwtypes.h"
#include "glottal_types.h"
#define NUMERO_PARAMETROS_ANALISIS		46
#define NUMERO_MAXIMO_MUESTRAS_VOZ		200000
#define NUMERO_MAXIMO_TRAMAS_ANALISIS	1000
#define NUMERO_MAXIMO_PUNTOS_FFT		16384
#define LONGITUD_CADENA_TEXTO			60

extern void glottal(int32_T long_si, const int16_T si[], int32_T sf, int32_T left_lim, int32_T righ_lim, int8_T sw22,
	int8_T *exec_status, TDimensiones *Dimensiones, int8_T out_ban[],
	real_T OutArr[][NUMERO_PARAMETROS_ANALISIS], char_T OutNam[][LONGITUD_CADENA_TEXTO],
	real_T par_vec_med[], real_T par_vec_std[], real_T ResultsPointRef[], real_T GloSourceRefPoint[], real_T GloFlowRefPoint[],
	real_T ugn_lvl[], int32_T *ref_pnt, int32_T arg_min[], real_T *dist_vfolds);
#endif
