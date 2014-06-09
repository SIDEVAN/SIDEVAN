#include <windows.h>
#include <stdio.h>
#include "glottal_types.h"
#include "glottal.h"


#define _CRT_SECURE_NO_WARNINGS
#pragma comment(lib,"user32.lib")


__declspec(dllexport) void __cdecl sidevan_parm(int32_T long_si, const int16_T si[], int32_T sf, int32_T 
		  left_lim, int32_T righ_lim, int8_T sw22, int8_T *exec_status,
                  TDimensiones *Dimensiones, int8_T out_ban[], real_T
                  OutArr[][NUMERO_PARAMETROS_ANALISIS], char_T OutNam[][LONGITUD_CADENA_TEXTO], real_T par_vec_med[],
                  real_T par_vec_std[], real_T ResultsPointRef[], real_T
                  GloSourceRefPoint[], real_T GloFlowRefPoint[],
                  real_T ugn_lvl[], int32_T *ref_pnt, int32_T arg_min[],
                  real_T *dist_vfolds);

void sidevan_parm(int32_T long_si, const int16_T si[], int32_T sf, int32_T
                  left_lim, int32_T righ_lim, int8_T sw22, int8_T *exec_status,
                  TDimensiones *Dimensiones, int8_T out_ban[], real_T
                  OutArr[][NUMERO_PARAMETROS_ANALISIS], char_T OutNam[][LONGITUD_CADENA_TEXTO], real_T par_vec_med[],
                  real_T par_vec_std[], real_T ResultsPointRef[], real_T
                  GloSourceRefPoint[], real_T GloFlowRefPoint[],
                  real_T ugn_lvl[], int32_T *ref_pnt, int32_T arg_min[],
                  real_T *dist_vfolds)
{
	glottal(long_si, si, sf, left_lim, righ_lim, sw22, exec_status, Dimensiones, out_ban, OutArr, OutNam, par_vec_med, par_vec_std, ResultsPointRef, GloSourceRefPoint, GloFlowRefPoint, ugn_lvl, ref_pnt, arg_min, dist_vfolds);
}

