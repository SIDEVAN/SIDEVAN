#ifndef __FFT_H__
#define __FFT_H__

#include "rtwtypes.h"
void FastFourierTransform(int32_T longitudTrama, const creal_T trama[], creal_T salidaFFT[]);

void FFT_32(const creal_T x[], creal_T y[32]);
void FFT_64(const creal_T x[], creal_T y[64]);
void FFT_128(const creal_T x[], creal_T y[128]);
void FFT_256(const creal_T x[], creal_T y[256]);
void FFT_512(const creal_T x[], creal_T y[512]);
void FFT_1024(const creal_T x[], creal_T y[1024]);
void FFT_2048(const creal_T x[], creal_T y[2048]);
void FFT_4096(const creal_T x[], creal_T y[4096]);
void FFT_8192(const creal_T x[], creal_T y[8192]);
void FFT_16384(const creal_T x[], creal_T y[16384]);
#endif
