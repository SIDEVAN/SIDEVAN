#include "glottal.h"
#include "fft.h"

#define NUMERO_MAXIMO_RASGOS_CEPSTRALES		16
#define NUMERO_MAXIMO_ETAPAS_CELOSIA		128


static void CalcularModulosFFT(int32_T longitudTrama, const real_T trama[],
  real_T modulosFFT[]);
static void CalcularModulosFFTConRelleno(int32_T longitudFFT, int32_T
  longitudTrama, const real_T trama[], real_T modulosFFT[]);
static void CalcularModulosLogaritmicosCompletosFFT(int32_T longitudTrama, const
  real_T trama[], real_T modulosFFT[]);
static void CalcularModulosLogaritmicosCompletosFFT_VentanaIntervalo(int32_T
  longitudFFT, const real_T ventana[], int32_T inicioTrama, int32_T
  finTrama, const real_T trama[], real_T modulosFFT_log[]);
static void CalcularModulosLogaritmicosFFT(int32_T longitudFFT, int32_T
  longitudTrama, const real_T trama[], real_T modulosFFT[], real_T
  modulosFFT_log[]);
static void CambiarSignoSerie(int32_T numeroElementos, const real_T serie[],
  real_T serieSignoCambiado[]);
static void ConfeccionarListaParametrosSalida(char_T OutNam[][LONGITUD_CADENA_TEXTO]);
static void ConvolucionarSeries(int32_T n1, const real_T vector1[],
  int32_T n2, const real_T vector2[], real_T vectorConvolucion[]);
static void CopiarFragmentoVectorEntero(int32_T inicioIntervalo, int32_T finIntervalo,
  const int32_T vectorOrigen[], int32_T vectorDestino[]);
static void CopiarFragmentoVector(int32_T inicioIntervalo, int32_T
  finIntervalo, const real_T vectorOrigen[], real_T vectorDestino[]);
static void CopiarFragmentoVectorSinMediaConversionDesdeEntero(int32_T
  inicioIntervalo, int32_T finIntervalo, int32_T longitudVectorOrigen, const
  int16_T vectorOrigen[], real_T vectorDestino[]);
static void CopiarVector(int32_T longitudVector, const real_T vectorOrigen
  [], real_T vectorDestino[]);
static void CorregirSignoSerie(int32_T numeroElementos, const real_T serie
  [], real_T serieCorregida[]);
static int32_T DivisionEntera(int32_T dividendo, int32_T divisor);
static void EliminarMediaSerie(int32_T numeroElementos, const real_T serie
  [], real_T serieSinMedia[]);
static creal_T ExponencialComplejaSoloParteImaginaria(real_T z_imag);
static void FFT_radix2(int32_T longitudTrama, const creal_T trama[],
  creal_T salidaFFT[]);
static void FFT_radix_no2(int32_T longitudTrama, const real_T trama[],
  int32_T numeroPuntosEnFrecuencia, creal_T salidaFFT[]);
static int32_T IndicePrimerValorMayorEntero(int32_T numeroElementos, const int32_T
  serie[], int32_T valorReferencia);
static int32_T IndicePrimerValorMayor(int32_T numeroElementos, const real_T
  serie[], real_T valorReferencia);
static int32_T IndiceUltimoValorMenorEntero(int32_T numeroElementos, const int32_T
  serie[], int32_T valorReferencia);
static int32_T IndiceUltimoValorMenor(int32_T numeroElementos, const real_T
  serie[], real_T valorReferencia);
static real_T LogaritmoBase10(real_T numero);
static void MaximoIntervaloSerie(int32_T inicioIntervalo, int32_T finIntervalo,
  const real_T serie[], int32_T *indice, real_T *maximo);
static void MaximoSerie(int32_T numeroElementos, const real_T serie[], int32_T *
  indice, real_T *maximo);
static void MaximoSerieEntera(int32_T numeroElementos, const int32_T serie[],
  int32_T *indice, int32_T *maximo);
static void MaximoValorAbsolutoIntervaloSerie(int32_T inicioIntervalo, int32_T
  finIntervalo, const real_T serie[], int32_T *indice, real_T *maximo);
static void MinimoIntervaloSerie(int32_T inicioIntervalo, int32_T finIntervalo,
  const real_T serie[], int32_T *indice, real_T *minimo);
static void MinimoSerie(int32_T numeroElementos, const real_T serie[],
  int32_T *indice, real_T *minimo);
static void MinimoSerieEmpezandoPorFinal(int32_T numeroElementos, const real_T
  serie[], int32_T *indice, real_T *minimo);
static void ObtenerFraccionesSenoPi(int32_T N, real_T serie[]);
static void ObtenerFragmentoTrama(int32_T longitudFFT, int32_T inicioTrama,
  int32_T finTrama, const real_T trama[], real_T tramaFFT[]);
static real_T ObtenerMediaSerieEnteraCorta(int32_T numeroElementos, const int16_T serie[]);
static real_T ObtenerMediaSerieEntera(int32_T numeroElementos, const int32_T serie[]);
static real_T ObtenerMediaSerie(int32_T numeroElementos, const real_T serie[]);
static real_T ObtenerMediaSerieIntervalo(int32_T inicio, int32_T fin, const
  real_T serie[]);
static void ObtenerMediayDesviacionEstandarEntera(int32_T numeroElementos, const
  int32_T serie[], real_T *media, real_T *desviacion);
static void ObtenerMediayDesviacionEstandar(int32_T numeroElementos, const
  real_T serie[], real_T *media, real_T *desviacion);
static void ObtenerMinimoMaximoSerie(int32_T numeroElementos, const real_T
  serie[], real_T *minimo, real_T *maximo);
static void ObtenerSeparacionesComienzoSerie(int32_T numeroElementos, const
  int32_T serie[], int32_T lista[]);
static int32_T ObtenerSiguientePotencia2(int32_T numero);
static real_T ObtenerSumaValoresCuadradoIntervaloSerie(int32_T inicioIntervalo,
  int32_T finIntervalo, const real_T serie[]);
static real_T ObtenerSumaValoresIntervaloSerie(int32_T inicioIntervalo, int32_T
  finIntervalo, const real_T serie[]);
static void ObtenerTramaEnventanada(int32_T longitudFFT, const real_T ventana
  [], int32_T inicioTrama, int32_T finTrama, const real_T trama[],
  real_T tramaEnventanada[]);
static int32_T ParteEntera(real_T numero);
static int32_T Redondear(real_T numero);
static int8_T Signo(real_T numero);
static int32_T ValorMayorEntero(int32_T numero1, int32_T numero2);
static real_T ValorMayor(real_T numero1, real_T numero2);
static int32_T ValorMenor(int32_T numero1, int32_T numero2);
static void ValoresAbsolutosSerie(int32_T numeroElementos, const real_T serie
  [], real_T serieCorregida[]);
static void VentanaHamming(int32_T N, real_T hamming[]);
static void cel2(int32_T long_s, const real_T s[], const real_T r[],
                   int32_T K, real_T u[], real_T c[]);
static void cel1(int32_T long_s, const real_T s[], int32_T K, real_T x[]);
static real_T cepstralpitch(int32_T long_InputVec, const real_T InputVector
  [], int32_T SamplingFrequency, int8_T SwitchSignal);
static void det_sing_mw4(int32_T long_mwpsd, const real_T mwpsd[], real_T
  *val_end, int32_T *pos_end, real_T mx_1[], int32_T pos_mx_1[], real_T
  mn_1[], int32_T pos_mn_1[], real_T mx_2[], int32_T pos_mx_2[],
  real_T nsf[], char_T *rel_ban);
static void glosotime(int32_T long_GloSource, const real_T GloSource[],
                      real_T MeanGloSource, int32_T sf, real_T *tR1r, real_T
                      *tR2r, real_T *tO1r, real_T *tO2r, real_T *tMr, real_T
                      *AR1r, real_T *AR2r, real_T *AO1r, real_T *AO2r, real_T
                      *RT, real_T *OT, real_T *CT, real_T GloFlow[]);
static void glottalclipping(int32_T long_GloSource, const real_T GloSource
  [], real_T SelectThreshold, real_T ExpectedPitch, int32_T *ResNoMinima,
  int32_T ResGloSourceMinArg[]);
static void highpass(int32_T long_Vi, const real_T Vi[], int32_T fils,
                     real_T Vo[]);
static void highpassfast(int32_T long_inp_vec, const real_T inp_vec[],
  int32_T tamavpas, real_T out_vec[]);
static void integ(int32_T long_inp_vec, const real_T inp_vec[], real_T
                  leakage, int32_T fl, real_T out_vec[]);
static void lowpass(int32_T long_Vi, const real_T Vi[], int32_T fils,
                    real_T Vo[]);
static void maxmin1d(int32_T long_Fun, const real_T Fun[], int8_T mxmn,
                     int32_T *dim, real_T val_min_sal[], int32_T
                     arg_min_sal[]);
static void unbias(int32_T long_ugn, const real_T ugn[], int32_T num_min,
                   const int32_T arg_min[], real_T ugndbl[]);
static void zeroxing(int32_T long_InputVector, const real_T InputVector[],
                     int32_T ZeroXingPoints[], int32_T *NoZeroXings);


static void CalcularModulosFFT(int32_T longitudTrama, const real_T trama[], real_T modulosFFT[])
{
  creal_T salidaFFT[NUMERO_MAXIMO_PUNTOS_FFT];
  int32_T nModulos;
  int32_T i;

  FFT_radix_no2(longitudTrama, trama, longitudTrama, salidaFFT);
  nModulos = DivisionEntera(longitudTrama, 2);
  for(i=0;i< nModulos; i++)
		modulosFFT[i]= sqrt(salidaFFT[i].re*salidaFFT[i].re + salidaFFT[i].im*salidaFFT[i].im);
}

static void CalcularModulosFFTConRelleno(int32_T longitudFFT, int32_T longitudTrama, const real_T trama[], real_T modulosFFT[])
{
  real_T tramaFFT[NUMERO_MAXIMO_PUNTOS_FFT];
  creal_T salidaFFT[NUMERO_MAXIMO_PUNTOS_FFT];
  int32_T nModulos;
  int32_T i;

  ObtenerFragmentoTrama(longitudFFT, 1, longitudTrama, trama, tramaFFT);
  FFT_radix_no2(longitudFFT, tramaFFT, longitudFFT, salidaFFT);
  nModulos = DivisionEntera(longitudFFT, 2);
  for(i=0;i< nModulos; i++)
    modulosFFT[i] = sqrt(salidaFFT[i].re*salidaFFT[i].re + salidaFFT[i].im*salidaFFT[i].im);
}

static void CalcularModulosLogaritmicosCompletosFFT(int32_T longitudTrama, const real_T trama[], real_T modulosFFT[])
{
  creal_T salidaFFT[NUMERO_MAXIMO_PUNTOS_FFT];
  int32_T i;

  FFT_radix_no2(longitudTrama, trama, longitudTrama, salidaFFT);
  for(i=0;i< longitudTrama; i++)
    modulosFFT[i] = 10.0 * LogaritmoBase10(salidaFFT[i].re*salidaFFT[i].re + salidaFFT[i].im*salidaFFT[i].im);
}

static void CalcularModulosLogaritmicosCompletosFFT_VentanaIntervalo(int32_T
  longitudFFT, const real_T ventana[], int32_T inicioTrama, int32_T
  finTrama, const real_T trama[], real_T modulosFFT_log[])
{
  real_T tramaFFT[NUMERO_MAXIMO_PUNTOS_FFT];
  creal_T salidaFFT[NUMERO_MAXIMO_PUNTOS_FFT];
  int32_T i;

  ObtenerTramaEnventanada(longitudFFT, ventana, inicioTrama, finTrama, trama, tramaFFT);
  FFT_radix_no2(longitudFFT, tramaFFT, longitudFFT, salidaFFT);
  for(i=0;i< longitudFFT; i++)
    modulosFFT_log[i] = 10.0 * LogaritmoBase10(salidaFFT[i].re*salidaFFT[i].re + salidaFFT[i].im*salidaFFT[i].im);
}

static void CalcularModulosLogaritmicosFFT(int32_T longitudFFT, int32_T
  longitudTrama, const real_T trama[], real_T modulosFFT[], real_T modulosFFT_log[])
{
  creal_T salidaFFT[NUMERO_MAXIMO_PUNTOS_FFT];
  int32_T nModulos;
  int32_T i;
  real_T modulo2;

  FFT_radix_no2(longitudTrama, trama, longitudFFT, salidaFFT);
  nModulos = DivisionEntera(longitudFFT, 2);
  for(i=0;i<nModulos;i++)
  {
	modulo2 = salidaFFT[i].re*salidaFFT[i].re + salidaFFT[i].im*salidaFFT[i].im;
    modulosFFT[i] = sqrt(modulo2);
    modulosFFT_log[i] = 10.0* LogaritmoBase10(modulo2);
  }
}

static void CambiarSignoSerie(int32_T numeroElementos, const real_T serie[], real_T serieSignoCambiado[])
{
  int32_T i;

  for(i=0;i<numeroElementos;i++)
    serieSignoCambiado[i] = -serie[i];
}

static void ConfeccionarListaParametrosSalida(char_T OutNam[][LONGITUD_CADENA_TEXTO])
{
  const char_T nombreParametro_01[]= " 1. Absolute Pitch           ";
  const char_T nombreParametro_02[]= " 2. Abs. Norm. Jitter        ";
  const char_T nombreParametro_03[]= " 3. Abs. Norm. Ar. Shimmer   ";
  const char_T nombreParametro_04[]= " 4. Abs. Norm. Min. Sharp.   ";
  const char_T nombreParametro_05[]= " 5. Noise-Harm. Ratio (NHR)  ";
  const char_T nombreParametro_06[]= " 6. Muc./AvAc. Energy (MAE)  ";
  const char_T nombreParametro_07[]= " 7. MWC Cepstral 1           ";
  const char_T nombreParametro_08[]= " 8. MWC Cepstral 2           ";
  const char_T nombreParametro_09[]= " 9. MWC Cepstral 3           ";
  const char_T nombreParametro_10[]= "10. MWC Cepstral 4           ";
  const char_T nombreParametro_11[]= "11. MWC Cepstral 5           ";
  const char_T nombreParametro_12[]= "12. MWC Cepstral 6           ";
  const char_T nombreParametro_13[]= "13. MWC Cepstral 7           ";
  const char_T nombreParametro_14[]= "14. MWC Cepstral 8           ";
  const char_T nombreParametro_15[]= "15. MWC Cepstral 9           ";
  const char_T nombreParametro_16[]= "16. MWC Cepstral 10          ";
  const char_T nombreParametro_17[]= "17. MWC Cepstral 11          ";
  const char_T nombreParametro_18[]= "18. MWC Cepstral 12          ";
  const char_T nombreParametro_19[]= "19. MWC Cepstral 13          ";
  const char_T nombreParametro_20[]= "20. MWC Cepstral 14          ";
  const char_T nombreParametro_21[]= "21. MW PSD 1st Max. ABS.     ";
  const char_T nombreParametro_22[]= "22. MW PSD 1st Min. rel.     ";
  const char_T nombreParametro_23[]= "23. MW PSD 2nd Max. rel.     ";
  const char_T nombreParametro_24[]= "24. MW PSD 2nd Min. rel.     ";
  const char_T nombreParametro_25[]= "25. MW PSD 3rd Max. rel.     ";
  const char_T nombreParametro_26[]= "26. MW PSD End Val. rel.     ";
  const char_T nombreParametro_27[]= "27. MW PSD 1st Max. Pos. ABS.";
  const char_T nombreParametro_28[]= "28. MW PSD 1st Min. Pos. rel.";
  const char_T nombreParametro_29[]= "29. MW PSD 2nd Max. Pos. rel.";
  const char_T nombreParametro_30[]= "30. MW PSD 2nd Min. Pos. rel.";
  const char_T nombreParametro_31[]= "31. MW PSD 3rd Max. Pos. rel.";
  const char_T nombreParametro_32[]= "32. MW PSD End Val. Pos. rel.";
  const char_T nombreParametro_33[]= "33. MW PSD 1st Min NSF       ";
  const char_T nombreParametro_34[]= "34. MW PSD 2nd Min NSF       ";
  const char_T nombreParametro_35[]= "35. Rel. Recov. 1 Time       ";
  const char_T nombreParametro_36[]= "36. Rel. Recov. 2 Time       ";
  const char_T nombreParametro_37[]= "37. Rel. Open 1 Time         ";
  const char_T nombreParametro_38[]= "38. Rel. Open 2 Time         ";
  const char_T nombreParametro_39[]= "39. Rel. Max. Ampl. Time     ";
  const char_T nombreParametro_40[]= "40. Rel. Recov. 1 Ampl.      ";
  const char_T nombreParametro_41[]= "41. Rel. Recov. 2 Ampl.      ";
  const char_T nombreParametro_42[]= "42. Rel. Open 1 Ampl.        ";
  const char_T nombreParametro_43[]= "43. Rel. Open 2 Ampl.        ";
  const char_T nombreParametro_44[]= "44. Rel. Stop Flow Time      ";
  const char_T nombreParametro_45[]= "45. Rel. Start Flow Time     ";
  const char_T nombreParametro_46[]= "46. Rel. Closing Time        ";

  memcpy(&OutNam[0][0], &nombreParametro_01[0], sizeof(nombreParametro_01));
  memcpy(&OutNam[1][0], &nombreParametro_02[0], sizeof(nombreParametro_02));
  memcpy(&OutNam[2][0], &nombreParametro_03[0], sizeof(nombreParametro_03));
  memcpy(&OutNam[3][0], &nombreParametro_04[0], sizeof(nombreParametro_04));
  memcpy(&OutNam[4][0], &nombreParametro_05[0], sizeof(nombreParametro_05));
  memcpy(&OutNam[5][0], &nombreParametro_06[0], sizeof(nombreParametro_06));
  memcpy(&OutNam[6][0], &nombreParametro_07[0], sizeof(nombreParametro_07));
  memcpy(&OutNam[7][0], &nombreParametro_08[0], sizeof(nombreParametro_08));
  memcpy(&OutNam[8][0], &nombreParametro_09[0], sizeof(nombreParametro_09));
  memcpy(&OutNam[9][0], &nombreParametro_10[0], sizeof(nombreParametro_10));
  memcpy(&OutNam[10][0], &nombreParametro_11[0], sizeof(nombreParametro_11));
  memcpy(&OutNam[11][0], &nombreParametro_12[0], sizeof(nombreParametro_12));
  memcpy(&OutNam[12][0], &nombreParametro_13[0], sizeof(nombreParametro_13));
  memcpy(&OutNam[13][0], &nombreParametro_14[0], sizeof(nombreParametro_14));
  memcpy(&OutNam[14][0], &nombreParametro_15[0], sizeof(nombreParametro_15));
  memcpy(&OutNam[15][0], &nombreParametro_16[0], sizeof(nombreParametro_16));
  memcpy(&OutNam[16][0], &nombreParametro_17[0], sizeof(nombreParametro_17));
  memcpy(&OutNam[17][0], &nombreParametro_18[0], sizeof(nombreParametro_18));
  memcpy(&OutNam[18][0], &nombreParametro_19[0], sizeof(nombreParametro_19));
  memcpy(&OutNam[19][0], &nombreParametro_20[0], sizeof(nombreParametro_20));
  memcpy(&OutNam[20][0], &nombreParametro_21[0], sizeof(nombreParametro_21));
  memcpy(&OutNam[21][0], &nombreParametro_22[0], sizeof(nombreParametro_22));
  memcpy(&OutNam[22][0], &nombreParametro_23[0], sizeof(nombreParametro_23));
  memcpy(&OutNam[23][0], &nombreParametro_24[0], sizeof(nombreParametro_24));
  memcpy(&OutNam[24][0], &nombreParametro_25[0], sizeof(nombreParametro_25));
  memcpy(&OutNam[25][0], &nombreParametro_26[0], sizeof(nombreParametro_26));
  memcpy(&OutNam[26][0], &nombreParametro_27[0], sizeof(nombreParametro_27));
  memcpy(&OutNam[27][0], &nombreParametro_28[0], sizeof(nombreParametro_28));
  memcpy(&OutNam[28][0], &nombreParametro_29[0], sizeof(nombreParametro_29));
  memcpy(&OutNam[29][0], &nombreParametro_30[0], sizeof(nombreParametro_30));
  memcpy(&OutNam[30][0], &nombreParametro_31[0], sizeof(nombreParametro_31));
  memcpy(&OutNam[31][0], &nombreParametro_32[0], sizeof(nombreParametro_32));
  memcpy(&OutNam[32][0], &nombreParametro_33[0], sizeof(nombreParametro_33));
  memcpy(&OutNam[33][0], &nombreParametro_34[0], sizeof(nombreParametro_34));
  memcpy(&OutNam[34][0], &nombreParametro_35[0], sizeof(nombreParametro_35));
  memcpy(&OutNam[35][0], &nombreParametro_36[0], sizeof(nombreParametro_36));
  memcpy(&OutNam[36][0], &nombreParametro_37[0], sizeof(nombreParametro_37));
  memcpy(&OutNam[37][0], &nombreParametro_38[0], sizeof(nombreParametro_38));
  memcpy(&OutNam[38][0], &nombreParametro_39[0], sizeof(nombreParametro_39));
  memcpy(&OutNam[39][0], &nombreParametro_40[0], sizeof(nombreParametro_40));
  memcpy(&OutNam[40][0], &nombreParametro_41[0], sizeof(nombreParametro_41));
  memcpy(&OutNam[41][0], &nombreParametro_42[0], sizeof(nombreParametro_42));
  memcpy(&OutNam[42][0], &nombreParametro_43[0], sizeof(nombreParametro_43));
  memcpy(&OutNam[43][0], &nombreParametro_44[0], sizeof(nombreParametro_44));
  memcpy(&OutNam[44][0], &nombreParametro_45[0], sizeof(nombreParametro_45));
  memcpy(&OutNam[45][0], &nombreParametro_46[0], sizeof(nombreParametro_46));
}

static void ConvolucionarSeries(int32_T n1, const real_T vector1[],
  int32_T n2, const real_T vector2[], real_T vectorConvolucion[])
{
  int32_T longitudVectorConvolucion;
  int32_T k;
  int32_T lim_b;
  real_T vC_k;
  int32_T j;
  longitudVectorConvolucion = (n1 + n2) - 1;

  for (k=1;k<=longitudVectorConvolucion;k++)
  {
    lim_b = ValorMenor(k, n1);
    vC_k = 0.0;
    for (j=ValorMayorEntero(1, (k - n2) + 1);j<=lim_b;j++)
      vC_k += vector1[j - 1] * vector2[k - j];

    vectorConvolucion[k - 1] = vC_k;
  }
}

static void CopiarFragmentoVectorEntero(int32_T inicioIntervalo, int32_T finIntervalo,
  const int32_T vectorOrigen[], int32_T vectorDestino[])
{
	memcpy(&vectorDestino[0], &vectorOrigen[inicioIntervalo-1], sizeof(int32_T)*(finIntervalo-inicioIntervalo+1));
}

static void CopiarFragmentoVector(int32_T inicioIntervalo, int32_T finIntervalo,
  const real_T vectorOrigen[], real_T vectorDestino[])
{
	memcpy(&vectorDestino[0], &vectorOrigen[inicioIntervalo-1], sizeof(real_T)*(finIntervalo-inicioIntervalo+1));
}


static void CopiarFragmentoVectorSinMediaConversionDesdeEntero(int32_T
  inicioIntervalo, int32_T finIntervalo, int32_T longitudVectorOrigen, const
  int16_T vectorOrigen[], real_T vectorDestino[])
{
  real_T media;
  int32_T longitud;
  int32_T inicio;
  int32_T i;

  media = ObtenerMediaSerieEnteraCorta(longitudVectorOrigen, vectorOrigen);
  longitud = (finIntervalo - inicioIntervalo) + 1;
  inicio= inicioIntervalo-1;
  for(i=0;i<longitud;i++)
    vectorDestino[i] = (real_T)vectorOrigen[inicio+i]- media;
}


static void CopiarVector(int32_T longitudVector, const real_T vectorOrigen[], real_T vectorDestino[])
{
	memcpy(&vectorDestino[0], &vectorOrigen[0], sizeof(real_T)*longitudVector);
}

static void CorregirSignoSerie(int32_T numeroElementos, const real_T serie[], real_T serieCorregida[])
{
  real_T maximo;
  real_T x;
  int32_T i;
  
  ObtenerMinimoMaximoSerie(numeroElementos, serie, &x, &maximo);
  if (Signo(fabs(x) - fabs(maximo)) == 1)
  {
    for(i=0;i<numeroElementos;i++)
      serieCorregida[i]= -serie[i];
  }
  else
  {
    for(i=0;i<numeroElementos;i++)
      serieCorregida[i]= serie[i];
  }
}

static int32_T DivisionEntera(int32_T dividendo, int32_T divisor)
{
  return ParteEntera((real_T)dividendo / (real_T)divisor);
}

static void EliminarMediaSerie(int32_T numeroElementos, const real_T serie[], real_T serieSinMedia[])
{
  real_T media;
  int32_T i;
  
  media = ObtenerMediaSerie(numeroElementos, serie);
  for(i=0;i< numeroElementos; i++)
    serieSinMedia[i] = serie[i] - media;
}

static creal_T ExponencialComplejaSoloParteImaginaria(real_T z_imag)
{
  creal_T exp_z;

  exp_z.re = cos(z_imag);
  exp_z.im = sin(z_imag);

  return exp_z;
}

static void FFT_radix2(int32_T longitudTrama, const creal_T trama[], creal_T salidaFFT[])
{
  FastFourierTransform(longitudTrama, trama, salidaFFT);
}

static void FFT_radix_no2(int32_T longitudTrama, const real_T trama[], int32_T numeroPuntosEnFrecuencia, creal_T salidaFFT[])
{
  static creal_T a[NUMERO_MAXIMO_PUNTOS_FFT];
  static creal_T b[NUMERO_MAXIMO_PUNTOS_FFT];
  static creal_T ch_p[NUMERO_MAXIMO_PUNTOS_FFT];
  static creal_T x_a[NUMERO_MAXIMO_PUNTOS_FFT];
  static creal_T x_ch[NUMERO_MAXIMO_PUNTOS_FFT];
  static creal_T phi[NUMERO_MAXIMO_PUNTOS_FFT];
  int32_T lp;
  int32_T P;
  real_T inv_M;
  real_T inv_M2;
  int32_T i;
  real_T ch_p_im;

  lp = longitudTrama + numeroPuntosEnFrecuencia;
  P = ObtenerSiguientePotencia2(lp - 2);
  inv_M = 1.0 / (real_T)numeroPuntosEnFrecuencia;
  inv_M2 = 0.5 * inv_M;
  for (i = 1; i <= longitudTrama; i++)
    a[i - 1] = ExponencialComplejaSoloParteImaginaria(-6.2831853071795862*
      (inv_M2 * (real_T)(i * i) - inv_M * (real_T)i));

  for (i = 1; i <= numeroPuntosEnFrecuencia; i++)
    b[i - 1] = ExponencialComplejaSoloParteImaginaria(-6.2831853071795862 *
      (inv_M2 * (real_T)(i * i) - inv_M * ((real_T)i - 1.0)));

  for (i = 1; i < lp; i++)
    ch_p[i - 1] = ExponencialComplejaSoloParteImaginaria(6.2831853071795862 *
      (inv_M2 * (real_T)((i - longitudTrama) * (i - longitudTrama))));

  while (lp<= P)
  {
    ch_p[lp - 1].re = 0.0;
    ch_p[lp - 1].im = 0.0;
    lp++;
  }

  FFT_radix2(P, ch_p, phi);
  for(i=0;i< longitudTrama; i++)
  {
    x_a[i].re = trama[i] * a[i].re;
    x_a[i].im = trama[i] * a[i].im;
  }

  for (i = longitudTrama; i<P; i++)
  {
    x_a[i].re = 0.0;
    x_a[i].im = 0.0;
  }

  FFT_radix2(P, x_a, ch_p);
  for(i=0;i< P; i++)
  {
    x_ch[i].re = ch_p[i].re * phi[i].re - ch_p[i].im * phi[i].im;
    x_ch[i].im = -(ch_p[i].re * phi[i].im + ch_p[i].im * phi[i].re);
  }

  FFT_radix2(P, x_ch, ch_p);
  inv_M = 1.0 / (real_T)P;
  for(i=0;i< numeroPuntosEnFrecuencia; i++)
  {
    inv_M2 = inv_M * ch_p[(i + longitudTrama) - 1].re;
    ch_p_im = inv_M * -ch_p[(i + longitudTrama) - 1].im;
    salidaFFT[i].re = inv_M2 * b[i].re - ch_p_im * b[i].im;
    salidaFFT[i].im = inv_M2 * b[i].im + ch_p_im * b[i].re;
  }
}

static int32_T IndicePrimerValorMayorEntero(int32_T numeroElementos, const int32_T serie[], int32_T valorReferencia)
{
  int32_T indice;
  int32_T i;
  boolean_T exitg1;
  indice = 1;
  i = 1;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (i <= numeroElementos)) {
    if (serie[i - 1] > valorReferencia) {
      indice = i;
      exitg1 = TRUE;
    } else {
      i++;
    }
  }

  return indice;
}

static int32_T IndicePrimerValorMayor(int32_T numeroElementos, const real_T serie[], real_T valorReferencia)
{
  int32_T indice;
  int32_T i;
  boolean_T exitg1;
  
  indice = 1;
  i = 1;
  exitg1 = FALSE;
  while((exitg1 == FALSE) && (i <= numeroElementos))
  {
    if (serie[i - 1] > valorReferencia)
	{
      indice = i;
      exitg1 = TRUE;
    }
	else
      i++;
  }

  return indice;
}

static int32_T IndiceUltimoValorMenor(int32_T numeroElementos, const real_T serie[], real_T valorReferencia)
{
  int32_T indice;
  int32_T i;
  boolean_T exitg1;

  indice = numeroElementos;
  i = numeroElementos;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (i > 0))
  {
    if (serie[i - 1] < valorReferencia)
	{
      indice = i;
      exitg1 = TRUE;
    }
	else
      i--;
  }

  return indice;
}

static int32_T IndiceUltimoValorMenorEntero(int32_T numeroElementos, const int32_T serie[], int32_T valorReferencia)
{
  int32_T indice;
  int32_T i;
  boolean_T exitg1;

  indice = numeroElementos;
  i = numeroElementos;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (i > 0))
  {
    if (serie[i - 1] < valorReferencia)
	{
      indice = i;
      exitg1 = TRUE;
    }
	else
      i--;
  }

  return indice;
}

static real_T LogaritmoBase10(real_T numero)
{
  real_T logaritmo10;

  if (numero > 1.0E-10)
    logaritmo10 = log10(numero);
  else
    logaritmo10 = -7.0;

  return logaritmo10;
}

static void MaximoIntervaloSerie(int32_T inicioIntervalo, int32_T finIntervalo, const real_T serie[], int32_T *indice, real_T *maximo)
{
  int32_T i;
  *indice = inicioIntervalo;
  *maximo = serie[inicioIntervalo - 1];

  for (i=inicioIntervalo;i<finIntervalo;i++)
  {
    if (*maximo < serie[i])
	{
      *maximo = serie[i];
      *indice = i + 1;
    }
  }
}

static void MaximoSerie(int32_T numeroElementos, const real_T serie[], int32_T *indice, real_T *maximo)
{
  real_T max;
  int32_T ind;
  int32_T i;

  ind= 0;
  max= serie[0];
  for(i=1;i<numeroElementos;i++)
  {
    real_T elemento= serie[i];

	if(max<elemento)
	{
      max= elemento;
      ind= i;
    }
  }
  *maximo= max;
  *indice= ind+1;
}

static void MaximoSerieEntera(int32_T numeroElementos, const int32_T serie[], int32_T *indice, int32_T *maximo)
{
  int32_T max;
  int32_T ind;
  int32_T i;

  ind= 0;
  max= serie[0];
  for (i=1;i<numeroElementos;i++)
  {
	int32_T elemento= serie[i];
	
	if(max<elemento)
	{
      max= elemento;
      ind= i;
    }
  }
  *maximo= max;
  *indice= ind+1;
}

static void MaximoValorAbsolutoIntervaloSerie(int32_T inicioIntervalo, int32_T finIntervalo, const real_T serie[], int32_T *indice, real_T *maximo)
{
  int32_T i;
  real_T maximo_abs;
  
  *indice = inicioIntervalo;
  *maximo = fabs(serie[inicioIntervalo - 1]);
  for(i=inicioIntervalo+1;i<=finIntervalo;i++)
  {
    maximo_abs = fabs(serie[i - 1]);
    if (*maximo < maximo_abs)
	{
      *maximo = maximo_abs;
      *indice = i;
    }
  }
}

static void MinimoIntervaloSerie(int32_T inicioIntervalo, int32_T finIntervalo, const real_T serie[], int32_T *indice, real_T *minimo)
{
  int32_T i;
  
  *indice = inicioIntervalo;
  *minimo = serie[inicioIntervalo-1];
  for(i=inicioIntervalo;i<finIntervalo;i++)
  {
    if(*minimo > serie[i])
	{
      *minimo= serie[i];
      *indice= i+1;
    }
  }
}

static void MinimoSerie(int32_T numeroElementos, const real_T serie[],
  int32_T *indice, real_T *minimo)
{
  int32_T i;
  
  *indice = 1;
  *minimo = serie[0];
  for(i=2;i<=numeroElementos;i++)
  {
    if(*minimo>serie[i - 1])
	{
      *minimo= serie[i - 1];
      *indice= i;
    }
  }
}

static void MinimoSerieEmpezandoPorFinal(int32_T numeroElementos, const real_T serie[], int32_T *indice, real_T *minimo)
{
  int32_T i;
  
  *indice = numeroElementos;
  *minimo = serie[numeroElementos - 1];
  for(i=numeroElementos-1;i>0;i--)
  {
    if (*minimo>serie[i - 1])
	{
      *minimo = serie[i - 1];
      *indice = i;
    }
  }
}

static void ObtenerFraccionesSenoPi(int32_T N, real_T serie[])
{
  int32_T i;

  for (i = 1; i <= N; i++)
    serie[i - 1] = sin(3.1415926535897931 * (real_T)i / (real_T)N);
}

static void ObtenerFragmentoTrama(int32_T longitudFFT, int32_T inicioTrama, int32_T finTrama, const real_T trama[], real_T tramaFFT[])
{
  int32_T longitudTrama;
  int32_T i;

  longitudTrama = finTrama - inicioTrama;
  for(i=0;i<longitudTrama+1;i++)
    tramaFFT[i] = trama[(inicioTrama + i) - 1];

  for(i=longitudTrama+1;i<longitudFFT;i++)
    tramaFFT[i] = 0.0;
}

static real_T ObtenerMediaSerieEnteraCorta(int32_T numeroElementos, const int16_T serie[])
{
  real_T media;
  int32_T i;
  
  media = 0.0;
  for(i=0;i<numeroElementos;i++)
    media+= (real_T)serie[i];

  media/= (real_T)numeroElementos;
  
  return media;
}

static real_T ObtenerMediaSerieIntervalo(int32_T inicio, int32_T fin, const real_T serie[])
{
  real_T media;
  int32_T i;
  
  media= 0.0;
  for(i=inicio-1;i<fin;i++)
    media += serie[i];

  media/= (real_T)((fin - inicio) + 1);
  
  return media;
}

static real_T ObtenerMediaSerie(int32_T numeroElementos, const real_T serie[])
{
  real_T media;
  int32_T i;
  
  media= 0.0;
  for(i=0;i<numeroElementos;i++)
    media+= serie[i];

  media /= (real_T)numeroElementos;
  return media;
}

static real_T ObtenerMediaSerieEntera(int32_T numeroElementos, const int32_T serie[])
{
  real_T media;
  int32_T i;
  
  media = 0.0;
  for(i=0;i<numeroElementos;i++)
    media+= (real_T)serie[i];

  media /= (real_T)numeroElementos;
  
  return media;
}


static void ObtenerMediayDesviacionEstandar(int32_T numeroElementos, const real_T serie[], real_T *media, real_T *desviacion)
{
  int32_T i;
  real_T dif;
  
  *media = ObtenerMediaSerie(numeroElementos, serie);
  *desviacion = 0.0;
  for (i = 1; i <= numeroElementos; i++)
  {
    dif = serie[i - 1] - *media;
    *desviacion += dif * dif;
  }

  *desviacion = sqrt(*desviacion / (real_T)(numeroElementos - 1));
}

static void ObtenerMediayDesviacionEstandarEntera(int32_T numeroElementos, const int32_T serie[], real_T *media, real_T *desviacion)
{
   real_T med;
   real_T acum_desv;
   int32_T i;
  
  med= ObtenerMediaSerieEntera(numeroElementos, serie);
  acum_desv= 0.0;
  for(i=0;i<numeroElementos;i++)
  {
	real_T dif;
    
	dif = (real_T)serie[i]- med;
    acum_desv+= dif*dif;
  }

  *media= med;
  *desviacion = sqrt(acum_desv/ (real_T)(numeroElementos-1));
}

static void ObtenerMinimoMaximoSerie(int32_T numeroElementos, const real_T serie[], real_T *minimo, real_T *maximo)
{
  real_T max, min;
  int32_T i;
  
  min= serie[0];
  max= serie[0];
  for(i=1;i<numeroElementos;i++)
  { 
    real_T valor= serie[i];

	if(min>valor)
      min= valor;
	else
	{
      if(max<valor)
        max= valor;
    }
  }
  *maximo= max;
  *minimo= min;
}

static void ObtenerSeparacionesComienzoSerie(int32_T numeroElementos, const
  int32_T serie[], int32_T lista[])
{
  int32_T i;
  
  for (i=0;i<(numeroElementos-1);i++)
    lista[i]= serie[i+1]-serie[i];
}

static int32_T ObtenerSiguientePotencia2(int32_T numero)
{
  int32_T exponente;
  real_T potencia2;
  
  frexp(numero, &exponente);
  potencia2= pow(2.0, (real_T)exponente);
  return (int32_T)ParteEntera(potencia2);
}

static real_T ObtenerSumaValoresCuadradoIntervaloSerie(int32_T inicioIntervalo,int32_T finIntervalo, const real_T serie[])
{
  real_T suma;
  int32_T i;
  
  suma = 0.0;
  for(i=inicioIntervalo-1;i<finIntervalo;i++)
    suma+= serie[i]*serie[i];

  return suma;
}

static real_T ObtenerSumaValoresIntervaloSerie(int32_T inicioIntervalo, int32_T finIntervalo, const real_T serie[])
{
  real_T suma;
  int32_T i;
  
  suma = 0.0;
  for(i=inicioIntervalo-1;i<finIntervalo;i++)
    suma+= serie[i];

  return suma;
}

static void ObtenerTramaEnventanada(int32_T longitudFFT, const real_T ventana[], int32_T inicioTrama, int32_T finTrama,
	const real_T trama[], real_T tramaEnventanada[])
{
  static real_T tramaFFT[NUMERO_MAXIMO_MUESTRAS_VOZ];
  int32_T i;
  
  ObtenerFragmentoTrama(longitudFFT, inicioTrama, finTrama, trama, tramaFFT);
  for(i=0;i<longitudFFT;i++)
    tramaEnventanada[i]= ventana[i]*tramaFFT[i];
}

static int32_T ParteEntera(real_T numero)
{
  return (int32_T)floor(numero);
}

static int32_T Redondear(real_T numero)
{
	real_T y;
	
	if(numero>= 0.5)
		y= floor(numero+ 0.5);
	else
	{
		if(numero> -0.5)
			y= 0.0;
		else
			y= ceil(numero- 0.5);
	}

	return((int32_T)y);
}

static int8_T Signo(real_T numero)
{
  int8_T signoNumero;

  if(numero> 0.0)
    signoNumero= 1;
  else if (numero< 0.0)
    signoNumero= -1;
  else
    signoNumero= 0;

  return signoNumero;
}

static real_T ValorMayor(real_T numero1, real_T numero2)
{
  real_T mayor;
  
  if(numero1>numero2)
    mayor= numero1;
  else
    mayor= numero2;

  return mayor;
}

static int32_T ValorMayorEntero(int32_T numero1, int32_T numero2)
{
  int32_T mayor;
  
  if(numero1>numero2)
    mayor= numero1;
  else
	mayor= numero2;

  return mayor;
}

static int32_T ValorMenor(int32_T numero1, int32_T numero2)
{
  int32_T menor;
  
  if(numero1<numero2)
    menor = numero1;
  else
    menor = numero2;

  return menor;
}

static void ValoresAbsolutosSerie(int32_T numeroElementos, const real_T serie[], real_T serieCorregida[])
{
  int32_T i;

  for(i=0;i<numeroElementos;i++)
	serieCorregida[i] = fabs(serie[i]);
}

static void VentanaHamming(int32_T N, real_T hamming[])
{
  real_T inv_N1;
  int32_T i;
  
  inv_N1= 1.0 / (real_T)(N - 1);
  for (i=0;i<N;i++)
    hamming[i] = 0.54 - 0.46 * cos(6.2831853071795862 * (real_T)i * inv_N1);
}


static void cel2(int32_T long_s, const real_T s[], const real_T r[],
                   int32_T K, real_T u[], real_T c[])
{
  static real_T f[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T p[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T f_sig[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T p_sig[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T g[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T q[NUMERO_MAXIMO_MUESTRAS_VOZ];
  real_T ccor;
  real_T varf;
  real_T varg;
  int32_T i;

  CopiarVector(long_s, s, f_sig);
  CopiarVector(long_s, r, p_sig);
  CopiarVector(long_s, s, g);
  CopiarVector(long_s, r, q);
  for(i=0;i<K+1;i++)
  {
    int32_T j;
  
	CopiarVector(long_s, p_sig, p);
    CopiarVector(long_s, f_sig, f);
    ccor = 0.0;
    varf = f[long_s-1] * f[long_s-1];
    varg = g[long_s-1] * g[long_s-1];
    for(j=long_s-2;j>=0;j--)
	{
      ccor += f[j + 1] * g[j];
      varf += f[j] * f[j];
      varg += g[j] * g[j];
    }

    c[i] = -ccor / sqrt(varf * varg);
    f_sig[0] = f[0];
    for(j=1;j<long_s;j++)
      f_sig[j] = f[j] + c[i] * g[j - 1];

    for(j=long_s-1;j>0;j--)
	{
      g[j] = c[i] * f[j] + g[j - 1];
      p_sig[j] = p[j] + c[i] * q[j - 1];
      q[j] = c[i] * p[j] + q[j - 1];
    }

    g[0] = c[i] * f[0];
    p_sig[0] = p[0];
    q[0] = c[i] * p[0];
  }

  CopiarVector(long_s, q, u);
}

static void cel1(int32_T long_s, const real_T s[], int32_T K, real_T x[])
{
  static real_T f_sig[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T f[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T b[NUMERO_MAXIMO_MUESTRAS_VOZ];
  real_T ccor;
  real_T varf;
  real_T varb;
  int32_T i;

  CopiarVector(long_s, s, f_sig);
  CopiarVector(long_s, s, b);
  for(i=0;i<K;i++)
  {
	int32_T j;

    CopiarVector(long_s, f_sig, f);
    ccor = 0.0;
    varf = f[long_s - 1] * f[long_s - 1];
    varb = b[long_s - 1] * b[long_s - 1];
    for(j=long_s-2;j>=0;j--)
	{
      ccor += f[j + 1] * b[j];
      varf += f[j] * f[j];
      varb += b[j] * b[j];
    }

    ccor = -ccor / sqrt(varf * varb);
    f_sig[0] = f[0];
    for (j=long_s-1;j>0;j--)
	{
      f_sig[j] = f[j] + ccor * b[j - 1];
      b[j] = ccor * f[j] + b[j - 1];
    }

    b[0] = ccor * f[0];
  }

  CopiarVector(long_s, b, x);
}

static real_T cepstralpitch(int32_T long_InputVec, const real_T InputVector[], int32_T SamplingFrequency, int8_T SwitchSignal)
{
  static real_T InterVector_l[NUMERO_MAXIMO_MUESTRAS_VOZ];
  real_T LogSpectralDens[NUMERO_MAXIMO_PUNTOS_FFT];
  real_T LogSpectralDensHigh[NUMERO_MAXIMO_PUNTOS_FFT];
  real_T CepstrumHigh[NUMERO_MAXIMO_PUNTOS_FFT];
  real_T vacio;
  int32_T longitudTramaFFT;
  int32_T TopLimit;
  int32_T tam;
  int32_T MinSize;
  int32_T argMaxCepstrumHigh;

  lowpass(long_InputVec, InputVector, DivisionEntera(long_InputVec, NUMERO_MAXIMO_TRAMAS_ANALISIS), InterVector_l);
  longitudTramaFFT = ValorMenor(long_InputVec, 8192);
  CalcularModulosLogaritmicosCompletosFFT(longitudTramaFFT, &InterVector_l[0], LogSpectralDens);
  TopLimit = DivisionEntera(longitudTramaFFT, 2);
  if (SwitchSignal == 1)
  {
    tam = DivisionEntera(longitudTramaFFT, 20);
    MinSize = DivisionEntera(TopLimit, 40);
  }
  else
  {
    tam = DivisionEntera(longitudTramaFFT, 80);
    MinSize = DivisionEntera(TopLimit, 100);
  }

  highpass(longitudTramaFFT, LogSpectralDens, tam, LogSpectralDensHigh);
  CalcularModulosFFT(longitudTramaFFT, LogSpectralDensHigh, CepstrumHigh);
  MaximoIntervaloSerie(MinSize, TopLimit, CepstrumHigh, &argMaxCepstrumHigh, &vacio);
  
  return (real_T)SamplingFrequency / (real_T)(argMaxCepstrumHigh-1);
}


static void det_sing_mw4(int32_T long_mwpsd, const real_T mwpsd[], real_T
  *val_end, int32_T *pos_end, real_T mx_1[], int32_T pos_mx_1[], real_T
  mn_1[], int32_T pos_mn_1[], real_T mx_2[], int32_T pos_mx_2[],
  real_T nsf[], char_T *rel_ban)
{
  real_T val_max_sal[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  real_T val_min_sal[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int32_T arg_max_sal[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int32_T dim_pro;
  int32_T arg_min_sal[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int32_T dim_min;
  int32_T dim_pat;

  maxmin1d(long_mwpsd, mwpsd, 1, &dim_pro, val_max_sal, arg_max_sal);
  maxmin1d(long_mwpsd, mwpsd, 0, &dim_min, val_min_sal, arg_min_sal);
  if (dim_min <= 1)
  {
    dim_pat = 1;
    *val_end = 1.0;
    *pos_end = 1;
    mx_1[0] = 1.0;
    pos_mx_1[0] = 1;
    mn_1[0] = 1.0;
    pos_mn_1[0] = 1;
    mx_2[0] = 1.0;
    pos_mx_2[0] = 1;
    nsf[0] = 1.0;
  }
  else
  {
    dim_pat = 0;
    if (dim_min < dim_pro)
      dim_pro = dim_min;

    for (dim_min = 0; dim_min + 1 < dim_pro; dim_min++)
	{
      mx_1[dim_pat] = val_max_sal[dim_min];
      pos_mx_1[dim_pat] = arg_max_sal[dim_min];
      mn_1[dim_pat] = val_min_sal[dim_min];
      pos_mn_1[dim_pat] = arg_min_sal[dim_min];
      mx_2[dim_pat] = val_max_sal[dim_min + 1];
      pos_mx_2[dim_pat] = arg_max_sal[dim_min + 1];
      nsf[dim_pat] = fabs(((val_max_sal[dim_min] + val_max_sal[dim_min + 1]) /
                           2.0 - val_min_sal[dim_min]) / (2.0 * (real_T)
        (arg_max_sal[dim_min + 1] - arg_max_sal[dim_min]) / (real_T)
        (arg_max_sal[dim_min + 1] + arg_max_sal[dim_min])));
      dim_pat++;
    }

    *val_end = mwpsd[long_mwpsd - 1];
    *pos_end = long_mwpsd;
  }

  if(dim_pat < 2)
  {
    *rel_ban = '\x00';
    *val_end = 1.0;
    *pos_end = 1;
    mx_1[0] = 1.0;
    pos_mx_1[0] = 1;
    mn_1[0] = 1.0;
    pos_mn_1[0] = 1;
    mx_2[0] = 1.0;
    pos_mx_2[0] = 1;
    nsf[0] = 1.0;
  }
  else
    *rel_ban = '\x01';
}

static void glosotime(int32_T long_GloSource, const real_T GloSource[],
                      real_T MeanGloSource, int32_T sf, real_T *tR1r, real_T
                      *tR2r, real_T *tO1r, real_T *tO2r, real_T *tMr, real_T
                      *AR1r, real_T *AR2r, real_T *AO1r, real_T *AO2r, real_T
                      *RT, real_T *OT, real_T *CT, real_T GloFlow[])
{
  real_T MinMaxResid[NUMERO_MAXIMO_PUNTOS_FFT];
  real_T TotalCycleTime;
  real_T inv_dif;
  real_T MaxGloSource;
  real_T Recov1Val;
  real_T Open2Val;
  real_T AmpRange;
  int32_T NoZeroXings;
  int32_T GloSoZeroXings[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int32_T MaxGloSourcePoint;
  int32_T Recov1Point;
  int32_T long_Recov1MaxResid;
  int32_T long_MaxEndResid;
  int32_T Open2Point;
  int32_T Open1Point;
  int32_T Recov2Point;
  int32_T i;

  TotalCycleTime = (real_T)(long_GloSource - 1) / (real_T)sf;
  GloFlow[0] = 0.0;
  for (i=1;i<long_GloSource;i++)
    GloFlow[i] = GloFlow[i - 1] + 0.99 * GloSource[i];

  ObtenerMinimoMaximoSerie(long_GloSource, GloFlow, &Open2Val, &inv_dif);
  inv_dif = 1.0 / (inv_dif - Open2Val);
  for(i=0;i<long_GloSource;i++)
    GloFlow[i] = (GloFlow[i] - Open2Val) * inv_dif;

  zeroxing(long_GloSource, GloSource, GloSoZeroXings, &NoZeroXings);
  MaximoSerie(long_GloSource, GloSource, &MaxGloSourcePoint, &MaxGloSource);
  AmpRange = MaxGloSource - GloSource[0];
  for(i=0;i<MaxGloSourcePoint;i++)
  {
    MinMaxResid[i] = GloSource[i] - ((real_T)i * (MaxGloSource- GloSource[0])/
		((real_T)MaxGloSourcePoint-1.0) + GloSource[0]);
  }

  if (ObtenerMediaSerie(MaxGloSourcePoint, MinMaxResid) > 0.0)
  {
    MaximoSerie(MaxGloSourcePoint, MinMaxResid, &Recov1Point, &inv_dif);
    long_Recov1MaxResid = (MaxGloSourcePoint - Recov1Point) + 1;
    Recov1Val = GloSource[Recov1Point - 1];
    for (i=0;i<long_Recov1MaxResid;i++)
	{
      MinMaxResid[i] = GloSource[(i + Recov1Point) - 1] - (((real_T)i)
        * (MaxGloSource - GloSource[Recov1Point - 1]) / (real_T)
        (MaxGloSourcePoint - Recov1Point) + GloSource[Recov1Point-1]);
    }

    MinimoSerieEmpezandoPorFinal(long_Recov1MaxResid, MinMaxResid,
      &long_MaxEndResid, &inv_dif);
    Open2Point = (long_MaxEndResid + Recov1Point) - 2;
    Open2Val = GloSource[Open2Point];
  }
  else
  {
    MinimoSerieEmpezandoPorFinal(MaxGloSourcePoint, MinMaxResid,
      &long_MaxEndResid, &inv_dif);
    Open2Point = long_MaxEndResid - 1;
    Open2Val = GloSource[long_MaxEndResid - 1];
    for (i=0;i<long_MaxEndResid;i++)
	{
      MinMaxResid[i] = GloSource[i] - ((real_T)(i-1- long_MaxEndResid) *
        (GloSource[long_MaxEndResid - 1] - GloSource[0]) / ((real_T)
        long_MaxEndResid - 1.0) + GloSource[0]);
    }

    MaximoSerie(long_MaxEndResid, MinMaxResid, &Recov1Point, &inv_dif);
    Recov1Val = GloSource[Recov1Point - 1];
  }

  if (Recov1Point != Open2Point + 1)
  {
    long_Recov1MaxResid = (Open2Point - Recov1Point) + 2;
    for (i=0;i<long_Recov1MaxResid;i++)
	{
      MinMaxResid[i] = GloSource[(i + Recov1Point) - 1] - ((real_T)i * (Open2Val - Recov1Val) / (real_T)((Open2Point - Recov1Point) + 1) +
        Recov1Val);
    }

    MinimoSerie(long_Recov1MaxResid, MinMaxResid, &Open1Point, &inv_dif);
    Open1Point = (Open1Point + Recov1Point) - 2;
    MaximoSerie(long_Recov1MaxResid, MinMaxResid, &Recov2Point, &inv_dif);
    Recov2Point = (Recov2Point + Recov1Point) - 2;
  }
  else
  {
    Open1Point = Open2Point;
    Recov2Point = Open2Point;
  }

  *tO2r = (real_T)(Open2Point + 1) / (real_T)sf / TotalCycleTime;
  *AO2r = (GloSource[Open2Point] + MeanGloSource) / AmpRange;
  *tR1r = (real_T)Recov1Point / (real_T)sf / TotalCycleTime;
  *AR1r = (Recov1Val + MeanGloSource) / AmpRange;
  long_MaxEndResid = (long_GloSource - MaxGloSourcePoint) + 1;
  for(i=0;i<long_MaxEndResid;i++)
  {
    MinMaxResid[i] = GloSource[(i + MaxGloSourcePoint) - 1] - ((real_T)i * (GloSource[long_GloSource - 1] - MaxGloSource) / (real_T)
      (long_GloSource - MaxGloSourcePoint) + MaxGloSource);
  }

  MinimoSerie(long_MaxEndResid, MinMaxResid, &long_Recov1MaxResid, &inv_dif);
  *tR2r = (real_T)(Recov2Point + 1) / (real_T)sf / TotalCycleTime;
  *tO1r = (real_T)(Open1Point + 1) / (real_T)sf / TotalCycleTime;
  *tMr = (real_T)MaxGloSourcePoint / (real_T)sf / TotalCycleTime;
  *AR2r = (GloSource[Recov2Point] + MeanGloSource) / AmpRange;
  *AO1r = (GloSource[Open1Point] + MeanGloSource) / AmpRange;
  long_Recov1MaxResid = IndiceUltimoValorMenorEntero(NoZeroXings, GloSoZeroXings, MaxGloSourcePoint);
  
  *RT = (real_T)GloSoZeroXings[0] / (real_T)sf / TotalCycleTime;
  *OT = (real_T)GloSoZeroXings[long_Recov1MaxResid - 1] / (real_T)sf / TotalCycleTime;
  *CT = (real_T)GloSoZeroXings[IndicePrimerValorMayorEntero(NoZeroXings, GloSoZeroXings, MaxGloSourcePoint) - 1] / (real_T)sf / TotalCycleTime;
}

static void glottalclipping(int32_T long_GloSource, const real_T GloSource[], real_T SelectThreshold,
	real_T ExpectedPitch, int32_T *ResNoMinima, int32_T ResGloSourceMinArg[])
{
  static real_T RectGloSource[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T VerySlowGloSource[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T SlowGloSource[NUMERO_MAXIMO_MUESTRAS_VOZ];
  real_T GloSourceMinVal[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  real_T MaxSpikeMark;
  real_T AreaThreshold;
  real_T SpikeMark[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int32_T LeftArgument[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int32_T RightArgument[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int32_T GloSourceMinArg[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int32_T LeftFlankIndex;
  int32_T RightFlankIndex;
  int32_T MinimaIndex;
  int32_T LeftMinIndex;
  int32_T j;
  int8_T IntervalDet;

  for (LeftFlankIndex=0;LeftFlankIndex<long_GloSource;LeftFlankIndex++)
  {
    if (GloSource[LeftFlankIndex] < 0.0)
	{
      RectGloSource[LeftFlankIndex] = GloSource[LeftFlankIndex] +
        GloSource[LeftFlankIndex];
    }
	else
      RectGloSource[LeftFlankIndex] = 0.0;
  }

  lowpass(long_GloSource, RectGloSource, ParteEntera(0.8 * (real_T)long_GloSource / ExpectedPitch), VerySlowGloSource);
  lowpass(long_GloSource, RectGloSource, ParteEntera(0.4 * (real_T)long_GloSource / ExpectedPitch), SlowGloSource);
  for(LeftFlankIndex=0;LeftFlankIndex<long_GloSource;LeftFlankIndex++)
  {
    if (SlowGloSource[LeftFlankIndex] > 0.0)
      SlowGloSource[LeftFlankIndex] = 0.0;

    if (VerySlowGloSource[LeftFlankIndex] > 0.0)
      VerySlowGloSource[LeftFlankIndex] = 0.0;
  }

  IntervalDet = 0;
  LeftFlankIndex = -1;
  RightFlankIndex = 0;
  MinimaIndex = 0;

  for(j=0;j<long_GloSource;j++)
  {
    if (SlowGloSource[j] < VerySlowGloSource[j])
	{
      if (IntervalDet == 0)
	  {
        IntervalDet = 1;
        LeftArgument[LeftFlankIndex+1] = j-1;
        LeftFlankIndex++;
      }
    }
	else
	{
      if (IntervalDet == 1)
	  {
        IntervalDet = 0;
        RightArgument[RightFlankIndex] = j-1;
        RightFlankIndex++;
        MinimoIntervaloSerie(LeftArgument[LeftFlankIndex], RightArgument[LeftFlankIndex], RectGloSource,
                             &GloSourceMinArg[MinimaIndex], &GloSourceMinVal[MinimaIndex]);
        MinimaIndex++;
      }
    }
  }

  ValoresAbsolutosSerie(MinimaIndex, GloSourceMinVal, SpikeMark);
  ObtenerMinimoMaximoSerie(MinimaIndex, SpikeMark, &AreaThreshold,
    &MaxSpikeMark);
  AreaThreshold += (MaxSpikeMark - AreaThreshold) * SelectThreshold / 2.0;
  LeftMinIndex= 0;
  for (LeftFlankIndex = 0; LeftFlankIndex<MinimaIndex; LeftFlankIndex++)
  {
    if (SpikeMark[LeftFlankIndex] >= AreaThreshold)
	{
      ResGloSourceMinArg[LeftMinIndex] = GloSourceMinArg[LeftFlankIndex];
      LeftMinIndex++;
    }
  }
	*ResNoMinima= LeftMinIndex;
}

static void highpass(int32_T long_Vi, const real_T Vi[], int32_T fils, real_T Vo[])
{
  int32_T i;

  for (i = 1; i <= fils; i++)
    Vo[i - 1] = Vi[i - 1] - ObtenerMediaSerieIntervalo(1, i, Vi);

  for (i=fils+1;i<=(long_Vi - fils);i++)
    Vo[i - 1] = Vi[i - 1] - ObtenerMediaSerieIntervalo(i - fils, i + fils, Vi);

  for(i=long_Vi-fils;i<long_Vi;i++)
    Vo[i] = Vi[i] - ObtenerMediaSerieIntervalo(i + 1, long_Vi, Vi);
}

static void highpassfast(int32_T long_inp_vec, const real_T inp_vec[], int32_T tamavpas, real_T out_vec[])
{
  real_T fnorm;
  real_T meansca;
  int32_T n;
 
  fnorm = 1.0 / (real_T)((tamavpas << 1) + 1);
  meansca = 0.0;
  for(n = 1; n <= long_inp_vec; n++)
  {
    if (n <= tamavpas + 1)
      meansca = 1.0 / (real_T)(n + tamavpas) * ObtenerSumaValoresIntervaloSerie(1, n + tamavpas, inp_vec);
	else if(n > long_inp_vec - tamavpas)
	{
      meansca = 1.0 / (real_T)(((long_inp_vec - n) + tamavpas) + 1) *
        ObtenerSumaValoresIntervaloSerie(n - tamavpas, tamavpas, inp_vec);
    }
	else
      meansca += fnorm * (inp_vec[(n + tamavpas) - 1] - inp_vec[(n - tamavpas) -2]);

    out_vec[n - 1] = inp_vec[n - 1] - meansca;
  }
}

static void integ(int32_T long_inp_vec, const real_T inp_vec[], real_T leakage, int32_T fl, real_T out_vec[])
{
  static real_T int_vec[NUMERO_MAXIMO_MUESTRAS_VOZ];
  int32_T i;

  int_vec[0] = inp_vec[0];
  for(i=1;i<long_inp_vec;i++)
    int_vec[i]= inp_vec[i] + leakage * int_vec[i - 1];

  if (fl!= 0)
    highpassfast(long_inp_vec, int_vec, fl, out_vec);
  else
    CopiarVector(long_inp_vec, int_vec, out_vec);
}

static void lowpass(int32_T long_Vi, const real_T Vi[], int32_T fils, real_T Vo[])
{
  int32_T i;

  for (i = 1; i <= fils; i++)
    Vo[i - 1] = ObtenerMediaSerieIntervalo(1, i, Vi);

  for (i = fils + 1; i <=(long_Vi - fils); i++)
    Vo[i - 1] = ObtenerMediaSerieIntervalo(i - fils, i + fils, Vi);

  for (i = (long_Vi - fils) + 1; i <= long_Vi; i++)
    Vo[i - 1] = ObtenerMediaSerieIntervalo(i, long_Vi, Vi);
}

static void maxmin1d(int32_T long_Fun, const real_T Fun[], int8_T mxmn,
                     int32_T *dim, real_T val_min_sal[], int32_T arg_min_sal[])
{
  int32_T j;

  *dim = 0;
  for(j=1;j+1<long_Fun;j++)
  {
    switch (mxmn)
	{
     case 1:
      if ((Fun[j] - Fun[j - 1] > 0.0) && (Fun[j] - Fun[j + 1] > 0.0) && (fabs
           ((Fun[j] - Fun[j - 1]) * (Fun[j + 1] - Fun[j])) > 0.0))
	  {
        arg_min_sal[*dim] = j + 1;
        val_min_sal[*dim] = Fun[j];
        (*dim)++;
      }
      break;

     case 0:
      if ((Fun[j] - Fun[j - 1] < 0.0) && (Fun[j] - Fun[j + 1] < 0.0) && (fabs
           ((Fun[j] - Fun[j - 1]) * (Fun[j + 1] - Fun[j])) > 0.0))
	  {
        arg_min_sal[*dim] = j + 1;
        val_min_sal[*dim] = Fun[j];
        (*dim)++;
      }
      break;
    }

    if (*dim + 1 == 1)
	{
      val_min_sal[0] = Fun[0];
      val_min_sal[long_Fun - 1] = Fun[long_Fun - 1];
      arg_min_sal[0] = 1;
      arg_min_sal[1] = long_Fun;
    }
  }
}


static void unbias(int32_T long_ugn, const real_T ugn[], int32_T num_min,
                   const int32_T arg_min[], real_T ugndbl[])
{
  real_T slop[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int32_T j;
  int32_T i;
  int32_T int_siz;
  real_T newval;

  for (j=0;j<arg_min[0];j++)
    ugndbl[j] = ugn[j] - ugn[arg_min[0] - 1];

  for(i=0;i+1<num_min;i++)
  {
    int_siz = arg_min[i + 1] - arg_min[i];
    slop[i] = (ugn[arg_min[i + 1] - 1] - ugn[arg_min[i] - 1]) / (real_T)int_siz;
    for (j = -1; j + 2 <= int_siz; j++)
	{
      newval = (ugn[arg_min[i] + j] - ugn[arg_min[i] - 1]) - slop[i] * (real_T)(j + 1);
      if (newval <= 0.0)
        newval = -newval;

      ugndbl[arg_min[i] + j] = newval;
    }
  }

  for(j=arg_min[num_min - 1];j<=long_ugn;j++)
    ugndbl[j - 1] = ugn[j - 1] - ugn[arg_min[num_min - 1] - 1];
}

static void zeroxing(int32_T long_InputVector, const real_T InputVector[],
                     int32_T ZeroXingPoints[], int32_T *NoZeroXings)
{
  int8_T signoAnterior;
  int32_T i;
  int8_T signoActual;
  
  *NoZeroXings = 0;
  signoAnterior = Signo(InputVector[0]);
  for(i=1;i<long_InputVector;i++)
  {
    signoActual = Signo(InputVector[i]);
    if (signoActual != signoAnterior)
	{
      (*NoZeroXings)++;
      ZeroXingPoints[*NoZeroXings - 1] = i;
    }

    signoAnterior = signoActual;
  }
}


void glottal(int32_T long_si, const int16_T si[], int32_T sf, int32_T
                  left_lim, int32_T righ_lim, int8_T sw22, int8_T *exec_status,
                  TDimensiones *Dimensiones, int8_T out_ban[], real_T
                  OutArr[][NUMERO_PARAMETROS_ANALISIS], char_T OutNam[][LONGITUD_CADENA_TEXTO], real_T par_vec_med[],
                  real_T par_vec_std[], real_T ResultsPointRef[], real_T
                  GloSourceRefPoint[], real_T GloFlowRefPoint[],
                  real_T ugn_lvl[], int32_T *ref_pnt, int32_T arg_min[],
                  real_T *dist_vfolds)
{
  static real_T se[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T s_conv[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T s_m[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T sl[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T sl_m[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T s_s[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T sv0[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T sv0_m[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T fg_m[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T fg[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T dugn[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T dugn_m[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T xa[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T ugn[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T rugn[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T ugn_r[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T mucwav[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T avacwv[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T vacio[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T gsc[NUMERO_MAXIMO_MUESTRAS_VOZ];
  static real_T h[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T QualityFactor[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_ugn_lvl[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T men_min[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T pitch[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T CyclePer[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T jitter_abs_nor[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T shimmer_area_nor[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T acu_min_nor[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T t_tMr[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T t_tO1r[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T t_tO2r[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T t_tR1r[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T t_tR2r[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_AO1r[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_AO2r[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_AR1r[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_AR2r[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_CT[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_OT[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_RT[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T p_max_r_2[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T p_max_r_4[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T p_min_r_1[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T p_min_r_2[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_end_r[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_max_r_1[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_max_r_2[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_max_r_4[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_min_r_1[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T v_min_r_2[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T saw[NUMERO_MAXIMO_PUNTOS_FFT];
  static real_T gap[NUMERO_MAXIMO_PUNTOS_FFT];
  static real_T gap_m[NUMERO_MAXIMO_PUNTOS_FFT];
  static real_T mucwavnor[NUMERO_MAXIMO_PUNTOS_FFT];
  static real_T selwinc[NUMERO_MAXIMO_PUNTOS_FFT];
  static real_T GlFlRf[NUMERO_MAXIMO_PUNTOS_FFT];
  static real_T Cpstr[NUMERO_MAXIMO_RASGOS_CEPSTRALES][NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T MFmw[NUMERO_MAXIMO_TRAMAS_ANALISIS][NUMERO_MAXIMO_PUNTOS_FFT];
  static real_T MFmwdB[NUMERO_MAXIMO_TRAMAS_ANALISIS][NUMERO_MAXIMO_PUNTOS_FFT];
  static real_T MFuglp[NUMERO_MAXIMO_TRAMAS_ANALISIS][NUMERO_MAXIMO_PUNTOS_FFT];
  static real_T NHR[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T rel_pot_1[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T IndexRoughDistanceToLips[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T nsf_1[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T nsf_2[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T p_end_r[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  static real_T PrOutA[NUMERO_MAXIMO_TRAMAS_ANALISIS][NUMERO_PARAMETROS_ANALISIS];
  real_T cref[NUMERO_MAXIMO_ETAPAS_CELOSIA];
  real_T min_1[2];
  real_T pos_min_1[2];
  real_T max_2[2];
  real_T pos_max_2[2];
  real_T media_gap;
  real_T CT, OT, RT, AO2r, AO1r, AR2r, AR1r, tMr, tR2r, tR1r, tO2r, tO1r;
  real_T inv_media_men_min;
  real_T min_dugK, max_dugK;
  real_T GloCepPitch, SigCepPitch;
  real_T fd;
  real_T inv_dif;
  real_T SizeTrialStep;
  real_T pitch_ave;
  real_T media_v_ugn_lvl;
  real_T min_gsc, max_gsc;
  real_T g0;
  real_T saw2;
  real_T separacion;
  real_T ExpectedPitch;
  int32_T int_per[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int32_T arg_min_salr[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int32_T CycleDuration[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int32_T p_max_r_1[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int32_T long_win;
  int32_T num_min;
  int32_T long_ugn;
  int32_T n_factores;
  int32_T OptimIteration;
  int32_T ActualNoExpCycles;
  int32_T NoTrialSteps;
  int32_T long_se;
  int32_T long_s;
  int32_T long_h;
  int32_T st_tra;
  int32_T nd_tra;
  int32_T p;
  int32_T K;
  int32_T tama;
  int32_T med_win;
  int32_T nd_pnt;
  int32_T st_pnt;
  int32_T noven;
  int32_T max_per;
  int32_T vacioEntero;
  int32_T i;
  int32_T MinNoExpCycles, MaxNoExpCycles;
  boolean_T exitg1;
  int8_T rel_flg[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int8_T rel_ban_b[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int8_T rel_ban_c[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  int8_T rel_ban_d[NUMERO_MAXIMO_TRAMAS_ANALISIS];
  char_T rel_ban;
  int32_T j;

  
  Dimensiones->n_Param = 0;
  Dimensiones->n_OutFlags = 0;
  Dimensiones->n_ResultsPointRef = 0;
  Dimensiones->n_GloSourcePoints = 0;
  Dimensiones->n_GlottalSourceLevd = 0;
  Dimensiones->n_MinArg = 0;
  memset(&out_ban[0], 0, NUMERO_MAXIMO_TRAMAS_ANALISIS * sizeof(int8_T));
  memset(&OutArr[0][0], 0, NUMERO_PARAMETROS_ANALISIS*NUMERO_MAXIMO_TRAMAS_ANALISIS * sizeof(real_T));
  memset(&ResultsPointRef[0], 0, NUMERO_MAXIMO_TRAMAS_ANALISIS * sizeof(real_T));
  for (i=0;i<NUMERO_PARAMETROS_ANALISIS;i++)
  {
    par_vec_med[i] = 0.0;
    par_vec_std[i] = 0.0;
  }
  for (i=0;i<NUMERO_MAXIMO_PUNTOS_FFT;i++)
  {
    GloSourceRefPoint[i] = 0.0;
    GloFlowRefPoint[i] = 0.0;
  }

  memset(&ugn_lvl[0], 0, NUMERO_MAXIMO_MUESTRAS_VOZ * sizeof(real_T));
  *ref_pnt = 0;
  memset(&arg_min[0], 0, NUMERO_MAXIMO_TRAMAS_ANALISIS * sizeof(int32_T));
  *dist_vfolds = 0.0;
  ConfeccionarListaParametrosSalida(OutNam);



  if ((real_T)long_si / (real_T)sf < 0.05)
  {
    st_tra = 1;
    nd_tra = long_si;
  }
  else
  {
    st_tra = left_lim;
    nd_tra = righ_lim;
  }
  long_se= nd_tra-st_tra+1;

  CopiarFragmentoVectorSinMediaConversionDesdeEntero(st_tra, nd_tra, long_si, si, se);
  SigCepPitch = (real_T)sf / 2.0;
  fd = 0.9 * SigCepPitch;
  p= Redondear(3.3 * SigCepPitch / (SigCepPitch - 0.8 * SigCepPitch));
  long_h= p+1;
  for (i=1;i<=long_h; i++) {
    h[i - 1] = (0.54 - 0.46 * cos(6.2831853071795862 * ((real_T)i - 1.0) /
      (real_T)p)) * sin(6.2831853071795862 * fd / (real_T)sf * (((real_T)i
      - 1.0) - (real_T)p / 2.0)) / (3.1415926535897931 * (((real_T)i - 1.0)
      - (real_T)p / 2.0));
  }

  ConvolucionarSeries(long_se, se, long_h, h, s_conv);
  long_s= long_se+long_h-1;
  EliminarMediaSerie(long_s, s_conv, s_m);
  CorregirSignoSerie(long_s, s_m, s_s);
  for (i=1;i<long_s; i++) {
    sl[i]= s_s[i]- 0.95*s_s[i-1];
  }
  sl[0]= sl[1];
  EliminarMediaSerie(long_s, sl, sl_m);

  K= 4;
  cel1(nd_tra, sl_m, K, sv0);
  EliminarMediaSerie(long_s, sv0, sv0_m);
  cel2(long_s, sv0_m, sl_m, K, fg, cref);
  EliminarMediaSerie(long_s, fg, fg_m);

  tama= long_se-1;
  ObtenerMinimoMaximoSerie(tama, fg_m, &min_dugK, &max_dugK);
  inv_dif= 1.0 / (max_dugK - min_dugK);
  for(i=0;i<tama; i++) {
    dugn[i] = fg_m[i]* inv_dif;
  }

  EliminarMediaSerie(tama, dugn, dugn_m);
  integ(tama, dugn_m, 0.99, DivisionEntera(sf, 400), xa);
  
  if (sw22== 1)
    CambiarSignoSerie(tama, xa, ugn);
  else
    CopiarVector(tama, xa, ugn);

  NoTrialSteps= 20;
  SizeTrialStep= 2.0/(real_T)NoTrialSteps;

  MinNoExpCycles= DivisionEntera(50 * tama, sf);
  MaxNoExpCycles= DivisionEntera(300 * tama, sf);
  EliminarMediaSerie(tama, ugn, rugn);
  SigCepPitch= cepstralpitch(long_s, sl_m, sf, 0);
  GloCepPitch= cepstralpitch(tama, rugn, sf, 1);
  ExpectedPitch = 0.5 * ValorMayor(SigCepPitch, GloCepPitch);
  
  n_factores= 0;
  exitg1 = FALSE;
  for(i=0;!exitg1 && i<NoTrialSteps;i++)
  {
    real_T LevelTrialStep;
	real_T MeanCycleDur, StandCycleDur;

	LevelTrialStep= SizeTrialStep*(real_T)i;
	glottalclipping(tama, rugn, LevelTrialStep, ExpectedPitch, &ActualNoExpCycles, arg_min_salr);
    if (ActualNoExpCycles < DivisionEntera(MaxNoExpCycles, 8))
      exitg1 = TRUE;
    else
	{
      real_T DistributionFactor;
	  int32_T NoCycles;

	  MaxNoExpCycles = ActualNoExpCycles;
      ObtenerSeparacionesComienzoSerie(ActualNoExpCycles, arg_min_salr, CycleDuration);
      ObtenerMediayDesviacionEstandarEntera(ActualNoExpCycles-1, CycleDuration, &MeanCycleDur, &StandCycleDur);
	  DistributionFactor=StandCycleDur/(MeanCycleDur*MeanCycleDur);
      NoCycles=ActualNoExpCycles-1;
      if(NoCycles<6)
            exitg1= TRUE;
      else
	  {
        QualityFactor[i]= DistributionFactor/(real_T)(NoCycles*NoCycles);
        n_factores++;
        if(ActualNoExpCycles<MinNoExpCycles)
          exitg1 = TRUE;
      }
    }
  }

  MinimoSerieEmpezandoPorFinal(n_factores, QualityFactor, &OptimIteration, vacio);
  glottalclipping(tama, rugn, (real_T)(OptimIteration-1)*SizeTrialStep, ExpectedPitch, &ActualNoExpCycles, arg_min_salr);
  if(ActualNoExpCycles< 6)
    *exec_status = 0;
  else
  {
    int32_T arg_min_0;
    int32_T limsupcep;
	int32_T noceps;

	CopiarFragmentoVectorEntero(2, ActualNoExpCycles-1, arg_min_salr, arg_min);
    for (i=0;i<ActualNoExpCycles-3; i++)
      pitch[i] = (real_T)sf/(real_T)(arg_min_salr[i+2]- arg_min_salr[i+1]);

    pitch_ave= ObtenerMediaSerie(ActualNoExpCycles-3, pitch);
    num_min= ActualNoExpCycles-2;
    Dimensiones->n_MinArg = num_min;
    for (i=0;i<num_min;i++)
      CyclePer[i] = (real_T)(arg_min[i+1]-arg_min[i]) / (real_T)sf;

    med_win= 5;
	nd_pnt= arg_min[num_min-1]+ med_win;
    arg_min_0= arg_min[0];
	if(arg_min_0> med_win)
	{
	  st_pnt= arg_min_0- med_win;
      CopiarFragmentoVector(st_pnt, nd_pnt, ugn, ugn_r);
      for(i=0;i<ActualNoExpCycles;i++)
        arg_min[i]+= (med_win+1-arg_min_0);
    }
	else
	{
      int32_T arg_min_1= arg_min[1]; 

	  st_pnt= arg_min_0- med_win;
      CopiarFragmentoVector(st_pnt, nd_pnt, ugn, ugn_r);
      for(i=0;i<ActualNoExpCycles;i++)
        arg_min[i]+= (med_win+1-arg_min_1);
    }

	long_ugn= nd_pnt-st_pnt+1;
    unbias(long_ugn, ugn_r, num_min, arg_min, ugn_lvl);
    Dimensiones->n_GlottalSourceLevd = long_ugn;

	ObtenerSeparacionesComienzoSerie(num_min, arg_min, int_per);
    for (i=0;i<ActualNoExpCycles-4;i++)
      jitter_abs_nor[i] = fabs((pitch[i+1]-pitch[i])/pitch_ave);

    for (i=0;i<num_min-1; i++)
      v_ugn_lvl[i]= ObtenerSumaValoresIntervaloSerie(arg_min[i], arg_min[i+1], ugn_lvl);

    media_v_ugn_lvl= ObtenerMediaSerie(num_min-1, v_ugn_lvl);
    for (i=0;i<num_min-2;i++)
      shimmer_area_nor[i]= fabs(v_ugn_lvl[i+1]-v_ugn_lvl[i])/media_v_ugn_lvl;

    for (i=0;i<num_min-2;i++)
      men_min[i]= ObtenerMediaSerieIntervalo(arg_min[i]-med_win, arg_min[i]+med_win, ugn_lvl);

    inv_media_men_min= 1.0/ObtenerMediaSerie(num_min-2, men_min);
    for (i=0;i<num_min-2;i++)
      acu_min_nor[i] = men_min[i]* inv_media_men_min;

    noven= num_min;
	*ref_pnt = Redondear(((real_T)noven)/ 2.0);
    MaximoSerieEntera(num_min-1, int_per, &vacioEntero, &max_per);
    long_win = 0;
    for(i=0;i<arg_min[0];i++)
    {
      mucwav[i] = 0.0;
      avacwv[i] = 0.0;
    }

    for(i=arg_min[noven-1]-1;i<long_ugn;i++)
	{
      mucwav[i] = 0.0;
      avacwv[i] = 0.0;
    }

    for(j=0;j<noven-1;j++)
	{
	  int32_T str_win, end_win;
	  int32_T x;

	  str_win= arg_min[j];
	  end_win= arg_min[j+1];
	  long_win= end_win-str_win+1;
      
	  CopiarFragmentoVector(str_win, end_win, ugn_lvl, gsc);
      ObtenerMinimoMaximoSerie(long_win, gsc, &min_gsc, &max_gsc);
      separacion= max_gsc-min_gsc;
      if (separacion!= 0.0)
	  {
        ObtenerFraccionesSenoPi(long_win, saw);
        inv_dif= 1.0/separacion;
        for(x=0;x<long_win;x++)
          gap[x]= gsc[x]* inv_dif;

        g0 = 0.0;
        saw2 = 0.0;
        for(x=0;x<long_win;x++)
		{
          g0+= gap[x]* saw[x];
          saw2+= saw[x]*saw[x];
        }

        g0/= saw2;
        for(x=0;x<long_win;x++)
		{
          real_T aaw;

		  aaw= g0*saw[x];
          mucwav[str_win-1+x]= gap[x]- aaw;
          avacwv[str_win-1+x]= aaw;
        }

        EliminarMediaSerie(long_win, gap, gap_m);
        media_gap= ObtenerMediaSerie(long_win, gap);
        if (j+1 == *ref_pnt)
		{
          glosotime(long_win, gap_m, media_gap, sf, &tR1r, &tR2r,
                    &tO1r, &tO2r, &tMr, &AR1r, &AR2r, &AO1r, &AO2r, &RT, &OT,
                    &CT, GlFlRf);
          ResultsPointRef[0] = tR1r;
          ResultsPointRef[1] = tR2r;
          ResultsPointRef[2] = tO1r;
          ResultsPointRef[3] = tO2r;
          ResultsPointRef[4] = tMr;
          ResultsPointRef[5] = AR1r;
          ResultsPointRef[6] = AR2r;
          ResultsPointRef[7] = AO1r;
          ResultsPointRef[8] = AO2r;
          ResultsPointRef[9] = RT;
          ResultsPointRef[10] = OT;
          ResultsPointRef[11] = CT;

	  Dimensiones->n_ResultsPointRef = 12;
          Dimensiones->n_GloSourcePoints = long_win;
          CopiarVector(long_win, gap, GloSourceRefPoint);
          CopiarVector(long_win, GlFlRf, GloFlowRefPoint);
        }
		else
		{
          glosotime(long_win, gap_m, media_gap, sf, &tR1r, &tR2r,
                    &tO1r, &tO2r, &tMr, &AR1r, &AR2r, &AO1r, &AO2r, &RT, &OT,
                    &CT, GlFlRf);
        }

        t_tR1r[j]= tR1r;
        t_tO1r[j]= tO1r;
        t_tR2r[j]= tR2r;
        t_tO2r[j]= tO2r;
        t_tMr[j]= tMr;
        v_AR1r[j]= AR1r;
        v_AR2r[j]= AR2r;
        v_AO1r[j]= AO1r;
        v_AO2r[j]= AO2r;
        v_RT[j]= RT;
        v_OT[j]= OT;
        v_CT[j]= CT;
      }
    }

    limsupcep= Redondear((real_T)max_per / 2.0);
	noceps= 14;
    for(j=0;j<noven-1;j++)
	{
	  static real_T avacwvnor[NUMERO_MAXIMO_PUNTOS_FFT];
	  static real_T MFugdbldB[NUMERO_MAXIMO_PUNTOS_FFT];
	  static real_T MFgscpr[NUMERO_MAXIMO_PUNTOS_FFT];
	  real_T valnor, inv_valnor;
	  real_T media_mucwavnor;
	  real_T factorMedia;
	  real_T nhr_num, nhr_den;
	  int32_T str_win, end_win, str_win_m1;
	  int32_T liminfcep;
	  int32_T min_elementos;
	  int32_T x;

	  str_win= arg_min[j];
	  end_win= arg_min[j+1]-1;
	  long_win= end_win-str_win+1;

	  VentanaHamming(int_per[j], selwinc);
      MaximoValorAbsolutoIntervaloSerie(str_win, end_win, avacwv, &vacioEntero, &valnor);
      
	  str_win_m1= str_win-1;
	  inv_valnor= 1.0/valnor;
	  for(x=0;x<long_win;x++)
		mucwavnor[x] = mucwav[str_win_m1+x]* inv_valnor;

	  media_mucwavnor= ObtenerMediaSerie(long_win, mucwavnor);
      factorMedia= 0.99 * media_mucwavnor;
      for(x=0;x<long_win;x++)
        mucwavnor[x]-= factorMedia;

      CalcularModulosLogaritmicosFFT(max_per, long_win, mucwavnor, &MFmw[j][0], &MFmwdB[j][0]);
      CalcularModulosLogaritmicosCompletosFFT_VentanaIntervalo(max_per, selwinc, str_win, end_win, ugn_lvl, MFugdbldB);
      inv_valnor= 1.0 / valnor;
      for(x=0;x<long_win;x++)
        avacwvnor[x] = avacwv[str_win_m1+x]* inv_valnor;

      CalcularModulosFFTConRelleno(max_per, long_win, avacwvnor, &MFuglp[j][0]);

      liminfcep= DivisionEntera(limsupcep, 25);
      if(liminfcep<1)
        liminfcep= 1;

      CalcularModulosFFT(max_per, MFugdbldB, MFgscpr);
	  nhr_num= ObtenerSumaValoresIntervaloSerie(liminfcep+1, limsupcep, MFgscpr);
	  nhr_den= ObtenerSumaValoresIntervaloSerie(1, limsupcep, MFgscpr);
      NHR[j] = nhr_num/nhr_den;

	  if(limsupcep<noceps)
        min_elementos= limsupcep;
	  else
        min_elementos= noceps;

      for(x=0;x<min_elementos;x++)
        Cpstr[x][j]= MFgscpr[x];
    }

    if ((*ref_pnt < 3) || (long_win-1 < *ref_pnt+2)) {
      *exec_status = 0;
    }
	else
	{
	  const real_T csound= 35400.0;
	  real_T S[NUMERO_MAXIMO_TRAMAS_ANALISIS];
	  int32_T inf_ind_lip, sup_ind_lip;
	  int32_T posvocfol;
	  int32_T long_MF;
	  int32_T size_arr;

	  for(j=0;j<num_min-1;j++)
	  {
		real_T rel_pot_1_num, rel_pot_1_den;
		int32_T str_win, end_win;
		
		str_win= arg_min[j];
		end_win= arg_min[j+1]-1;

	    rel_pot_1_num= ObtenerSumaValoresCuadradoIntervaloSerie(str_win, end_win, mucwav);
	    rel_pot_1_den= ObtenerSumaValoresCuadradoIntervaloSerie(str_win, end_win, avacwv);
		rel_pot_1[j]= rel_pot_1_num/rel_pot_1_den;
      }

      S[0]= 1.0;
      for(i=1;i<K;i++)
        S[i]= S[i-1]* (1.0- cref[i])/ (1.0 + cref[i]);

	  for (i=0;i<K;i++)
        IndexRoughDistanceToLips[i]= (real_T)i* csound/(real_T)sf;

      inf_ind_lip= IndicePrimerValorMayor(K, IndexRoughDistanceToLips, (real_T)5.0);
	  sup_ind_lip= IndiceUltimoValorMenor(K, IndexRoughDistanceToLips, (real_T)22.0);
	  MinimoIntervaloSerie(inf_ind_lip, sup_ind_lip, S, &posvocfol, vacio);
	  *dist_vfolds = IndexRoughDistanceToLips[posvocfol-1];

	  long_MF= DivisionEntera(max_per, 2);
      rel_ban= 0;
      for(j=0;j<noven-1;j++)
	  {
		real_T mx_1[NUMERO_MAXIMO_TRAMAS_ANALISIS], mn_1[NUMERO_MAXIMO_TRAMAS_ANALISIS];
		real_T mx_2[NUMERO_MAXIMO_TRAMAS_ANALISIS], nsf[NUMERO_MAXIMO_TRAMAS_ANALISIS];
		real_T val_end;
		int32_T pos_mx_1[NUMERO_MAXIMO_TRAMAS_ANALISIS], pos_mn_1[NUMERO_MAXIMO_TRAMAS_ANALISIS], pos_mx_2[NUMERO_MAXIMO_TRAMAS_ANALISIS];
		int32_T pos_end;

		det_sing_mw4(long_MF, &MFmwdB[j][0], &val_end, &pos_end, mx_1, pos_mx_1, mn_1, pos_mn_1, mx_2, pos_mx_2, nsf, &rel_ban);
        if(rel_ban== 1)
		{
		  real_T mx1;
          int32_T pmx1;

	      rel_flg[j] = 1;
		  mx1=mx_1[0];
          pmx1=pos_mx_1[0];

          for (i=0;i<2;i++)
		  {
            min_1[i]= mn_1[i]- mx1;
            pos_min_1[i]= (real_T)pos_mn_1[i] / (real_T)pmx1;
            max_2[i]= mx_2[i]- mx1;
            pos_max_2[i] = (real_T)pos_mx_2[i] / (real_T)pmx1;
          }

          p_max_r_1[j]= pos_mx_1[0];
          p_min_r_1[j] = pos_min_1[0];
          p_max_r_2[j] = pos_max_2[0];
          p_min_r_2[j] = pos_min_1[1];
          p_max_r_4[j] = pos_max_2[1];
          p_end_r[j] = (real_T)pos_end/ (real_T)pmx1;
          v_max_r_1[j] = mx_1[0];
          v_min_r_1[j] = min_1[0];
          v_max_r_2[j] = max_2[0];
          v_min_r_2[j] = min_1[1];
          v_max_r_4[j] = max_2[1];
          v_end_r[j] = val_end- mx1;
          nsf_1[j]= nsf[0];
          nsf_2[j]= nsf[1];
        }
		else
		{
          rel_flg[j] = 0;
          p_max_r_1[j] = 1;
          p_min_r_1[j] = 1.0;
          p_max_r_2[j] = 1.0;
          p_min_r_2[j] = 1.0;
          p_max_r_4[j] = 1.0;
          p_end_r[j] = 1.0;
          v_max_r_1[j] = 0.0;
          v_min_r_1[j] = 0.0;
          v_max_r_2[j] = 0.0;
          v_min_r_2[j] = 0.0;
          v_max_r_4[j] = 0.0;
          v_end_r[j] = 0.0;
          nsf_1[j] = 0.0;
          nsf_2[j] = 0.0;
        }
      }

      size_arr= noven- 2;
      for(i=0;i<size_arr-1;i++)
	  {
        int32_T im1, im2;

		im1= i+1;
        im2= i+2;
        PrOutA[0][i] = pitch[im1];
        PrOutA[1][i] = jitter_abs_nor[im1];
        PrOutA[2][i] = shimmer_area_nor[im1];
        PrOutA[3][i] = acu_min_nor[im1];
        PrOutA[4][i] = NHR[im1];
        PrOutA[5][i] = rel_pot_1[im1];
        PrOutA[6][i] = Cpstr[0][im1];
        PrOutA[7][i] = Cpstr[1][im1];
        PrOutA[8][i] = Cpstr[2][im1];
        PrOutA[9][i] = Cpstr[3][im1];
        PrOutA[10][i] = Cpstr[4][im1];
        PrOutA[11][i] = Cpstr[5][im1];
        PrOutA[12][i] = Cpstr[6][im1];
        PrOutA[13][i] = Cpstr[7][im1];
        PrOutA[14][i] = Cpstr[8][im1];
        PrOutA[15][i] = Cpstr[9][im1];
        PrOutA[16][i] = Cpstr[10][im1];
        PrOutA[17][i] = Cpstr[11][im1];
        PrOutA[18][i] = Cpstr[12][im1];
        PrOutA[19][i] = Cpstr[13][im1];
        PrOutA[20][i] = v_max_r_1[im1];
        PrOutA[21][i] = v_min_r_1[im1];
        PrOutA[22][i] = v_max_r_2[im1];
        PrOutA[23][i] = v_min_r_2[im1];
        PrOutA[24][i] = v_max_r_4[im1];
        PrOutA[25][i] = v_end_r[im1];
        PrOutA[26][i] = (real_T)p_max_r_1[im1];
        PrOutA[27][i] = p_min_r_1[im1];
        PrOutA[28][i] = p_max_r_2[im1];
        PrOutA[29][i] = p_min_r_2[im1];
        PrOutA[30][i] = p_max_r_4[im1];
        PrOutA[31][i] = p_end_r[im1];
        PrOutA[32][i] = nsf_1[im1];
        PrOutA[33][i] = nsf_2[im1];
        PrOutA[34][i] = t_tR1r[im1];
        PrOutA[35][i] = t_tR2r[im1];
        PrOutA[36][i] = t_tO1r[im1];
        PrOutA[37][i] = t_tO2r[im1];
        PrOutA[38][i] = t_tMr[im1];
        PrOutA[39][i] = v_AR1r[im1];
        PrOutA[40][i] = v_AR2r[im1];
        PrOutA[41][i] = v_AO1r[im1];
        PrOutA[42][i] = v_AO2r[im1];
        PrOutA[43][i] = v_RT[im1];
        PrOutA[44][i] = v_OT[im1];
        PrOutA[45][i] = v_CT[im1];
      }

	for(i=0;i<size_arr;i++)
	{
		rel_ban_d[i]=1;
		for (j=0;j<NUMERO_PARAMETROS_ANALISIS;j++)
        {
			real_T valor= PrOutA[j][i];
			if(_isnan(valor))
			{
				OutArr[i][j]= (real_T)0.0;
				rel_ban_d[i]=0;
			}
			else
				OutArr[i][j]= valor;
		}
	}


      *exec_status = 1;
      Dimensiones->n_Param = NUMERO_PARAMETROS_ANALISIS;
      Dimensiones->n_OutFlags = size_arr-1;
      for (i=0;i<(size_arr-1);i++) {
        if ((rel_flg[i+1] == 1) && (rel_ban_b[i+1] == 1) && (rel_ban_c[i+1] == 1) && (rel_ban_d[i] == 1)) {
          out_ban[i]= 1;
        }
		else
		{
          out_ban[i]= 0;
          *exec_status= 0;
        }
      }

      for (i=0;i<NUMERO_PARAMETROS_ANALISIS;i++)
	  {
		
		real_T listaOutArr[NUMERO_MAXIMO_TRAMAS_ANALISIS];
		int32_T long_listaOutArr;
		int32_T x;


        long_listaOutArr= size_arr-1;
	    for(x=0;x<long_listaOutArr;x++)
			listaOutArr[x]= OutArr[x][i];
        ObtenerMediayDesviacionEstandar(long_listaOutArr, listaOutArr, &par_vec_med[i], &par_vec_std[i]);
      }
    }
  }
}

