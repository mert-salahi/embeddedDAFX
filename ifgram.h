#include "lookups.h"
#ifndef __chroma_features__ifgram__
#define __chroma_features__ifgram__

#include <stdio.h>


#define _USE_MATH_DEFINES

#endif

void yell(void);

//Sum of array
float sumArray(float *a, int numElements);

//For the fft function
void swap(float *v1, float *v2);
void fftshift(float *data, int count);

//FFT function
void four1(float data[], int nn, int isign);

void cfftr2_dit(float* x, float* w, short n);
//x	Pointer to Array of Dimension 2*N elements holding
//*			Input to and Outputs from function cfftr2_dit()
//*		w	Pointer to an array holding the coefficient (Dimension
//*			n/2 complex numbers)
//*		N	Number of complex points in x

void ifgram(float *X, float IF[][nhops], float DFT[][nhops]);
