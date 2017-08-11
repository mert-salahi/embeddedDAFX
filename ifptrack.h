#ifndef __chroma_features__ifptrack__
#define __chroma_features__ifptrack__
#define nhops 1

#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "ifgram.h"
#include "ifptrack.h"


#endif /* defined(__chroma_features__ifptrack__) */

float dot_product(float* u, float* v, int n);

float max(float x, float y);

void ifptrack(float *X, float p[][nhops], float m[][nhops], float IF[][nhops], float DFT[][nhops], float fminl, float fminu, float fmaxl, float fmaxu);
