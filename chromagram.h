#ifndef __chromagram__chromagram__
#define __chromagram__chromagram__

#define nhops 1

#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "ifgram.h"
#include "ifptrack.h"


#endif /* defined(__chromagram__chromagram__) */

//Converts a frequency in Hz to a real number counting the octaves above A0 = 27.5
float hz2octs(float freq, float A0);

//Convert a real-number octave above A0 = 27.5 into a frequency in Hz.
float octs2hz(float octs, float A0);

//Find the maximum element of an array
float find_max(float* array, short numElements);

//Find index of max of an array - integer array
short find_index_max(short* array, short numElements);

//Find the minimum element of an array
float find_min(float* array, short numElements);

//Figure out best tuning alignment
//[counts,centers] = hist(x) so hx is the frequencies and hn is the number of occurrences in that frequency bin

//Matlab hist function
//void hist(int* counts, float* centers, float* tuning, int tune_count, int num_bins);

//Full chromagram_IF function
void chromagram_IF(float *X, float C[][nhops], float p[][nhops], float m[][nhops], int nbin, int f_ctr, int f_sd);
