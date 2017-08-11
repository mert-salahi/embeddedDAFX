#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "ifgram.h"
#include "ifptrack.h"
#include "chromagram.h"

float A0 = 27.5;
float Pocts[maxbin][nhops];

short i, j, k, num_count;
short counts[num_bins];
float centers[num_bins];

float tuning[maxbin];
short tune_count;

short index_max;
float centsoff;
int Pmapc[maxbin][nhops];
int PoctsQ_int;
float PoctsQ_float;



//Converts a frequency in Hz to a real number counting the octaves above A0 = 27.5
float hz2octs(float freq, float A0){
    float octs = (log(freq/A0))/(log(2));
    return octs;
}

//Convert a real-number octave above A0 = 27.5 into a frequency in Hz.
float octs2hz(float octs, float A0) {
    float hz = A0*(pow(2, octs));
    return hz;
}

//Find the maximum element of a 1-D array
float find_max(float* array, short numElements) {
	int i;
    float max = -10000;
    for (i=0; i<numElements; i++) {
        if (array[i] > max) {
            max = array[i];
        }
    }
    return max;
}


//Find the index of the maximum element of a 1-D array
short find_index_max(short* array, short numElements) {
	int i;
    int max = -10000;
    int index_max = 0;
    for (i=0; i<numElements; i++) {
        if (array[i] > max) {
            max = array[i];
            index_max = i;
        }
    }
    return index_max;
}

//Find the minimum element of a 1-D array
float find_min(float* array, short numElements) {
	int i;
    float min = 10000;
    for (i=0; i<numElements; i++) {
        if (array[i] < min) {
            min = array[i];
        }
    }
    return min;
}

void hist(short* counts, float* centers, float* tuning, short tune_count) {
    int i, j;
    float max = find_max(tuning, tune_count);
    float min = find_min(tuning, tune_count);
    float bin_step = (max-min)/(num_bins);
    float start_center = min+(bin_step/2.0);

    for (i=0; i<num_bins; i++){
    	counts[i] = 0;
        centers[i] = start_center+(bin_step*i);
        float left_end = min+(bin_step*i);
        float right_end = min+(bin_step*(i+1));
        for (j=0; j<tune_count; j++) {
            if ((tuning[j] >= left_end) && (tuning[j] <= right_end)) {
                counts[i] = counts[i] + 1;
            }
        }
    }
}


//Full chromagram_IF function declaration
void chromagram_IF(float *X, float C[][nhops], float p[][nhops], float m[][nhops], int nbin, int f_ctr, int f_sd) {

    //perform the conversion into octaves above A0 = 27.5 Hz for ONLY the regions that matter.

    num_count = 0;
    tune_count = 0;

    for (i=0; i<maxbin; i++){
        for (j=0; j<nhops; j++) {
            //Discarding the frequency elements with a magnitude of 0
            if (p[i][j] < 0.0001) {
                Pocts[i][j] = 0;
            }
            else {
                //Compute the octaves of each non-zero element above A0=27.5Hz
                num_count = num_count + 1;
                Pocts[i][j] = hz2octs(p[i][j], A0);
                if (Pocts[i][j] > 0.0001) {
                    tuning[tune_count] = (nchr*Pocts[i][j]) - round(nchr*Pocts[i][j]);
                    tune_count = tune_count+1;
                }
            }
        }
    }



    //Figure out best tuning alignment
    //[counts,centers] = hist - centers is center of the bins and counts is the number of occurrences in that bin
    hist(counts, centers, tuning, tune_count);
    index_max = find_index_max(counts, num_bins);
    centsoff = centers[index_max];


    //Quantize to chroma bins
    for (i=0; i<maxbin; i++) {
        for (j=0; j<nhops; j++) {
            if (Pocts[i][j] > 0.0001) {
                PoctsQ_float = round(nchr*(Pocts[i][j]-(centsoff/nchr)))/nchr;
                PoctsQ_int = round(nchr*(Pocts[i][j]-(centsoff/nchr)))/nchr;
                Pmapc[i][j] = round(nchr*(PoctsQ_float-PoctsQ_int)); //Mapping octaves into chroma bins
                //Integer number of octaves means element falls into A bin
                if (Pmapc[i][j] == nchr) {
                    Pmapc[i][j] = 0;
                }
            }
            else {
                Pmapc[i][j] = -1;
            }
        }
    }


    //Calculate the chroma matrix!! Shows the strength of each pitch class at every windowed time by adding all of the magnitudes of the harmonics together!
    short y;
    for (k=0; k<nhops; k++) {
        for (i=0; i<(int)nchr; i++) {
            float sum = 0;
            for (j=0; j<maxbin; j++){
                if (i==Pmapc[j][k]) {
                    y=1;
                }
                else {
                    y=0;
                }
                sum = sum+(y*m[j][k]);
            }
            C[i][k] = sum;
        }
    }
}
