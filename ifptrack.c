#include "ifptrack.h"
#include "ifgram.h"
#include <math.h>

int counter = 0;

//For fmin = 50 and fmax = 1050

//Bin indexes for comparing neighboring frequencies and finding plateaus in the instantaneous frequency matrix
int d1[maxbin] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,133};
int d2[maxbin] = {0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132};

//initialize starts and ends matrices here
int st[maxbin];
int en[maxbin];
int starts[maxbin];
int ends[maxbin];
int average_bin;

float frqs[maxbin];
float mags[maxbin];
float s_bump[maxbin];
float i_bump[maxbin];
float ddif[maxbin][nhops];
bool dgood[maxbin][nhops];

#pragma DATA_SECTION(ddif, ".EXT_RAM")
#pragma DATA_SECTION(dgood, ".EXT_RAM")
#pragma DATA_SECTION(frqs, ".EXT_RAM")
#pragma DATA_SECTION(mags, ".EXT_RAM")
#pragma DATA_SECTION(s_bump, ".EXT_RAM")
#pragma DATA_SECTION(i_bump, ".EXT_RAM")

float dot_product(float* u, float* v, int n)
{

	int i;
    float result = 0;
    for (i = 0; i < n; i++)
        result += v[i]*u[i];
    return result;
}

float max(float x, float y){
    float result;
    if (x>y) {
        result = x;
    }
    else {
        result = y;
    }
    return result;
}


//Full ifptrack function!
void ifptrack(float *X, float p[][nhops], float m[][nhops], float IF[][nhops], float DFT[][nhops], float fminl, float fminu, float fmaxl, float fmaxu) {
    int i, j , t;
//float fftlen = 2048; float w = 1024; float sr = 16000; float fminl = 50; float fmaxu = 1050;
    //maxbin is the index of the bin that corresponds to the frequency of 1050Hz

    //Find plateaus in ifgram - stretches where delta IF is < thr

    //Calculating the difference in magnitude between neighboring frequency bins
    for (i=0; i<maxbin; i++){
        //expected increment per bin = sr/N = Hz of separation between frequency bins
        //if this is true, that means, we can group these neighboring frequency bins together

        for (j=0; j<nhops; j++) {
            ddif[i][j] = IF[d1[i]][j]-IF[d2[i]][j];
            if (fabsf(ddif[i][j]) < binIncrement) {
                dgood[i][j] = 1;
            }
            else {
                dgood[i][j] = 0;
            }
            //at every time t, we are grouping frequencies together and weighting the magnitudes accordingly
    		//delete any single bins (both above and below are zero);
        }
    }

 	for(i=0; i<maxbin; i++){
 		for(j=0; j<nhops; j++){
 	        if ((dgood[d1[i]][j] > 0) || (dgood[d2[i]][j] > 0)) {
                dgood[i][j] = dgood[i][j];
            }
            else {
                dgood[i][j] = 0;
            }
 		}
 	}

    //Find non-zero regions for every time t - group together close plateus to account for noisy channel
    for (t=0; t<nhops; t++) {
        int num_starts = 0;
        int num_ends=0; //number of non-zero regions
        int bin = 0;
        if (dgood[bin][t] > 0) {
            st[bin] = 1;
            starts[num_starts] = 0;
            num_starts = num_starts+1;
        }
        else {
            st[bin] = 0;
        }

        //Compute start and end arrays for non-zero regions (regions with peak frequencies)
        for (bin=0; bin<(maxbin-1); bin++) {
            if((dgood[bin][t]==0) && (dgood[bin+1][t] > 0)) {
                st[bin+1] = 1;
                starts[num_starts] = bin+1;
                num_starts = num_starts+1;
            }
            else {
                st[bin+1] = 0;
            }

            if((dgood[bin][t]>0) && (dgood[bin+1][t] == 0)) {
                en[bin] = 1;
                ends[num_ends] = bin;
                num_ends = num_ends+1;
            }
            else {
                en[bin]=0;
            }
        }

        if ((dgood[maxbin-1][t] > 0) && (dgood[0][t] == 0)){
            en[maxbin-1] = 1;
            ends[num_ends] = maxbin-1;
            num_ends = num_ends+1;
        }
        else {
            en[maxbin-1] = 0;
        }


        //Finding the DFT and IF values at the local maximums
        for (i=0; i<num_starts; i++) {
            float s_bump[maxbin];
            float i_bump[maxbin];
            counter = 0;
            for (j=(starts[i]); j<(ends[i]+1); j++) {
                s_bump[counter] = DFT[j][t];
                i_bump[counter] = IF[j][t];
                counter = counter+1;
            }

            int zero_sum = 0;
            if (sumArray(s_bump, counter) < 0.00001) {
                zero_sum = 1;
            }

            //compute the magnitude of the peak frequencies
            float sum_bump = sumArray(s_bump, counter);
            frqs[i] = (dot_product(s_bump, i_bump, counter))/(sum_bump+zero_sum);
            mags[i] = sum_bump;

            //If the frequency is above the max, we are discarding it!
            if (frqs[i] > fmaxu) {
                frqs[i] = 0;
                mags[i] = 0;
            }
            //Scale the magnitude down if above fmaxl -> 1 octave fadeout
            else if (frqs[i] > fmaxl) {
                mags[i] = mags[i]*(max(0, (fmaxu-frqs[i])/(fmaxu-fmaxl)));
            }

            //Same procedure for frequencies in the lowest range of observation
            if (frqs[i] < fminl) {
                frqs[i] = 0;
                mags[i] = 0;
            }
            else if (frqs[i] < fminu) {
                mags[i] = mags[i]*((frqs[i]-fminl)/(fminu-fminl));
            }

            if (frqs[i] < 0) {
                frqs[i] = 0;
                mags[i] = 0;
            }
        }

        for (i=0; i<maxbin; i++) {
        	p[i][t] = 0;
        	m[i][t] = 0;
        }

        //average of start and end bin and store into pitch and magnitude arrays
        for (i=0; i<num_starts; i++) {
        	if ((starts[i]+ends[i]) == 1) {
        		average_bin = 0;
        	}
        	else {
            	average_bin = ((starts[i]+ends[i])+1)>>1;
        	}
            p[average_bin][t] = frqs[i];
            m[average_bin][t] = mags[i];
        }
    }

}
