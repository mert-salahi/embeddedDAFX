#include "filters.h"
#include <math.h>
#include <stdio.h>
#include "vibratolu.h"
#include "pwrap.h"


//void vibrato(int n, int sw, short *signalDelay, int count)
short Vibrato(int n, int sw, short *signalDelay, int count,int j, int index, short input)
{
	#define width 1
	#define modFreq 4
	#define Delay 22
	#define L 66
	float Modfreq = 0.000090703;

	//length of entire delay
	//#define D 88

	//short signalDelay[D];
	//int j;
	//int count;//sample number
	//int sw;

//	short sample_data,Data;
//	float dMod,iFrac,frac,sinMod;
//	int  k,i,WIDTH,p;
//	WIDTH = Delay;

	short signal, sample_data,dNew,dOld,dCount,output,dOut, d3, d4, d5, d6, d23,debugVal,Data;
	float dMod,iFrac,frac,out_data,sinMod,Vali,Vali2;
	int  k,i,calc,WIDTH,p;
	//double i;
	WIDTH = Delay;
		sample_data = input; //input sample
		debugVal = sample_data;
		Data = sample_data;
		dCount = count;
		if(sw ==0){	//if rest of delay buffer filled with 0s
			if (count > 1){	//more than one value store in buffer
				if(count <L){
					for(j = count; j>=1; j--){
						signalDelay[j] = signalDelay[j-1]; //shift to the right
					}
				}
				else{
					for(j = count-1; j>=1; j--){
						signalDelay[j] = signalDelay[j-1];
					}
				}
			}
		}
		//once cycled through delay buffer a full Length L
		//start shifting from complete end of buffer
		else{
			for(j = L-1; j>=1; j--){
					signalDelay[j] = signalDelay[j-1];
				}
		}

		signalDelay[0] = (float)Data;
		sinMod = (float)sin(Modfreq*2*3.14159*n);
		//sinMod = MOD[index];
		iFrac = 1 + WIDTH + (float)WIDTH*sinMod;
		i = floor(iFrac);
		frac = iFrac - i;
		out_data = signalDelay[0]*0.5 + 0.5*signalDelay[i+1]*frac + signalDelay[i]*(1-frac)*0.5;
		//out_data = signalDelay[0]>>1 + (signalDelay[i+1]*frac)>>1 + (signalDelay[i]*(1-frac))*0.5;

		//out_data = signalDelay[i+1]*frac + signalDelay[i]*(1-frac);
		output = (short) out_data;
		//AIC23_data.channel[LEFT] = output;
		//AIC23_data.channel[RIGHT] =  output;
		//output_sample(AIC23_data.combo);
		return output;
}

short Tremolo(short ModTrem, short effRate,short input){
	//short depth = 1;
	float depth = 0.7;
	//short controller = 1;
	short sample_data;
	float out_data;
	//short ModTrem = 0;
	//short effRate = 4000;
	float offset = 1-depth;
	float m = (float)ModTrem*depth/effRate;
	sample_data = input_left_sample();
	out_data= (m + offset)*sample_data;
	//output_left_sample((short) out_data);
	return (short) out_data;

}

short *Delay1(short *buffer1, short dimension, short *p, short input){
	short output;
	short output_try;
	output = *pwrap(dimension,buffer1,p+dimension);
	//*p = input + 0.88*output;
	*p = input + output>>1;
	//p = pwrap;
	p = pwrap(dimension, buffer1, --p);
	//output_left_sample(*p);
	//output_left_sample(input);
	return p;
//	delayed = buffer[j];
//	output = input + delayed;
//	output_left_sample(output);
//	buffer[j] = input + delayed*GAIN;
//	if(++j >= BUF_SIZE) j = 0;
//	else (++j);
}
//
short *allPass1(short *buffer, short dimension, short **p, short input){
	short output;
	short s0;

	output = *pwrap(dimension,buffer,*p+dimension);

	s0 = input + output>>1;
	**p = s0;
	*p = pwrap(dimension, buffer, --*p);
//	output_left_sample(**p);

	return *p;
}

short Distortion(short x){

//Distortion Effect Variables
	#define a  0.67
	#define b  (1-a)
	float c = 1024;
	float y, xc = x/c;
	y = (1 - b*xc*xc)*x;
	if (x>c) y = a*c;
	if (x<-c)	y = -a*c;
	return ((short) y);
}
