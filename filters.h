#ifndef FILTERS_H_
#define FILTERS_H_

short Vibrato(int n, int sw, short *signalDelay, int count,int j,int index, short input);
short Tremolo(short ModTrem, short effRate,short input);
short *Delay1(short *buffer1,short dimension, short *p, short input);
short *allPass1(short *buffer, short dimension, short **p, short input);
short Distortion(short x);


#endif /*FILTERS_H_*/
