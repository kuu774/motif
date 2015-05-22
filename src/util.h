#ifndef _UTIL_
#define _UTIL_

#define MLEN 1000

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

#include "util_class.h"

using namespace std;

void chkFile(char file[]);
int countSeq(char file[]);
void getSeqSet(char file[], int nSeq, seqSet *seq);

void minimizeH(seqSet *seq, paramSet parm, int nSeq);

double initialH(seqSet *seq, paramSet parm, int nSeq,
		int DI[][MLEN], int MONO[][MLEN]);

void initDI(int di[][16], int len);

double calcMONOPWM(int mono[]);
double calcDIPWM(int di[][16], int len);

double calcEntropyMONO(seqSet *seq, paramSet parm, int nSeq,
		       int MONO[][MLEN]);
double calcEntropyDI(seqSet *seq, paramSet parm, int nSeq,
		     int DI[][MLEN]);

double calcDIST(seqSet *seq, int dmin, int dmax, int nSeq);
void showDIST(seqSet *seq, int dmin, int dmax, int nSeq);

void makeMONOTable(seqSet *seq, int MONO[][MLEN], int nSeq);
void makeDITable(seqSet *seq, int DI[][MLEN], int nSeq);

//void countMONO(string nuc, int mono[]);
//double calcMONOPWM(int mono[]);
//void countDI(string nuc, int di[][16]);

#endif
