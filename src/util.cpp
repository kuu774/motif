#include "util.h"

// check whether fasta exists
void chkFile(char file[]){
  ifstream myfile;

  myfile.open(file);
  if(myfile.fail()){
    cerr<<"ERROR: CANNOT open ["<<file<<"]!!"<<endl;
    exit(1);
  }
  myfile.close();
}


// get number of sequences
int countSeq(char file[]){
  int num=0;
  ifstream myfile;
  string buf;

  chkFile(file);

  myfile.open(file);
  while(myfile && getline(myfile, buf)){
    if(buf[0] == '>'){
      num++;
    }
  }
  myfile.close();
  return num;
  }

// get sequence set
void getSeqSet(char file[], int nSeq, seqSet *seq){
  int num=0;
  ifstream myfile;
  string buf;

  myfile.open(file);
  while(myfile && getline(myfile, buf)){
    if(buf[0] == '>'){
      buf.erase(0,1);
      seq[num].name=buf;
      num++;
    }
    else{
      seq[num-1].seq+=buf;
    }
  }
  myfile.close();
}

void minimizeH(seqSet *seq, paramSet parm, int nSeq){
  int i,j,k,num;
  int iter=parm.getITER();

  int left=parm.getLM(); int right=parm.getRM();
  int dmin=parm.getDMIN(); int dmax=parm.getDMAX();

  int motifLen, old_sPos, old_gap;

  seqSet *opt=new seqSet[nSeq];

  for(i=0; i<nSeq; i++){
    opt[i].name=seq[i].name; opt[i].seq=seq[i].seq;
  }

  string motif;
  double minimum=1e6;

  int DI[nSeq][MLEN]; int MONO[nSeq][MLEN];

  //
  // precount di, mono nucleotides

  makeDITable(seq, DI, nSeq);
  makeMONOTable(seq, MONO, nSeq);

  for(int a=0; a<iter; a++){
    double oldH=initialH(seq, parm, nSeq, DI, MONO);
    double newH, prevH;
    prevH=oldH;

    /*for(i=0; i<nSeq; i++){
      if(opt[i].flag==-1){
	continue;
      }
      opt[i].showMOTIF(left, right);
      }*/

    while(1){
      for(i=0; i<nSeq; i++){
	if(seq[i].flag==-1){
	  continue;
	}
	motif=seq[i].seq;
	motifLen=motif.length();

	for(j=dmin; j<=dmax; j++){
	  num=motifLen-left-right-j;
	  for(k=0; k<=num; k++){
	    old_sPos=seq[i].sPos; old_gap=seq[i].gap;
	    seq[i].sPos=k; seq[i].gap=j;
	    newH=calcEntropyDI(seq, parm, nSeq, DI)+calcDIST(seq, dmin, dmax, nSeq)+calcEntropyMONO(seq, parm, nSeq, MONO);

	    //cout<<i<<","<<j<<","<<k<<","<<newH<<","<<oldH<<endl;
	    if(newH < oldH){
	      oldH=newH;
	    }
	    else{
	      seq[i].sPos=old_sPos; seq[i].gap=old_gap;
	    }
	  }
	}
      }
      //cout<<prevH<<","<<oldH<<endl;
      if(oldH < prevH){
	prevH=oldH;
      }
      else if(prevH-oldH <= 1e-10){
	if(minimum > prevH){
	  minimum=prevH;
	  for(int m=0; m<nSeq; m++){
	    opt[m].sPos=seq[m].sPos; opt[m].flag=seq[m].flag;
	    opt[m].gap=seq[m].gap;
	  }
	}
	break;
      }
    }
    //cout<<a<<","<<prevH<<","<<minimum<<endl;
  }

  //cout<<minimum<<endl;

  double res=calcEntropyDI(opt, parm, nSeq, DI)+calcDIST(opt, dmin, dmax, nSeq)+calcEntropyMONO(opt, parm, nSeq, MONO);
  cout<<"minimum H="<<res<<endl;

  showDIST(opt, dmin, dmax, nSeq);

  for(i=0; i<nSeq; i++){
    if(opt[i].flag==-1){
      continue;
    }
    opt[i].showMOTIF(left, right);
  }
}

//
double initialH(seqSet *seq, paramSet parm, int nSeq, int DI[][MLEN],
		int MONO[][MLEN]){

  int i,len;
  int startL=-1;
  int left=parm.getLM(); int right=parm.getRM();
  int dmin=parm.getDMIN(); int dmax=parm.getDMAX();

  int motifLen=left+right+dmin;

  // initialize
  srand((unsigned) time(NULL));
  for(i=0; i<nSeq; i++){
    string eachSeq=seq[i].seq;
    len=eachSeq.length()-motifLen;

    if(len>0){
      startL=rand()%len;
      seq[i].sPos=startL;
      seq[i].gap=dmin;
      seq[i].flag=0;
    }
    else if(len==0){
      startL=0;
      seq[i].sPos=startL;
      seq[i].gap=dmin;
      seq[i].flag=0;
    }
    else{
      seq[i].sPos=-1;
      seq[i].gap=-1;
      seq[i].flag=-1;
    }
  }
  //return calcEntropyMONO(seq, parm, nSeq)+calcEntropyDI(seq, parm, nSeq);
  return calcEntropyDI(seq, parm, nSeq, DI)+calcDIST(seq, dmin, dmax, nSeq)+calcEntropyMONO(seq, parm, nSeq, MONO);
}

double calcEntropyMONO(seqSet *seq, paramSet parm, int nSeq, int MONO[][MLEN]){
  int i, startL, startR, count;

  int monoL[4]={0,0,0,0}; int monoR[4]={0,0,0,0};

  string nucL, nucR, motif;

  int left=parm.getLM();

  for(i=0; i<nSeq; i++){
    if(seq[i].flag==-1){
      continue;
    }

    startL=seq[i].sPos;
    startR=startL+seq[i].gap+left;

    count=MONO[i][startL]; monoL[count]++;
    count=MONO[i][startR]; monoR[count]++;
  }

  return calcMONOPWM(monoL)+calcMONOPWM(monoR);
}

double calcMONOPWM(int mono[]){
  int i,num=0;
  double p; double ent=0; double u=0.25;

  for(i=0; i<4; i++){
    num+=mono[i];
  }

  for(i=0; i<4; i++){
    p=((double) mono[i]+u)/(num+1);
    ent-=p*log(p)/log(2); //p*log2(p)
  }
  return ent;
}


double calcEntropyDI(seqSet *seq, paramSet parm, int nSeq, int DI[][MLEN]){
  int i,j, startL, endL, startR, endR, count;
  string nucL, nucR, motif;

  int left=parm.getLM(); int right=parm.getRM();

  int diL[left-1][16]; int diR[right-1][16];
  initDI(diL, left-1); initDI(diR, right-1);

  for(i=0; i<nSeq; i++){
    if(seq[i].flag==-1){
      continue;
    }
    motif=seq[i].seq;
    startL=seq[i].sPos;            endL=startL+left;
    startR=startL+seq[i].gap+left; endR=startR+right;

    for(j=startL+1; j<endL; j++){
      count=DI[i][j];
      diL[(j-startL-1)][count]++;
    }

    for(j=startR+1; j<endR; j++){
      count=DI[i][j];
      diR[(j-startR-1)][count]++;
    }
  }

  /*
  for(i=0; i<16; i++){
    for(j=0; j<right-1; j++){
      cout<<diR[j][i]<<",";
    }
    cout<<endl;
    }*/

  return calcDIPWM(diL, left-1)+calcDIPWM(diR, right-1);
}

double calcDIPWM(int di[][16], int len){
  int i,j,num=0;
  double p, ent=0, u=0.0625;
  
  for(i=0; i<16; i++){
    num+=di[0][i];
  }

  for(i=0; i<len; i++){
    for(j=0; j<16; j++){
      p=((double) di[i][j]+u)/(num+1);
      //diPWM[i][j]=ent;
      ent-=p*log(p)/log(2);
    }
  }
  return ent;
}

double calcDIST(seqSet *seq, int dmin, int dmax, int nSeq){
  int i, num=0;
  int d=dmax-dmin+1;
  int *freq=new int[d];
  double p, H=0.0;

  for(i=0; i<d; i++){
    freq[i]=0;
  }

  for(i=0; i<nSeq; i++){
    if(seq[i].flag==-1){
      continue;
    }
    int j=seq[i].gap-dmin;
    freq[j]++;
    num++;
  }

  for(i=0;i<d;i++){
    p=((double)freq[i]+1.0/d)/(num+1);
    H-=p*log(p)/log(2);
  }

  H+=1.0/log(2)*15.0/2.0/num;
  delete[] freq;
  return H;
}

void showDIST(seqSet *seq, int dmin, int dmax, int nSeq){
  int i, num=0;
  int d=dmax-dmin+1;
  int *freq=new int[d];

  for(i=0; i<d; i++){
    freq[i]=0;
  }

  for(i=0; i<nSeq; i++){
    if(seq[i].flag==-1){
      continue;
    }
    int j=seq[i].gap-dmin;
    freq[j]++;
    num++;
  }

  for(i=0; i<d; i++){
    cout<<(i+dmin)<<":"<<freq[i]<<endl;
  }


  delete[] freq;
}



void initDI(int di[][16], int len){
  int i,j;
  for(i=0; i<len; i++){
    for(j=0; j<16; j++){
      di[i][j]=0;
    }
  }
}

void makeMONOTable(seqSet *seq, int MONO[][MLEN], int nSeq){
  int i,j;
  string motif;

  for(i=0; i<nSeq; i++){
    motif=seq[i].seq;
    int length=motif.length();

    for(j=0; j<length; j++){
      string dinuc=motif.substr(j,1);
      if(dinuc=="A"){MONO[i][j]=0;}
      else if(dinuc=="C"){MONO[i][j]=1;}
      else if(dinuc=="G"){MONO[i][j]=2;}
      else if(dinuc=="T"){MONO[i][j]=3;}
    }
  }
}

void makeDITable(seqSet *seq, int DI[][MLEN], int nSeq){
  int i,j;
  string motif;

  for(i=0; i<nSeq; i++){
    motif=seq[i].seq;
    int length=motif.length();

    for(j=1; j<length; j++){
      string dinuc=motif.substr(j-1,2);
      if(dinuc=="AA"){DI[i][j]=0;}
      else if(dinuc=="AC"){DI[i][j]=1;}
      else if(dinuc=="AG"){DI[i][j]=2;}
      else if(dinuc=="AT"){DI[i][j]=3;}
      //
      else if(dinuc=="CA"){DI[i][j]=4;}
      else if(dinuc=="CC"){DI[i][j]=5;}
      else if(dinuc=="CG"){DI[i][j]=6;}
      else if(dinuc=="CT"){DI[i][j]=7;}
      //
      else if(dinuc=="GA"){DI[i][j]=8;}
      else if(dinuc=="GC"){DI[i][j]=9;}
      else if(dinuc=="GG"){DI[i][j]=10;}
      else if(dinuc=="GT"){DI[i][j]=11;}
      //
      else if(dinuc=="TA"){DI[i][j]=12;}
      else if(dinuc=="TC"){DI[i][j]=13;}
      else if(dinuc=="TG"){DI[i][j]=14;}
      else if(dinuc=="TT"){DI[i][j]=15;}
    }
  }
}



/*
void countMONO(string nuc, int mono[]){
  if(nuc=="A"){mono[0]++;}
  else if(nuc=="C"){mono[1]++;}
  else if(nuc=="G"){mono[2]++;}
  else if(nuc=="T"){mono[3]++;}
}

void countDI(string nuc, int di[][16]){
  int i;
  int len=nuc.length();

  for(i=1; i<len; i++){
    string dinuc=nuc.substr(i-1,2);
    if(dinuc=="AA"){di[i-1][0]++;}
    else if(dinuc=="AC"){di[i-1][1]++;}
    else if(dinuc=="AG"){di[i-1][2]++;}
    else if(dinuc=="AT"){di[i-1][3]++;}
    //
    else if(dinuc=="CA"){di[i-1][4]++;}
    else if(dinuc=="CC"){di[i-1][5]++;}
    else if(dinuc=="CG"){di[i-1][6]++;}
    else if(dinuc=="CT"){di[i-1][7]++;}
    //
    else if(dinuc=="GA"){di[i-1][8]++;}
    else if(dinuc=="GC"){di[i-1][9]++;}
    else if(dinuc=="GG"){di[i-1][10]++;}
    else if(dinuc=="GT"){di[i-1][11]++;}
    //
    else if(dinuc=="TA"){di[i-1][12]++;}
    else if(dinuc=="TC"){di[i-1][13]++;}
    else if(dinuc=="TG"){di[i-1][14]++;}
    else if(dinuc=="TT"){di[i-1][15]++;}
  }
}


double calcScore(string motif, pwm *matrix){

  char MONO[4]={'A', 'C', 'G', 'T'};
  string DI[16]={"AA","AC","AG","AT",
		 "CA","CC","CG","CT",
		 "GA","GC","GG","GT",
		 "TA","TC","TG","TT"};

  double pw=0.0;
  int i,j,k;

  for(k=0; k<4; k++){
    if(motif[0] == MONO[k]){
      pw+=log(matrix->first[k]);
      break;
    }
  }

  for(j=1; j<MOTIF; j++){
    string di=motif.substr(j-1, 2);
    for(k=0; k<16; k++){
      if(di == DI[k]){
	pw+=log(matrix->di[k][j]);
	break;
      }
    }
  }
  return pw;
}
*/

    /*
    res[i].assign(flg, startL, dmin);

    if(flg!=-1){
      // FOR LEFT MOTIF
      countMONO(eachSeq.substr(startL,1), monoL);
      countDI(eachSeq.substr(startL, left), diL);

      // FOR RIGHT MOTIF
      countMONO(eachSeq.substr(startR,1), monoR);
      countDI(eachSeq.substr(startR, right), diR);

      cout<<eachSeq.substr(startL, left)<<","<<eachSeq.substr(startR, right)<<endl;
    }
  }

  // inital entropy
  entropy=calcMONOPWM(monoL)+calcMONOPWM(monoR)+calcDIPWM(diL, left-1)+calcDIPWM(diR, right-1);

  // search minimum entropy by greedy algorithm
  while(iter>0){
    for(j=0; j<nSeq; j++){
      if(res[j].getFLG()==-1){
	continue;
      }

      int judge=0;
      string eachSeq=seq[j].seq;
      string monoL_old=eachSeq.substr(res[j].getSPOS(),1);
      string monoR_old=eachSeq.substr((res[j].getSPOS()+res[j].getGAP()+left),1);
      int end=eachSeq.length()-right;

      int rightS, dist=dmin;

      if(monoL_old=="A"){monoL[0]--;}
      else if(monoL_old=="C"){monoL[1]--;}
      else if(monoL_old=="G"){monoL[2]--;}
      else if(monoL_old=="T"){monoL[3]--;}

      if(monoR_old=="A"){monoR[0]--;}
      else if(monoR_old=="C"){monoR[1]--;}
      else if(monoR_old=="G"){monoR[2]--;}
      else if(monoR_old=="T"){monoR[3]--;}

      for(k=0; k<(eachSeq.length()-right-dmin); k++){
	for(l=dmin; l<=dmax; l++){
	  rightS=k+l;
	  if(rightS==end){
	    break;
	  }
	}
      }
    }
    iter--;
    }
*/
