#include "util_class.h"

void seqSet::showMOTIF(int left, int right){
  string L=seq.substr(sPos, left);
  string R=seq.substr(sPos+gap+left, right);

  printf("%3d", sPos+1);
  cout<<"  "<<L<<"  <";
  printf("%3d", gap);
  cout<<">  ";
  printf("%3d", sPos+left+gap+1);
  cout<<"  "<<R<<endl;
}

void paramSet::show(){
  cout<<"ITER:"<<ITER<<"   Left:"<<LM<<"   Right:"<<RM<<"   DMIN:"<<DMIN<<"   DMAX:"<<DMAX<<endl;
}

void paramSet::assign(int a_DMIN, int a_DMAX, int a_LM,
		      int a_RM, int a_ITER){
  DMIN=a_DMIN; DMAX=a_DMAX; LM=a_LM; RM=a_RM; ITER=a_ITER;
}

int paramSet::getLM(){return LM;}
int paramSet::getRM(){return RM;}
int paramSet::getDMIN(){return DMIN;}
int paramSet::getDMAX(){return DMAX;}
int paramSet::getITER(){return ITER;}

void resultSet::show(){
  cout<<"FLAG:"<<flg<<"  POS:"<<sPos<<"  GAP:"<<gap<<endl;
}

void resultSet::assign(int a_flg, int a_sPos, int a_gap){
  flg=a_flg; sPos=a_sPos; gap=a_gap;
}

int resultSet::getFLG(){return flg;}
int resultSet::getSPOS(){return sPos;}
int resultSet::getGAP(){return gap;}
