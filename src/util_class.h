#ifndef _UTIL_CLS_
#define _UTIL_CLS_

#include <string>
#include <iostream>
#include <stdio.h>

using namespace std;

class seqSet{
 public:
  string name;
  string seq;
  int sPos;
  int gap;
  int flag;
  void showMOTIF(int left, int right);
  ~seqSet(){};
};

class paramSet{
 private:
  int DMIN, DMAX, LM, RM, ITER;

 public:
  void show();
  void assign(int a_DMIN, int a_DMAX, int a_LM, int a_RM, int a_ITER);
  int getLM();
  int getRM();
  int getDMIN();
  int getDMAX();
  int getITER();
  ~paramSet(){};
};

class resultSet{
 private:
  int flg, sPos, gap;

 public:
  resultSet(){
    flg=0; sPos=0; gap=0;
  }
  void show();
  void assign(int a_flg, int a_sPos, int a_gap);
  int getFLG();
  int getSPOS();
  int getGAP();
  ~resultSet();
};

#endif
