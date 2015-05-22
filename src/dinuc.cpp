/*
  2015 4 21

  H. Takahashi

  Novel Motif Finding tool, dinucleotide independencies

*/

#include <iostream>
#include <string>

#include "util.h"
#include "util_class.h"

using namespace std;

static void usage(char* cmd){
  printf("[USAGE]\n");
  printf(" %s [OPTIONS] fasta_file \n\n", cmd);
  printf("[OPTIONS]\n");
  printf(" -h: display this help\n");
  printf(" -a: DMIN [default:2]\n");
  printf(" -b: DMAX [default:10]\n");
  printf(" -l: length of left motif [default:12]\n");
  printf(" -r: length of right motig [default:12]\n");
  printf(" -i: number of iterations [default:100]\n");
  printf("\n");
  exit(0);
}


int main(int argc, char **argv){

  char* fileName; //fileName=argv[1];

  int nSeq, result, dmin=2, dmax=10, left=12, right=12, iter=100;

  while((result=getopt(argc, argv, "ha:b:l:r:i:"))!=-1){
    switch(result){
    case 'h':
      usage(argv[0]);
    case 'a':
      dmin=atoi(optarg);
      break;
    case 'b':
      dmax=atoi(optarg);
      break;
    case 'i':
      iter=atoi(optarg);
      break;
    case 'l':
      left=atoi(optarg);
      break;
    case 'r':
      right=atoi(optarg);
      break;
    case '?':
      exit(1);
    }
  }

  /* get fasta file*/
  if(argc==optind){
    usage(argv[0]);
  }
  fileName=argv[optind];

  /* check DMIN, DMAX*/
  if(dmin>dmax){
    printf("DMIN should be less than equal DMAX!!\n");
    usage(argv[0]);
  }

  /* print arguments*/
  printf("%s -a %d -b %d -l %d -r %d -i %d %s\n",
	 argv[0], dmin, dmax, left, right, iter, fileName);

  paramSet param;
  param.assign(dmin, dmax, left, right, iter);

  param.show();

  nSeq=countSeq(fileName);
  cout<<"# of sequences: "<<nSeq<<endl;

  seqSet *seq=new seqSet[nSeq];
  getSeqSet(fileName, nSeq, seq);

  minimizeH(seq, param, nSeq);

  return 0;
}
