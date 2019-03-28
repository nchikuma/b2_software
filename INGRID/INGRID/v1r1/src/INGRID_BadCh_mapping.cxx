#ifndef _INGRID_BadCh_mapping_C
#define _INGRID_BadCh_mapping_C

#include"INGRID_BadCh_mapping.hxx"
using namespace std;
static const int TPL = 0;
static const int VETO = 1;

bool INGRID_BadCh_mapping::badchannel(string filename,int mod, int view, int pln,int ch){

  int number_of_bad;
  int bad[1000][4];
  this->readbad(filename,&number_of_bad,bad);

  for(int i=0;i<number_of_bad;i++){
    //std::cout << bad[i][0] << " " << bad[i][1] << " " << bad[i][2] << " " << bad[i][3] << "\n";
    if(bad[i][0]==mod){
      if(bad[i][1]==view){
        if(bad[i][2]==pln){
          if(bad[i][3]==ch){
            return true;
          }//channel
        }//planen
      }//TPL or VETO
    }//mod
  }//i
  return false;  

}

bool INGRID_BadCh_mapping::badchannel(int *mod,int *plane,int *ch,bool *tpl,bool *veto){
  int number_of_bad=8;
  int bad[8][4]={
    //{module, TPLorVETO, Plane, Channel}
    //### channel 0 ~ 23 -> view = 1, channel 24 ~ 48 -> view = 0
    {0, VETO, 2, 12},
    {1, TPL, 1, 14 },
    {5, TPL, 2, 0  },
    {5, TPL, 4, 18 },
    {5, TPL, 4,  2 },
    {7, TPL, 5, 24 },
    {9, TPL, 0, 12 },
    {5, TPL, 7, 9  } 
  };

  for(int i=0;i<number_of_bad;i++){
    if(bad[i][0]==*mod){
      if((*tpl&&bad[i][1]==IDTPL)||(*veto&&bad[i][1]==IDVETO)){
        if(bad[i][2]==*plane){
          if(bad[i][3]==*ch){
            return true;
          }//channel
        }//planen
      }//TPL or VETO
    }//mod
  }//i
  return false;  
}


bool INGRID_BadCh_mapping::readbad(string cardfilename,int *number_of_bad, int bad[][4]){

  std::ifstream ifs(cardfilename.c_str());
  if(!ifs){
    std::cout << "There is no 'card.txt': "
     << cardfilename << endl;
    *number_of_bad = 0;
    return false;
  }
  char c_temp[1000];
  char c_bad[1000]="#bad";

  int badnum=0;
  int mod,view,pln,ch;

  while(ifs >> c_temp){
    if(strcmp(c_temp,c_bad)==0){
      ifs >> mod >> view >> pln >> ch;
      //std::cout << mod << " " << view << " " << pln << " " << ch << "\n";
      bad[badnum][0]=mod;
      bad[badnum][1]=view;
      bad[badnum][2]=pln;
      bad[badnum][3]=ch;
      badnum++;
    }
  }
  //std::cout << "badch num " << badnum << "\n";
  *number_of_bad = badnum;
  ifs.close();
  return true;
}


#endif
