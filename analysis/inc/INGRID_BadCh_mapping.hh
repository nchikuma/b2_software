#ifndef _INGRID_BadCh_mapping_H
#define _INGRID_BadCh_mapping_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<fstream>
#include<string.h>

#define MaxNum_BadCh 1000
#define CARD_FILE "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/ingrid_info/card.txt"


using namespace std;
static const int IDTPL = 0;
static const int IDVETO = 1;


class INGRID_BadCh_mapping{

  private:
    int BadCh_Map[MaxNum_BadCh][4];
    int Num_BadCh;

  public:
    INGRID_BadCh_mapping(){};
    ~INGRID_BadCh_mapping(){};
    bool badchannel(string filename,int mod, int view, int pln, int ch);
    bool badchannel(int *mod,int *plane,int *ch, bool *tpl, bool *veto);
    bool readbad(string cardfilename,int *number_of_bad, int bad[][4]);

    void set_BadCh(bool opt);
    bool is_BadCh(int mod,int view,int pln,int ch);

};



#endif

