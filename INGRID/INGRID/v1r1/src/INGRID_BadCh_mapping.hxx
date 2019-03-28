#ifndef _INGRID_BadCh_mapping_H
#define _INGRID_BadCh_mapping_H

#include<iostream>
#include<sstream>
#include<fstream>

#include "setup.hxx"

using namespace std;
static const int IDTPL = 0;
static const int IDVETO = 1;

class INGRID_BadCh_mapping{
private:

public:
  INGRID_BadCh_mapping(){};
  ~INGRID_BadCh_mapping(){};
  bool badchannel(string filename,int mod, int view, int pln, int ch);
  bool badchannel(int *mod,int *plane,int *ch, bool *tpl, bool *veto);
  bool readbad(string cardfilename,int *number_of_bad, int bad[][4]);

};
#endif

