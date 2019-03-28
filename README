<b2_software> := <full path to your b2_software (This directory)>

1. Compile INGRID library :
     $ cd <b2_software>/lib
     $ make
   This library path must be set in the following environments :
     - wagasci_software/Analysis/process/< wgSpillCorr, wgMerge_inglib, wgMerge_spill >/GNUmakefile
     - INGRID/ingrid_format/app/Makefile
     - INGRID/INGRID/v1r1/cmt/requirements
     - INGRID/INGRID/v1r1/app/DSTMaker.hxx
     - analysis/app/Makefile
     - b2mc/GNUmakefile

2. Compile WAGASCI software :
     $ cd <b2_software>/wagasci_software/Analysis
     $ make

3. Compile INGRID software :
- 3.1. Compile "ingrid_format"
     $ cd <b2_software>/INGRID/ingrid_format
     $ cd bsd; make
     $ cp lib/libBeamData.a <your library place>
     $ cp lib/so/libBeamData.so <your library place>
     $ cd ../app; make
- 3.2. Compile "INGRID"
     $ cd <b2_software>/INGRID/INGRID/v1r1/cmt
     $ cmt broadcast cmt config
   Check your software versions :
     - nd280Policy v2r21
     - testBase    v1r5
     - EXTERN      v3r1
     - MYSQL       v5r051a
     - ROOT        v5r24p00n02
     - oaEvent     v7r3
     - oaRawEvent  v3r5
   If the version is wrong, set the correct one :
     $ cmt create <package name> <version>
     $ cmt broadcast cmt config
     $ source setup.sh
     $ cmt broadcast make

4. Compile Analysis software
      $ cd <b2_software>/analysis/app
      $ make

5. Compile Monte Carlo Simulation
      $ cd <b2_software>/b2mc
      $ source setup.sh
      $ make
