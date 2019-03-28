#include "EventAction.hh"
#include "RunAction.hh"
#include "HLayerSD.hh"
#include "VLayerSD.hh"
#include "VetoSD.hh"

#include "DetectorDimension.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4VisManager.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Square.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"

#include <iostream>
#include <assert.h>
#include "Riostream.h"


EventAction::EventAction(RunAction* rac)
  :runaction(rac)
{
  hitsum    = new HitSummary();
  simhitsum = new SimHitSummary();
  Flag=0;
}
EventAction::EventAction(RunAction* rac,DetectorDimension* fdim)
  :runaction(rac)
{
  detdim = fdim;
  hitsum    = new HitSummary();
  simhitsum = new SimHitSummary();
  Flag=0;
}


EventAction::~EventAction()
{
  if(hitsum){
    delete hitsum;
  }
  if(simhitsum){
    delete simhitsum;
  }
}

void EventAction::BeginOfEventAction(const G4Event* anEvent)
{
#ifdef DEBUG_EVTACT
  G4cout << "==================================================="  << G4endl;
  G4cout << "========== EventAction::BeginOfEventAction ========"  << G4endl;
  G4cout << "           >>> Event " << anEvent->GetEventID() << G4endl;
  G4cout << "==================================================="  << G4endl;
#endif
}

void EventAction::EndOfEventAction(const G4Event* anEvent)
{
  
  if(Flag!=-1) {

    G4int event_id = anEvent->GetEventID();
    
    //Get number of stored trajectories
    G4TrajectoryContainer* trajectoryContainer = anEvent->GetTrajectoryContainer();
    G4int n_trajectories = -1;
    if(trajectoryContainer){
      n_trajectories = trajectoryContainer->entries();
    }
      
    //Periodic printing
    if(event_id>0&&event_id%1000 == 0){
      G4cout << ">>> Event " << anEvent->GetEventID() << G4endl;
      G4cout << "    " << n_trajectories 
	     << " trajectories stored in this event." << G4endl;
    }
    
    G4SDManager* SDManager= G4SDManager::GetSDMpointer();

    //Get Hit Collection of This Event
    G4HCofThisEvent* HCTE= anEvent-> GetHCofThisEvent();
    if(! HCTE){
      G4cout << ">> No hits in this events." << G4endl;
      return;
    }

    static G4int idcalx= -1;
    static G4int idcaly= -1;
    static G4int idcalv= -1;

    if(idcalx<0){
      idcalx= SDManager-> GetCollectionID("vlayerHitCollection");  
    }
    if(idcaly<0){
      idcaly= SDManager-> GetCollectionID("hlayerHitCollection");
    }
    if(idcalv<0){
      idcalv= SDManager-> GetCollectionID("vetoHitCollection");
    }
    
    VLayerHitsCollection* vlhitcol =  (VLayerHitsCollection*)HCTE-> GetHC(idcalx);
    HLayerHitsCollection* hlhitcol =  (HLayerHitsCollection*)HCTE-> GetHC(idcaly);
    VetoHitsCollection  * vehitcol =  (VetoHitsCollection*)  HCTE-> GetHC(idcalv);  

      
    int detid=0;
    int trkid=0;
    float edep=0.;
    double time=0.;
    float pos[3];
    int mod=-1, ch=-1, view=-1, pln=-1, pid=-1;
    double pe;
    double lope;

    int gridcell_id_x1=-1;
    int gridcell_id_x2=-1;
    int gridcell_id_y1=-1;
    int gridcell_id_y2=-1;

    int nsimhits    = 0;
    int nhitsvlayer = 0;
    int nhitshlayer = 0; 
    int nhitsveto   = 0; 

    if(vlhitcol){
      nhitsvlayer = vlhitcol->entries();
    }
    if(hlhitcol){
      nhitshlayer = hlhitcol->entries();
    }
    if(vehitcol){
      nhitsveto = vehitcol->entries();
    }

#ifdef SUPPRESS_FOR_EXT_BG  //To suppress entries for ExBG MC.
    if(nhitsvlayer<3&&nhitshlayer<3&&nhitsveto<3){return;} 
#endif

    //////// Action of hlayer Hits //////////
    if(hlhitcol){

      for(int l=0; l<nhitshlayer; l++){

        pid   = (*hlhitcol)[l]->GetParticle();
        detid = (*hlhitcol)[l]->GetDetID();
        trkid = (*hlhitcol)[l]->GetTrackID();
        edep  = (*hlhitcol)[l]->GetEdep();
	
        for(int i=0; i<3; i++){
	  pos[i]=((*hlhitcol)[l]->GetPosition())[i];
	}

        time  = (*hlhitcol)[l]->GetTime();
        mod   = (*hlhitcol)[l]->GetMod();
        pln   = (*hlhitcol)[l]->GetPln();
        view  = (*hlhitcol)[l]->GetView();
        ch    = (*hlhitcol)[l]->GetCh();
        pe    = (*hlhitcol)[l]->GetPE();
        lope  = (*hlhitcol)[l]->GetLOPE();

        gridcell_id_x1 = (*hlhitcol)[l]->GetCellidX1();
        gridcell_id_x2 = (*hlhitcol)[l]->GetCellidX2();
        gridcell_id_y1 = (*hlhitcol)[l]->GetCellidY1();
        gridcell_id_y2 = (*hlhitcol)[l]->GetCellidY2();
        
        double posx, posy, posz, posxy;
        detdim->GetPosInMod(mod,pln,view,ch,&posx,&posy,&posz);

        if(view==0){
	  posxy = posy;
	}
        else{
	  posxy = posx;
	}

        hitsum   -> Clear("C");
        hitsum   -> mod  = mod;
        hitsum   -> view = view;
        hitsum   -> pln  = pln;
        hitsum   -> ch   = ch;
        hitsum   -> xy   = posxy;
        hitsum   -> z    = posz;
        hitsum   -> time = time;
        hitsum   -> pe   = pe; 
        hitsum   -> lope = lope; 
        hitsum   -> gridcell_id_x1 = gridcell_id_x1;
        hitsum   -> gridcell_id_x2 = gridcell_id_x2;
        hitsum   -> gridcell_id_y1 = gridcell_id_y1;
        hitsum   -> gridcell_id_y2 = gridcell_id_y2;


        simhitsum-> Clear("C");
        simhitsum-> edeposit = edep;
        simhitsum-> trackid  = trkid;
        simhitsum-> pdg      = pid;
        
        if(hitsum -> pe > PE_THRESHOLD){
          assert( runaction->GetEvtSum()->AddSimHit(simhitsum) );
          nsimhits = runaction->GetEvtSum()->NSimHits();
          hitsum->AddSimHit( runaction->GetEvtSum()->GetSimHit( nsimhits-1 ) );
          runaction->GetEvtSum()->AddModHit(hitsum, mod, 4);
        }
        
#ifdef DEBUG_EVTACT
        G4cout << "\n=== hits in horizontal layer (cyan) ===\n";
        if( pe > PE_THRESHOLD ) 
          (*hlhitcol)[l]->Print();
#endif
        
        if( pe > PE_THRESHOLD ) 
          (*hlhitcol)[l]->Draw();
      }
    }
    

    //////// Action of vlayer Hits //////////
    if( vlhitcol ){

      for(int l=0; l<nhitsvlayer; l++){

        pid   = (*vlhitcol)[l]->GetParticle();
        detid = (*vlhitcol)[l]->GetDetID();
        trkid = (*vlhitcol)[l]->GetTrackID();
        edep  = (*vlhitcol)[l]->GetEdep();

        for(int i=0; i<3; i++){
	  pos[i]=((*vlhitcol)[l]->GetPosition())[i];
	}

        time  = (*vlhitcol)[l]->GetTime();
        mod   = (*vlhitcol)[l]->GetMod();
        pln   = (*vlhitcol)[l]->GetPln();
        view  = (*vlhitcol)[l]->GetView();
        ch    = (*vlhitcol)[l]->GetCh();
        pe    = (*vlhitcol)[l]->GetPE();
        lope  = (*vlhitcol)[l]->GetLOPE();

        gridcell_id_x1 = (*vlhitcol)[l]->GetCellidX1();
        gridcell_id_x2 = (*vlhitcol)[l]->GetCellidX2();
        gridcell_id_y1 = (*vlhitcol)[l]->GetCellidY1();
        gridcell_id_y2 = (*vlhitcol)[l]->GetCellidY2();

        double posx, posy, posz, posxy;
        detdim->GetPosInMod(mod,pln,view,ch,&posx,&posy,&posz);

        if(view==0){
	  posxy = posy;
	}
        else{
	  posxy = posx;
	}

        hitsum   -> Clear("C");
        hitsum   -> mod  = mod;
        hitsum   -> view = view;
        hitsum   -> pln  = pln;
        hitsum   -> ch   = ch;
        hitsum   -> xy   = posxy;
        hitsum   -> z    = posz;
        hitsum   -> time = time;
        hitsum   -> pe   = pe; 
        hitsum   -> lope = lope;
        hitsum   -> gridcell_id_x1 = gridcell_id_x1;
        hitsum   -> gridcell_id_x2 = gridcell_id_x2;
        hitsum   -> gridcell_id_y1 = gridcell_id_y1;
        hitsum   -> gridcell_id_y2 = gridcell_id_y2;

        simhitsum-> Clear("C");
        simhitsum-> edeposit = edep;
        simhitsum-> trackid  = trkid;
        simhitsum-> pdg      = pid;
        
        if(hitsum->pe > PE_THRESHOLD){
          assert( runaction->GetEvtSum()->AddSimHit(simhitsum) );
          nsimhits = runaction->GetEvtSum()->NSimHits();
          hitsum->AddSimHit( runaction->GetEvtSum()->GetSimHit(nsimhits-1) );
          runaction->GetEvtSum()->AddModHit(hitsum, mod, 4);
        }
        
        if(pe > PE_THRESHOLD){
          (*vlhitcol)[l]->Draw();
	}
      }
    }


    ////// Action of veto Hits ////////// 
    if(vehitcol){

      for(int l=0; l<nhitsveto; l++){

        pid   = (*vehitcol)[l]->GetParticle();
        detid = (*vehitcol)[l]->GetDetID();
        trkid = (*vehitcol)[l]->GetTrackID();
        edep  = (*vehitcol)[l]->GetEdep();

        for(int i=0; i<3; i++){
	  pos[i]=((*vehitcol)[l]->GetPosition())[i];
	}

        time  = (*vehitcol)[l]->GetTime();
        mod   = (*vehitcol)[l]->GetMod();
        pln   = (*vehitcol)[l]->GetPln();
        view  = (*vehitcol)[l]->GetView();
        ch    = (*vehitcol)[l]->GetCh();
        pe    = (*vehitcol)[l]->GetPE();
        lope  = (*vehitcol)[l]->GetLOPE();

        double posx, posy, posz, posxy;
        detdim->GetPosInMod(mod,pln,view,ch,&posx,&posy,&posz);
        if(view==0){
	  posxy = posy;
	}
        else{
	  posxy = posx;
	}

        hitsum   -> Clear("C");
        hitsum   -> mod  = mod;
        hitsum   -> view = view;
        hitsum   -> pln  = pln;
        hitsum   -> ch   = ch;
        hitsum   -> xy   = posxy;
        hitsum   -> z    = posz;
        hitsum   -> time = time;
        hitsum   -> pe   = pe ;
        hitsum   -> lope = lope ;

        simhitsum-> Clear("C");
        simhitsum-> edeposit = edep;
        simhitsum-> trackid  = trkid;
        simhitsum-> pdg      = pid;
        
        if( pe > VETO_PE_THRESHOLD) {
          assert( runaction->GetEvtSum()->AddSimHit( simhitsum ) );
          nsimhits = runaction->GetEvtSum()->NSimHits( );
          
          hitsum->AddSimHit( runaction->GetEvtSum()->GetSimHit( nsimhits-1 ) );
          runaction->GetEvtSum()->AddModHit( hitsum, mod, 4 );
        }
        
        // INGRID veto planes are shared in tow adjacent modules
        if( (0<=mod && mod<=5 && pln==12) ||
            (7<=mod && mod<=12 && pln==14) ) {
	  detdim->GetPosInMod(mod+1,pln-1,view,ch,&posx,&posy,&posz);
	  if(view==0) posxy = posy;
	  else        posxy = posx;     

	  hitsum -> mod = mod+1;
	  hitsum -> pln = pln-1;
	  hitsum -> xy  = posxy;
	  hitsum -> z   = posz;
                
	  if(pe > VETO_PE_THRESHOLD) {
	    runaction -> GetEvtSum() -> AddModHit(hitsum, mod+1, 4);
	  }
        }
        
        if(pe > PE_THRESHOLD){
          (*vehitcol)[l]->Draw();
	}
      }
    }
    
    //Extract the trajectories and draw them 
    //G4VVisManager* vis = G4VVisManager::GetConcreteInstance();

    //if(vis){

    //  for (G4int i=0; i<n_trajectories; i++){
    //     
    //    G4Trajectory* trj = (G4Trajectory*)
    //      ((*(anEvent->GetTrajectoryContainer()))[i]);
    //             
    //    if(trj->GetParentID()== 0){ //Particle created at neutrino interaction 
    //      //trj->DrawTrajectory(50);   
    //      trj->DrawTrajectory();   
    //    }
    //     
    //  } 

    //  //Draw vertex of neutrino interaction
    //  G4ThreeVector vertex; 

    //  for(int i=0; i<3; i++ ){
    //    vertex[i] = (runaction->vertex)[i]*cm;
    //  }

    //  G4Square square(vertex);
    //  square.SetScreenSize(6.5);
    //  square.SetFillStyle(G4Square::filled);
    //  G4Colour colour(1.,1.,1.); //White
    //  G4VisAttributes attribs(colour);
    //  square.SetVisAttributes(attribs);
    //  vis->Draw(square);
    //}
    
    //Fill Tree
    if(1){
      runaction->GetTree()  ->Fill();
      runaction->GetEvtSum()->Clear("C");
//#ifndef SUPPRESS_FOR_EXT_BG
//      runaction->GetSKTree()->Fill();//t2kreweight
//#endif
    }
    else{
      runaction->GetEvtSum()->Clear("C");
    }

  }
  
  else
    G4cout << "------- End of EventAction -->|\n" << G4endl;
  
}
