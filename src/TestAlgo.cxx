//e+e- --> e+e-/mu+ mu-/pi+ pi-/g g
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"

#include "EventModel/EventModel.h"
#include "EventModel/EventHeader.h"
#include "EventModel/Event.h"
#include "TrigEvent/TrigEvent.h"
#include "TrigEvent/TrigData.h"
#include "McTruth/McParticle.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"

#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/WTrackParameter.h"
#include "ParticleID/ParticleID.h"
#include "MdcRecEvent/RecMdcKalTrack.h"

#include "../TestAlgo/TestAlgo.h"
#include "../TestAlgo/Utils.h"
#include <vector>
//const double twopi = 6.2831853;

const double me  = 0.000511;
const double mpi = 0.13957;
const double mproton = 0.938272;
const double mmu = 0.105658;
const double mpsip = 3.686;
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
const double velc = 29.9792458;  // tof_path unit in cm.
const double PI = 3.1415926;
// const double velc = 299.792458;   // tof path unit in mm
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;

/////////////////////////////////////////////////////////////////////////////

TestAlgo::TestAlgo(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {
  //Declare the properties
  //declareProperty("Vr0cut", m_vr0cut=1.0);
  //declareProperty("Vz0cut", m_vz0cut=5.0);
  /*declareProperty("EMC_ENDCUP_MIN_COS_THETA", EMC_ENDCUP_MIN_COS_THETA=0.86);
  declareProperty("EMC_ENDCUP_MAX_COS_THETA", EMC_ENDCUP_MAX_COS_THETA=0.92);
  declareProperty("EMC_ENDCUP_MIN_ENERGY", EMC_ENDCUP_MIN_ENERGY=0.05);
  declareProperty("EMC_BARREL_MAX_COS_THETA", EMC_BARREL_MAX_COS_THETA=0.8);
  declareProperty("EMC_BARREL_MIN_ENERGY", EMC_BARREL_MIN_ENERGY=0.025);*/
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode TestAlgo::initialize()
{
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "in initialize()" << endmsg;

  StatusCode status;
  try
  {
    fEvent.make_tuple(this, "FILE1/event","Some signal");
  }
  catch(std::runtime_error & error)
  {
    log << MSG::ERROR << error.what() << endmsg;
    return StatusCode::FAILURE;
  }
  //
  //--------end of book--------
  //

  log << MSG::INFO << "successfully return from initialize()" <<endmsg;
  return StatusCode::SUCCESS;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode TestAlgo::execute()
{

  //std::cout << "execute()" << std::endl;

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;
  StatusCode sc=StatusCode::SUCCESS;
  //save the events passed selection to a new file
  setFilterPassed(false);

  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),  EventModel::MC::McParticleCol);
  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
  if(!eventHeader){
    log << MSG::ERROR << "EventHeader not found" << endreq;
    return StatusCode::SUCCESS;
  }
  int run = eventHeader->runNumber();
  int event = eventHeader->eventNumber();
  //if(event%1000==0) cout << "run: " << run << " event: " << event << endl;
  if(event%1000==0) cout << "run: " << run << " event: " << event << endl;

  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  if(!evtRecEvent){
    log << MSG::ERROR << "EvtRecEvent not found" << endreq;
    return StatusCode::SUCCESS;
  }
  log << MSG::DEBUG <<"ncharg, nneu, tottks = "
  << evtRecEvent->totalCharged() << " , "
  << evtRecEvent->totalNeutral() << " , "
  << evtRecEvent->totalTracks() <<endreq;

  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  if(!evtRecTrkCol){
    log << MSG::ERROR << "EvtRecTrackCol" << endreq;
    return StatusCode::SUCCESS;
  }

  if(evtRecEvent->totalCharged()!=0 || evtRecEvent->totalTracks()!=2) return sc;

  SelectionConfig cfg;
  cfg.EMC_ENDCUP_MIN_COS_THETA=0.86;
  cfg.EMC_ENDCUP_MAX_COS_THETA=0.92;
  cfg.EMC_ENDCUP_MIN_ENERGY=0.05;
  cfg.EMC_BARREL_MAX_COS_THETA=0.8;
  cfg.EMC_BARREL_MIN_ENERGY=0.025;

  std::list<EvtRecTrack*> nGood = createGoodNeutralTrackList(cfg,evtRecEvent,evtRecTrkCol);
  if (nGood.size()!=2) return sc;

  fEvent.run = eventHeader->runNumber();
  fEvent.event = eventHeader->eventNumber();
  fEvent.time = eventHeader->time();
  fEvent.ntrack = 2;

  for(int i = 0; i <2; i++)
  {
    std::list<EvtRecTrack*>::iterator itTrk=nGood.begin();
    std::advance(itTrk, i);
    fEvent.fill(i,*itTrk);
    //fEvent.Pid.fill(i,*itTrk);
    /*if(eventHeader->runNumber() < 0)
    {
      fEvent.McTruth.fill(i,*itTrk,mcParticleCol);
    }*/
  }
  fEvent.write();

  setFilterPassed(true);
  return StatusCode::SUCCESS;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode TestAlgo::finalize() {
  return StatusCode::SUCCESS;
}
