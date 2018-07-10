//e+e- --> e+e-/mu+ mu-/pi+ pi-/g-g
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

// counter for efficiency
static long m_cout_all(0), m_cout_col(0), m_cout_charge(0), m_cout_nGood(0), m_cout_mom(0), m_cout_everat(0);
/////////////////////////////////////////////////////////////////////////////

TestAlgo::TestAlgo(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator) {
  //Declare the properties
  declareProperty("Vr0cut", m_vr0cut=1.0);
  declareProperty("Vz0cut", m_vz0cut=5.0);
  declareProperty("TrackCosThetaCut",m_cosThetaCut=0.93);
  declareProperty("eeDangCut",m_ee_cdang_cut=-0.98);

  declareProperty("CheckDedx", m_checkDedx = true);
  declareProperty("CheckTof",  m_checkTof = true);

  declareProperty("DistinEMuon", m_distin_emuon=2.0);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode TestAlgo::initialize(){
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

  NTuplePtr nt1(ntupleSvc(), "FILE1/tr");
  if ( nt1 ) m_tuple1 = nt1;
  else {
    m_tuple1 = ntupleSvc()->book ("FILE1/tr", CLID_ColumnWiseTuple, "processed data");
    if ( m_tuple1 )    {
      status = m_tuple1->addItem ("vx0",   m_vx0);
      status = m_tuple1->addItem ("vy0",   m_vy0);
      status = m_tuple1->addItem ("vz0",   m_vz0);
      status = m_tuple1->addItem ("vr0",   m_vr0);
      status = m_tuple1->addItem ("cosdang", m_cdang);
      status = m_tuple1->addItem ("q1", m_q1);
      status = m_tuple1->addItem ("p1", m_p1);
      status = m_tuple1->addItem ("cost1", m_cost1);
      status = m_tuple1->addItem ("Eemc1", m_Eemc1);
      status = m_tuple1->addItem ("q2", m_q2);
      status = m_tuple1->addItem ("p2", m_p2);
      status = m_tuple1->addItem ("cost2", m_cost2);
      status = m_tuple1->addItem ("Eemc2", m_Eemc2);
      status = m_tuple1->addItem ("eveflag", m_event_flag); //-2: charged trk; -1: e+e-; 1: mu+mu-; 2: pi+pi-; 3: gg
      status = m_tuple1->addItem ("run", m_run);
      status = m_tuple1->addItem ("event", m_event);
    }
    else    {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
      return StatusCode::FAILURE;
    }
  }

  //
  //--------end of book--------
  //

  log << MSG::INFO << "successfully return from initialize()" <<endmsg;
  return StatusCode::SUCCESS;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode TestAlgo::execute() {

  //std::cout << "execute()" << std::endl;

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;
  m_cout_all ++;
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



  m_cout_col ++;
  if(evtRecEvent->totalCharged()==2) m_event_flag=-2;
  else if(evtRecEvent->totalCharged()==0 && evtRecEvent->totalNeutral()==2) m_event_flag=3;
  else return sc;


  //if(evtRecEvent->totalCharged()<3 || evtRecTrkCol->size()<3 || evtRecEvent->totalTracks()>99 || evtRecTrkCol->size()>99) return StatusCode::SUCCESS;
  m_cout_charge ++;

  // Asign four-momentum with KalmanTrack
  Vint iGood; iGood.clear();
  int m_pid[2]={0,0}; // particle id: 0: noid/-2:{charged-,chaged+}/-1:{e-, e+}/1:{mu-, mu+}/2:{pi-, pi+}/3: {g, g}
  int nCharge = 0;
  HepLorentzVector m_lv1, m_lv2;

  if(m_event_flag==-2)
  {
    for(int i = 0; i < evtRecEvent->totalCharged(); i++)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
      if(!(*itTrk)->isMdcKalTrackValid()) return sc;
      RecMdcKalTrack* mdcTrk = (*itTrk)->mdcKalTrack();
      if(!(*itTrk)->isEmcShowerValid()) return sc;

      m_vx0 = mdcTrk->x();
      m_vy0 = mdcTrk->y();
      m_vz0 = mdcTrk->z();
      m_vr0 = mdcTrk->r();

      if(fabs(m_vz0) >= m_vz0cut)  return sc;
      if(m_vr0 >= m_vr0cut) return sc;
      iGood.push_back(i);
      nCharge += mdcTrk->charge();
      RecEmcShower *emcTrk = (*itTrk)->emcShower();

      if(mdcTrk->charge()>0)
      {
       m_q2 = mdcTrk->charge();
       m_p2 = mdcTrk->p();
       m_pid[1]=-2;
       m_Eemc2 = emcTrk->energy();
       if(m_Eemc2/m_p2<1.2 && m_Eemc2/m_p2>0.8)
       {
         m_pid[1]=-1;
         m_lv2=mdcTrk->p4(me);
         m_cost2 = m_lv2.vect().cosTheta();
       }
      }
      else if(mdcTrk->charge()<0)
      {
        m_q1 = mdcTrk->charge();
        m_p1 = mdcTrk->p();
        m_pid[0]=-2;
        m_Eemc1 = emcTrk->energy();
        if(m_Eemc1/m_p1<1.2 && m_Eemc1/m_p1>0.8)
        {
          m_pid[0]=-1;
          m_lv1=mdcTrk->p4(me);
          m_cost1 = m_lv1.vect().cosTheta();
        }
        /*
          distinguish mu+ mu-; pi+ pi-
        */
      }
      else return sc;
     }
  }
  else
  {

  }





  int nGood = iGood.size();
  log << MSG::DEBUG << "With KalmanTrack, ngood, totcharge = " << nGood << " , " << nCharge << endreq;
  if(nGood!=2 || nCharge) return sc;
  m_cout_nGood ++;



/*  double m_ep_ratio = 0;
  for(int i=0; i< evtRecEvent->totalTracks(); i++){
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(!(*itTrk)->isEmcShowerValid()) continue;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    m_ep_ratio += emcTrk->energy();
  }*/

//  if(m_ep_ratio < m_distin_emuon) return sc;


  // dangle between leptons
  m_cdang = m_lv1.vect().cosTheta(m_lv2.vect());
  m_run = run;
  m_event = event;



  if(m_cdang<m_ee_cdang_cut) setFilterPassed(true);
  //cout << "passed" << endl;
    fEvent.run = eventHeader->runNumber();
    fEvent.event = eventHeader->eventNumber();
    fEvent.time = eventHeader->time();
    fEvent.ntrack = 2;
    for(int i = 0; i <2; i++)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
      //fEvent.Pid.fill(i,*itTrk);
      if(eventHeader->runNumber() < 0)
      {
        fEvent.McTruth.fill(i,*itTrk,mcParticleCol);
      }
    }

  m_tuple1->write();
    return StatusCode::SUCCESS;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode TestAlgo::finalize() {
  return StatusCode::SUCCESS;
}
