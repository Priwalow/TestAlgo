//e+e- --> e+e-
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

  declareProperty("Subsample", m_subsample_flag=false);
  declareProperty("Trigger", m_trigger_flag=false);
  declareProperty("DistinEMuon", m_distin_emuon=2.0);

  declareProperty("EventRate", m_eventrate=false);
  declareProperty("ChanDet", m_chan_det=1);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode TestAlgo::initialize(){
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "in initialize()" << endmsg;

  StatusCode status;

  NTuplePtr nt1(ntupleSvc(), "FILE1/bbsc");
  if ( nt1 ) m_tuple1 = nt1;
  else {
    m_tuple1 = ntupleSvc()->book ("FILE1/bbsc", CLID_ColumnWiseTuple, "processed data");
    if ( m_tuple1 )    {
      status = m_tuple1->addItem ("vx0",   m_vx0);
      status = m_tuple1->addItem ("vy0",   m_vy0);
      status = m_tuple1->addItem ("vz0",   m_vz0);
      status = m_tuple1->addItem ("vr0",   m_vr0);
      status = m_tuple1->addItem ("cosdang", m_elpos_cdang);
      status = m_tuple1->addItem ("elcharge", m_el_charge);
      status = m_tuple1->addItem ("elp", m_el_p);
      status = m_tuple1->addItem ("elcost", m_el_cTheta);
      status = m_tuple1->addItem ("elEemc", m_el_Eemc);
      status = m_tuple1->addItem ("poscharge", m_pos_charge);
      status = m_tuple1->addItem ("posp", m_pos_p);
      status = m_tuple1->addItem ("poscost", m_pos_cTheta);
      status = m_tuple1->addItem ("posEemc", m_pos_Eemc);
      status = m_tuple1->addItem ("eveflag", m_event_flag);
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

  if(m_trigger_flag)
  {
    SmartDataPtr<TrigData> trigData(eventSvc(),EventModel::Trig::TrigData);
    if (!trigData) {
      log << MSG::FATAL << "Could not find Trigger Data for physics analysis" << endreq;
      return StatusCode::FAILURE;
    }
    /// Print trigger information once:
    log << MSG::DEBUG << "Trigger conditions: " << endreq;
    for(int i=0; i < 48; i++){
      log << MSG::DEBUG << "i:" << i << "  name:" << trigData->getTrigCondName(i) << "  cond:" << trigData->getTrigCondition(i) << endreq;
    }
    // test event rate
    int m_trig_tot(0), m_trig_which(-1);
    if(m_eventrate){
      for(int j=0; j<16; j++){
	       if(trigData->getTrigChannel(j)){
	          m_trig_tot ++;
	          m_trig_which = j+1;
	       }
      }
      if(m_trig_tot==1 && m_trig_which==m_chan_det) m_cout_everat++;
      return sc;
    }
  }

  m_cout_col ++;
  if(evtRecEvent->totalCharged()!=2) return StatusCode::SUCCESS;
  //if(evtRecEvent->totalCharged()<3 || evtRecTrkCol->size()<3 || evtRecEvent->totalTracks()>99 || evtRecTrkCol->size()>99) return StatusCode::SUCCESS;
  m_cout_charge ++;

  // Asign four-momentum with KalmanTrack
  Vint iGood; iGood.clear();
  int m_num[2]={0,0}; // number of different particles: e-, e+
  int nCharge = 0, m_lep_matched = 0;
  HepLorentzVector m_lv_ele, m_lv_pos;

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
       m_pos_charge = mdcTrk->charge();
       m_pos_p = mdcTrk->p();
       m_pos_cTheta = m_lv_pos.vect().cosTheta();
       m_pos_Eemc = emcTrk->energy();
	     m_num[1] ++;
    }
    else if(mdcTrk->charge()<0)
    {
       m_el_charge = mdcTrk->charge();
       m_el_p = mdcTrk->p();
       m_el_cTheta = m_lv_ele.vect().cosTheta();
       m_el_Eemc = emcTrk->energy();
	     m_num[0] ++;
    }
    else return sc;
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
  m_elpos_cdang = m_lv_ele.vect().cosTheta(m_lv_pos.vect());
  m_run = run;
  m_event = event;



  if(m_subsample_flag) setFilterPassed(true);
  else if(m_elpos_cdang<m_ee_cdang_cut) setFilterPassed(true);
  //cout << "passed" << endl;

  m_tuple1->write();
    return StatusCode::SUCCESS;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode TestAlgo::finalize() {
  return StatusCode::SUCCESS;
}
