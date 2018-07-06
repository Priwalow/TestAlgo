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
      status = m_tuple1->addItem ("eveflag", m_event_flag);
      status = m_tuple1->addItem ("eveflag", m_run);
      status = m_tuple1->addItem ("eveflag", m_event);
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
  if(evtRecEvent->totalCharged()<3 || evtRecTrkCol->size()<3 || evtRecEvent->totalTracks()>99 || evtRecTrkCol->size()>99) return StatusCode::SUCCESS;
  m_cout_charge ++;

  // Asign four-momentum with KalmanTrack
  Vint iGood; iGood.clear();
  int m_num[2]={0,0}; // number of different particles: e-, e+
  int nCharge = 0, m_lep_matched = 0;
  HepLorentzVector m_lv_ele, m_lv_pos;

  for(int i = 0; i < evtRecEvent->totalCharged(); i++)
  {
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(!(*itTrk)->isMdcKalTrackValid()) continue;
    RecMdcKalTrack* mdcTrk = (*itTrk)->mdcKalTrack();

    m_vx0 = mdcTrk->x();
    m_vy0 = mdcTrk->y();
    m_vz0 = mdcTrk->z();
    m_vr0 = mdcTrk->r();

    if(fabs(m_vz0) >= m_vz0cut) continue;
    if(m_vr0 >= m_vr0cut) continue;
    iGood.push_back(i);
    nCharge += mdcTrk->charge();
    if(mdcTrk->p()<1.0){ if((*itTrk)->isEmcShowerValid()) return sc; }
    else{ if((*itTrk)->isEmcShowerValid()) m_lep_matched ++; }

    if(mdcTrk->charge()>0)
    {
      if(mdcTrk->p()<1.0)
      {
	     return sc;
      }
      else
      {
       mdcTrk->setPidType(RecMdcKalTrack::electron);
	     m_lv_pos = mdcTrk->p4(xmass[0]);
	     m_num[1] ++;
      }
    }
    else
    {
      if(mdcTrk->p()<1.0)
      {
	     return sc;
      }
      else
      {
	     mdcTrk->setPidType(RecMdcKalTrack::electron);
	     m_lv_ele = mdcTrk->p4(xmass[0]);
	     m_num[0] ++;
      }
    }
  }

  int nGood = iGood.size();
  log << MSG::DEBUG << "With KalmanTrack, ngood, totcharge = " << nGood << " , " << nCharge << endreq;
  if(nGood<1 || nGood>2) return sc;
  m_cout_nGood ++;

  double m_ep_ratio = 0;
  for(int i=0; i< evtRecEvent->totalTracks(); i++){
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(!(*itTrk)->isEmcShowerValid()) continue;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    m_ep_ratio += emcTrk->energy();
  }

  if(m_ep_ratio < m_distin_emuon) return sc;

  HepLorentzVector m_lv_lab(0.04,0,0,3.686);
  if(nGood==2){
    if(nCharge) return sc;
    m_event_flag = 2;
  }
  else{
    if(m_num[0]>1 || m_num[1]>1) return sc;
    if(m_num[0]==0){
      if(nCharge == -1) return sc;
      m_lv_ele = m_lv_lab - m_lv_pos;
      if(m_lv_ele.vect().cosTheta()>m_cosThetaCut) return sc;
      m_event_flag = 0;
    }
    if(m_num[1]==0){
      if(nCharge == 1) return sc;
      m_lv_pos = m_lv_lab - m_lv_ele;
      if(m_lv_pos.vect().cosTheta()>m_cosThetaCut) return sc;
      m_event_flag = 1;
    }
  }
  m_cout_mom ++;

  // dangle between leptons
  m_elpos_cdang = m_lv_ele.vect().cosTheta(m_lv_pos.vect());
  m_run = run;
  m_event = event;


  if(m_subsample_flag) setFilterPassed(true);
  else if(m_elpos_cdang<m_ee_cdang_cut) setFilterPassed(true);
  //cout << "passed" << endl;
  m_tuple1->write();
  //MC information
/*  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
  if(m_run<0){
    int m_numParticle(0), m_true_pid(0);
    if(!mcParticleCol){
      log << MSG::ERROR << "Could not retrieve McParticelCol" << endreq;
      return StatusCode::FAILURE;
    }
    else{
      bool psipDecay(false);
      int rootIndex(-1);
      Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
      for (; iter_mc != mcParticleCol->end(); iter_mc++){
        if ((*iter_mc)->primaryParticle()) continue;
        if (!(*iter_mc)->decayFromGenerator()) continue;
        //if ( ((*iter_mc)->mother()).trackIndex()<3 ) continue;
        if ((*iter_mc)->particleProperty()==100443){
          psipDecay = true;
          rootIndex = (*iter_mc)->trackIndex();
        }
        if (!psipDecay) continue;
        int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
        int pdgid = (*iter_mc)->particleProperty();
        m_pdgid[m_numParticle] = pdgid;
        m_motheridx[m_numParticle] = mcidx;
        m_numParticle ++;

	//if(!(*iter_mc)->leafParticle()) continue;
	if((*iter_mc)->particleProperty() == 211) m_true_pionp = (*iter_mc)->initialFourMomentum().vect().mag();
	if((*iter_mc)->particleProperty() == -211) m_true_pionm = (*iter_mc)->initialFourMomentum().vect().mag();
      }
      m_idxmc = m_numParticle;
    }
  }


  m_tuple8->write();


  // next is good photon selection
  Vint iGam;  iGam.clear();
  for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) {
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(!(*itTrk)->isEmcShowerValid()) continue;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
    // find the nearest charged track
    double dthe = 200.;
    double dphi = 200.;
    double dang = 200.;
    for(int j = 0; j < evtRecEvent->totalCharged(); j++) {
      EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
      if(!(*jtTrk)->isExtTrackValid()) continue;
      RecExtTrack *extTrk = (*jtTrk)->extTrack();
      if(extTrk->emcVolumeNumber() == -1) continue;
      Hep3Vector extpos = extTrk->emcPosition();
      //      double ctht = extpos.cosTheta(emcpos);
      double angd = extpos.angle(emcpos);
      double thed = extpos.theta() - emcpos.theta();
      double phid = extpos.deltaPhi(emcpos);
      thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;

      if(fabs(thed) < fabs(dthe)) dthe = thed;
      if(fabs(phid) < fabs(dphi)) dphi = phid;
      if(angd < dang) dang = angd;
    }
    if(dang>=200) continue;
    double eraw = emcTrk->energy();
    dthe = dthe * 180 / (CLHEP::pi);
    dphi = dphi * 180 / (CLHEP::pi);
    dang = dang * 180 / (CLHEP::pi);
    m_dthe = dthe;
    m_dphi = dphi;
    m_dang = dang;
    m_eraw = eraw;
    // if(eraw < m_energyThreshold) continue;
    // if((fabs(dthe) < m_gammaThetaCut) && (fabs(dphi)<m_gammaPhiCut) ) continue;
    // good photon cut will be set here
    iGam.push_back((*itTrk)->trackId());
  }
  // Finish Good Photon Selection
  m_nGam = iGam.size();
  log << MSG::DEBUG << "num Good Photon " << m_nGam  << " , " <<evtRecEvent->totalNeutral()<<endreq;
  m_tuple2->write();

  //
  // check dedx infomation
  //
  if(m_checkDedx) {
    int m_dedx_cout(0);
    for(int i = 0; i < nGood; i++) {
      EvtRecTrackIterator  itTrk = evtRecTrkCol->begin() + iGood[i];
      if(!(*itTrk)->isMdcDedxValid())continue;
      RecMdcKalTrack *mdcTrk = (*itTrk)->mdcKalTrack();
      RecMdcDedx *dedxTrk = (*itTrk)->mdcDedx();

      m_ptrk = mdcTrk->p();
      m_chie = dedxTrk->chiE();
      m_chimu = dedxTrk->chiMu();
      m_chipi = dedxTrk->chiPi();
      m_chik = dedxTrk->chiK();
      m_chip = dedxTrk->chiP();
      m_ghit = dedxTrk->numGoodHits();
      m_thit = dedxTrk->numTotalHits();
      m_probPH = dedxTrk->probPH();
      m_normPH = dedxTrk->normPH();

      m_tuple3->write();
    }
  }

  //
  // check TOF infomation
  //
  if(m_checkTof) {
    int m_endcap_cout(0), m_layer1_cout(0), m_layer2_cout(0);
    for(int i = 0; i < nGood; i++) {
      EvtRecTrackIterator  itTrk = evtRecTrkCol->begin() + iGood[i];
      if(!(*itTrk)->isTofTrackValid()) continue;

      RecMdcKalTrack *mdcTrk = (*itTrk)->mdcKalTrack();
      SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();

      double ptrk = mdcTrk->p();

      for( SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin() ;iter_tof != tofTrkCol.end(); iter_tof++ ) {
        TofHitStatus *status = new TofHitStatus;
        status->setStatus((*iter_tof)->status());
        if(!(status->is_barrel())){//endcap
          if( !(status->is_counter()) ) continue; // ?
          if( status->layer()!=0 ) continue;//layer1
          double path = (*iter_tof)->path(); // ? the unit is cm
          double tof  = (*iter_tof)->tof();  // the unit is ns/100
          double ph   = (*iter_tof)->ph();
          double rhit = (*iter_tof)->zrhit();
          double qual = 0.0 + (*iter_tof)->quality();
          double cntr = 0.0 + (*iter_tof)->tofID();
          double texp[5];
          for(int j = 0; j < 5; j++) {
            double gb = xmass[j]/ptrk;
            double beta = sqrt(1+gb*gb);
            texp[j] = path*beta/velc; // the unit here is ns
          }
          m_cntr_etof  = cntr;
          m_ptot_etof  = ptrk;
	  m_path_etof = path;
          m_ph_etof    = ph;
          m_rhit_etof  = rhit;
          m_qual_etof  = qual;
	  m_tof_etof = tof;
          m_te_etof    = tof - texp[0];
          m_tmu_etof   = tof - texp[1];
          m_tpi_etof   = tof - texp[2];
          m_tk_etof    = tof - texp[3];
          m_tp_etof    = tof - texp[4];

          m_tuple4->write();
        }
        else {//barrel
          if( !(status->is_counter()) ) continue; // ?
          if(status->layer()==1){ //layer1
            double path=(*iter_tof)->path(); // ?
            double tof  = (*iter_tof)->tof();
            double ph   = (*iter_tof)->ph();
            double rhit = (*iter_tof)->zrhit();
            double qual = 0.0 + (*iter_tof)->quality();
            double cntr = 0.0 + (*iter_tof)->tofID();
            double texp[5];
            for(int j = 0; j < 5; j++) {
            double gb = xmass[j]/ptrk;
            double beta = sqrt(1+gb*gb);
            texp[j] = path*beta/velc;
            }

            m_cntr_btof1  = cntr;
            m_ptot_btof1 = ptrk;
	    m_path_btof1 = path;
            m_ph_btof1   = ph;
            m_zhit_btof1  = rhit;
            m_qual_btof1  = qual;
	    m_tof_btof1 = tof;
            m_te_btof1    = tof - texp[0];
            m_tmu_btof1   = tof - texp[1];
            m_tpi_btof1   = tof - texp[2];
            m_tk_btof1    = tof - texp[3];
            m_tp_btof1    = tof - texp[4];

            m_tuple5->write();
          }

          if(status->layer()==2){//layer2
            double path=(*iter_tof)->path(); // ?
            double tof  = (*iter_tof)->tof();
            double ph   = (*iter_tof)->ph();
            double rhit = (*iter_tof)->zrhit();
            double qual = 0.0 + (*iter_tof)->quality();
            double cntr = 0.0 + (*iter_tof)->tofID();
            double texp[5];
            for(int j = 0; j < 5; j++) {
            double gb = xmass[j]/ptrk;
            double beta = sqrt(1+gb*gb);
            texp[j] = path*beta/velc;
            }

            m_cntr_btof2  = cntr;
            m_ptot_btof2 = ptrk;
	    m_path_btof2 = path;
            m_ph_btof2   = ph;
            m_zhit_btof2  = rhit;
            m_qual_btof2  = qual;
	    m_tof_btof2 = tof;
            m_te_btof2    = tof - texp[0];
            m_tmu_btof2   = tof - texp[1];
            m_tpi_btof2   = tof - texp[2];
            m_tk_btof2    = tof - texp[3];
            m_tp_btof2    = tof - texp[4];

	    m_tuple6->write();
          }
        }

        delete status;
      }
    } // loop all charged track
  }  // check tof


   */
    return StatusCode::SUCCESS;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode TestAlgo::finalize() {
  return StatusCode::SUCCESS;
}
