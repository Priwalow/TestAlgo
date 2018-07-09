//e+e- --> e+e-
#ifndef Physics_Analysis_PipiJpsi_H
#define Physics_Analysis_PipiJpsi_H

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "trkInfo.h"

class TestAlgo : public Algorithm {

public:
  TestAlgo(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private:

  // Declare r0, z0 cut for charged tracks
  double m_vr0cut;
  double m_vz0cut;

  //
  bool m_checkDedx;
  bool m_checkTof;
  bool m_eventrate;
  int m_chan_det;
  // declare the track angle cut and track's transverse momentum
  double m_cosThetaCut;

  // declare the cut angle between two leptons (e+ e-)
   double m_ee_cdang_cut;

  // declare whether pick up sub-sample
  bool m_subsample_flag, m_trigger_flag;

  // declare energy/momentum to distinguish e and muon
  double m_distin_emuon;

  // define Ntuples here
  NTuple::Tuple*  m_tuple1;      // charged track vertex
  NTuple::Item<double>  m_vx0;
  NTuple::Item<double>  m_vy0;
  NTuple::Item<double>  m_vz0;
  NTuple::Item<double>  m_vr0;
  NTuple::Item<double>  m_elpos_cdang;
  NTuple::Item<double>  m_el_charge;
  NTuple::Item<double>  m_el_p;
  NTuple::Item<double>  m_el_cTheta;
  NTuple::Item<double>  m_el_Eemc;
  NTuple::Item<double>  m_pos_charge;
  NTuple::Item<double>  m_pos_p;
  NTuple::Item<double>  m_pos_cTheta;
  NTuple::Item<double>  m_pos_Eemc;
  NTuple::Item<long> m_event_flag;
  NTuple::Item<long> m_run;
  NTuple::Item<long> m_event;
};


#endif
