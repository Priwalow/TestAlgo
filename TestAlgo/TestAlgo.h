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


  // declare energy/momentum to distinguish e and muon
  double m_distin_emuon;

  // define Ntuples here
  NTuple::Tuple*  m_tuple1;
  NTuple::Item<double>  m_vx0; // track vertex
  NTuple::Item<double>  m_vy0;
  NTuple::Item<double>  m_vz0;
  NTuple::Item<double>  m_vr0;
  NTuple::Item<double>  m_cdang;
  NTuple::Item<double>  m_q1; //charge
  NTuple::Item<double>  m_p1; //pulse
  NTuple::Item<double>  m_cost1; //cosTheta
  NTuple::Item<double>  m_Eemc1; //energy in electromagnetic calorimeter
  NTuple::Item<double>  m_q2;
  NTuple::Item<double>  m_p2;
  NTuple::Item<double>  m_cost2;
  NTuple::Item<double>  m_Eemc2;
  NTuple::Item<long> m_event_flag; //0: e+e-; 1: mu+mu-; 2: pi+pi-; 3: gg
  NTuple::Item<long> m_run;
  NTuple::Item<long> m_event;
};


#endif
