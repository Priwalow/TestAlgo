#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "trkInfo.h"

#include "RootEvent.h"

class TestAlgo : public Algorithm {

public:
  TestAlgo(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private:
  RootEvent fEvent;
  double EMC_ENDCUP_MIN_COS_THETA;
  double EMC_ENDCUP_MAX_COS_THETA;
  double EMC_ENDCUP_MIN_ENERGY;
  double EMC_BARREL_MAX_COS_THETA;
  double EMC_BARREL_MIN_ENERGY;
  // Declare r0, z0 cut for charged tracks
  //double m_vr0cut;
  //double m_vz0cut;

};
