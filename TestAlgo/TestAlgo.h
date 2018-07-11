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

  // Declare r0, z0 cut for charged tracks
  //double m_vr0cut;
  //double m_vz0cut;

};
