#ifndef __PhysicsTools_PatAlgos_SeedMvaEstimator__
#define __PhysicsTools_PatAlgos_SeedMvaEstimator__

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/ESInputTag.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/L1TCorrelator/interface/TkMuon.h"
#include "DataFormats/L1TCorrelator/interface/TkMuonFwd.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"

#include <memory>
#include <string>

// For Phase 2 variables
// -- for L1TkMu propagation
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTrackerBuilder.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

#include "MuonHLTTool/MuonHLTNtupler/interface/MuonHLTobjCorrelator.h"

typedef pair<const DetLayer*, TrajectoryStateOnSurface> LayerTSOS;
typedef pair<const DetLayer*, const TrackingRecHit*> LayerHit;

class GBRForest;

namespace edm {
  class FileInPath;
}

// Phase 2 SeedMvaEstimator
class SeedMvaEstimatorPhase2 {
public:
  SeedMvaEstimatorPhase2( const edm::FileInPath& weightsfile, std::vector<double> scale_mean, std::vector<double> scale_std );
  ~SeedMvaEstimatorPhase2();

  float computeMva( const TrajectorySeed&,
    GlobalVector,
    GlobalPoint,
    edm::Handle<l1t::TrackerMuonCollection>,
    edm::ESHandle<MagneticField>&,
    const Propagator&,
    GeometricSearchTracker*
  ) const;

private:
  std::unique_ptr<const GBRForest> gbrForest_;

  std::vector<double> scale_mean_;
  std::vector<double> scale_std_;

  vector< LayerTSOS > getTsosOnPixels(
    TTTrack<Ref_Phase2TrackerDigi_>,
    edm::ESHandle<MagneticField>&,
    const Propagator&,
    GeometricSearchTracker*
  ) const;

  vector< pair<LayerHit, LayerTSOS> > getHitTsosPairs(
    const TrajectorySeed&,
    edm::Handle<l1t::TrackerMuonCollection>,
    edm::ESHandle<MagneticField>&,
    const Propagator&,
    GeometricSearchTracker*
  ) const;

  void getL1TTVariables(   const TrajectorySeed&, GlobalVector, GlobalPoint, edm::Handle<l1t::TrackerMuonCollection>, float&, float& ) const;
  void getHitL1TkVatiables( const TrajectorySeed&, edm::Handle<l1t::TrackerMuonCollection>, edm::ESHandle<MagneticField>&, const Propagator&, GeometricSearchTracker*, float&, float&, float&) const;
};
#endif
