// -*- C++ -*-
//
// Package:    HLTrigger/MuonHLTSeedMVAClassifierPhase2
// Class:      MuonHLTSeedMVAClassifierPhase2
// 
/**\class MuonHLTSeedMVAClassifierPhase2 MuonHLTSeedMVAClassifierPhase2.cc HLTrigger/MuonHLTSeedMVAClassifierPhase2/plugins/MuonHLTSeedMVAClassifierPhase2.cc

 Description: [one line class summary]
*/
//
// Original Author:  OH Minseok
//         Created:  Mon, 08 Jun 2020 06:20:44 GMT
//
//


// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/ESInputTag.h"

// Geometry
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

// TrajectorySeed
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

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

#include "RecoMuon/TrackerSeedGenerator/interface/SeedMvaEstimatorPhase2.h"

// class declaration
bool sortByMvaScorePhase2(const std::pair<unsigned, double> &A, const std::pair<unsigned, double> &B) {
	return (A.second > B.second);
};

class MuonHLTSeedMVAClassifierPhase2 : public edm::stream::EDProducer<> {
	public:
		explicit MuonHLTSeedMVAClassifierPhase2(const edm::ParameterSet&);
		~MuonHLTSeedMVAClassifierPhase2() override = default;

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void produce(edm::Event&, const edm::EventSetup&) override;

		// ----------member data ---------------------------
		const edm::EDGetTokenT<TrajectorySeedCollection>             t_Seed_;
		const edm::EDGetTokenT<l1t::TrackerMuonCollection>           t_L1TkMu_;
		const edm::EDGetTokenT<reco::RecoChargedCandidateCollection> t_L2Muon_;

		typedef std::pair<std::unique_ptr<const SeedMvaEstimatorPhase2>, std::unique_ptr<const SeedMvaEstimatorPhase2>> pairSeedMvaEstimator;
		pairSeedMvaEstimator mvaEstimator_;

		const edm::FileInPath mvaFile_B_0_;
		const edm::FileInPath mvaFile_E_0_;

		const std::vector<double> mvaScaleMean_B_;
		const std::vector<double> mvaScaleStd_B_;
		const std::vector<double> mvaScaleMean_E_;
		const std::vector<double> mvaScaleStd_E_;

		const double etaEdge_;
		const double mvaCut_B_;
		const double mvaCut_E_;

		const bool doSort_;
		const int nSeedsMax_B_;
		const int nSeedsMax_E_;

		const double baseScore_;

		const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyESToken_;
		const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryESToken_;
		const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldESToken_;
		const edm::ESGetToken<GeometricDet, IdealGeometryRecord> geomDetESToken_;
		const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatorESToken_;

		double getSeedMva(const pairSeedMvaEstimator& pairMvaEstimator,
						  const TrajectorySeed& seed,
						  const GlobalVector& global_p,
						  const GlobalPoint&  global_x,
						  const edm::Handle<l1t::TrackerMuonCollection>& h_L1TkMu,
						  const edm::ESHandle<MagneticField>& magfieldH,
						  const Propagator& propagatorAlong,
						  const GeometricSearchTracker& geomTracker);
};

// constructors and destructor
MuonHLTSeedMVAClassifierPhase2::MuonHLTSeedMVAClassifierPhase2(const edm::ParameterSet& iConfig):
	t_Seed_(  consumes<TrajectorySeedCollection>            (iConfig.getParameter<edm::InputTag>("src"))),
	t_L1TkMu_(consumes<l1t::TrackerMuonCollection>          (iConfig.getParameter<edm::InputTag>("L1TkMu"))),
	t_L2Muon_(consumes<reco::RecoChargedCandidateCollection>(iConfig.getParameter<edm::InputTag>("L2Muon"))),

	mvaFile_B_0_   (iConfig.getParameter<edm::FileInPath>("mvaFile_B_0")),
	mvaFile_E_0_   (iConfig.getParameter<edm::FileInPath>("mvaFile_E_0")),

	mvaScaleMean_B_(iConfig.getParameter<std::vector<double>>("mvaScaleMean_B")),
	mvaScaleStd_B_ (iConfig.getParameter<std::vector<double>>("mvaScaleStd_B")),
	mvaScaleMean_E_(iConfig.getParameter<std::vector<double>>("mvaScaleMean_E")),
	mvaScaleStd_E_ (iConfig.getParameter<std::vector<double>>("mvaScaleStd_E")),

	etaEdge_       (iConfig.getParameter<double>("etaEdge")),
	mvaCut_B_      (iConfig.getParameter<double>("mvaCut_B")),
	mvaCut_E_      (iConfig.getParameter<double>("mvaCut_E")),

	doSort_        (iConfig.getParameter<bool>("doSort")),
	nSeedsMax_B_   (iConfig.getParameter<int>("nSeedsMax_B")),
	nSeedsMax_E_   (iConfig.getParameter<int>("nSeedsMax_E")),

	baseScore_     (iConfig.getParameter<double>("baseScore")),

	trackerTopologyESToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
	trackerGeometryESToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
	magFieldESToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
	geomDetESToken_(esConsumes<GeometricDet, IdealGeometryRecord>()),
	propagatorESToken_(esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("", "PropagatorWithMaterialParabolicMf")))
{
	produces<TrajectorySeedCollection>();

	mvaEstimator_ = std::make_pair(
		std::make_unique<SeedMvaEstimatorPhase2>(mvaFile_B_0_, mvaScaleMean_B_, mvaScaleStd_B_),
		std::make_unique<SeedMvaEstimatorPhase2>(mvaFile_E_0_, mvaScaleMean_E_, mvaScaleStd_E_) );
}

// member functions
// ------------ method called on each new Event  ------------
void MuonHLTSeedMVAClassifierPhase2::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	auto result = std::make_unique<TrajectorySeedCollection>();

	edm::ESHandle<TrackerGeometry> trkGeom = iSetup.getHandle(trackerGeometryESToken_);
	edm::ESHandle<GeometricDet> geomDet = iSetup.getHandle(geomDetESToken_);
	edm::ESHandle<TrackerTopology> trkTopo = iSetup.getHandle(trackerTopologyESToken_);

	GeometricSearchTrackerBuilder builder;
	GeometricSearchTracker geomTracker = *(builder.build(&(*geomDet), &(*trkGeom), &(*trkTopo)));

	edm::Handle<l1t::TrackerMuonCollection> h_L1TkMu;
	bool hasL1TkMu = iEvent.getByToken( t_L1TkMu_, h_L1TkMu);

	edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon;
	bool hasL2 = iEvent.getByToken( t_L2Muon_, h_L2Muon );

	edm::Handle<TrajectorySeedCollection> h_Seed;
	bool hasSeed = iEvent.getByToken( t_Seed_, h_Seed );

	edm::ESHandle<MagneticField> magfieldH = iSetup.getHandle(magFieldESToken_);
	edm::ESHandle<Propagator> propagatorAlongH = iSetup.getHandle(propagatorESToken_);
	std::unique_ptr<Propagator> propagatorAlong = SetPropagationDirection(*propagatorAlongH, alongMomentum);

	if( !( hasL1TkMu && hasL2 && hasSeed ) ) {
		std::cout << "MuonHLTSeedMVAClassifierPhase2::produce: !( hasL1 && hasL1TkMu && hasL2 && hasSeed )" << std::endl;
		return;
	}

	// -- sort seeds by MVA score and chooes top nSeedsMax_B_ / nSeedsMax_E_
	if( doSort_ ) {
		std::vector< std::pair<unsigned, double> > pairSeedIdxMvaScore_B = {};
		std::vector< std::pair<unsigned, double> > pairSeedIdxMvaScore_E = {};

		for( auto i=0U; i<h_Seed->size(); ++i ) {
			const auto& seed(h_Seed->at(i));

			GlobalVector global_p = trkGeom->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().momentum());
			GlobalPoint  global_x = trkGeom->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().position());

	
			bool isB = ( std::abs( global_p.eta() ) < etaEdge_ );

			if( isB && nSeedsMax_B_ < 0 ) {
				result->emplace_back( seed );
				continue;
			}

			if( !isB && nSeedsMax_E_ < 0 ) {
				result->emplace_back( seed );
				continue;
			}

			double mva = getSeedMva(
				mvaEstimator_,
				seed,
				global_p,
				global_x,
				h_L1TkMu,
				magfieldH,
				*(propagatorAlong.get()),
				geomTracker
			);

			double logistic = 1 / ( 1 + std::exp(-mva) );
			
			if(isB)  pairSeedIdxMvaScore_B.push_back( make_pair( i, logistic ) );
			else     pairSeedIdxMvaScore_E.push_back( make_pair( i, logistic ) );
		}

		std::sort(pairSeedIdxMvaScore_B.begin(), pairSeedIdxMvaScore_B.end(), sortByMvaScorePhase2 );
		std::sort(pairSeedIdxMvaScore_E.begin(), pairSeedIdxMvaScore_E.end(), sortByMvaScorePhase2 );

		for( auto i=0U; i<pairSeedIdxMvaScore_B.size(); ++i ) {
			if((int)i == nSeedsMax_B_)  break;
			const auto& seed(h_Seed->at( pairSeedIdxMvaScore_B.at(i).first ));
			result->emplace_back( seed );
		}

		for( auto i=0U; i<pairSeedIdxMvaScore_E.size(); ++i ) {
			if((int)i == nSeedsMax_E_)  break;
			const auto& seed(h_Seed->at( pairSeedIdxMvaScore_E.at(i).first ));
			result->emplace_back( seed );
		}
	}

	// -- simple fitering based on Mva threshold
	else {
		for( auto i=0U; i<h_Seed->size(); ++i ) {
			const auto& seed(h_Seed->at(i));

			GlobalVector global_p = trkGeom->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().momentum());
			GlobalPoint  global_x = trkGeom->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().position());

			bool isB = ( std::abs( global_p.eta() ) < etaEdge_ );

			if( isB && mvaCut_B_ <= 0. ) {
				result->emplace_back( seed );
				continue;
			}

			if( !isB && mvaCut_E_ <= 0. ) {
				result->emplace_back( seed );
				continue;
			}

			double mva = getSeedMva(
				mvaEstimator_,
				seed,
				global_p,
				global_x,
				h_L1TkMu,
				magfieldH,
				*(propagatorAlong.get()),
				geomTracker
			);

			double logistic = 1 / ( 1 + std::exp(-mva) );

			bool passMva = (
				(  isB && (logistic > mvaCut_B_) ) ||
				( !isB && (logistic > mvaCut_E_) )
			);

			if( passMva )  result->emplace_back( seed );
		}
	}

	iEvent.put(std::move(result));
}

double MuonHLTSeedMVAClassifierPhase2::getSeedMva(const pairSeedMvaEstimator& pairMvaEstimator,
												  const TrajectorySeed& seed,
												  const GlobalVector& global_p,
												  const GlobalPoint&  global_x,
												  const edm::Handle<l1t::TrackerMuonCollection>& h_L1TkMu,
												  const edm::ESHandle<MagneticField>& magfieldH,
												  const Propagator& propagatorAlong,
												  const GeometricSearchTracker& geomTracker) {
	double mva = 0.;

	if( fabs( global_p.eta() ) < etaEdge_ ) {
		mva = pairMvaEstimator.first->computeMva(
			seed,
			global_p,
			global_x,
			h_L1TkMu,
			magfieldH,
			propagatorAlong,
			geomTracker
		);
	}
	else {
		mva = pairMvaEstimator.second->computeMva(
			seed,
			global_p,
			global_x,
			h_L1TkMu,
			magfieldH,
			propagatorAlong,
			geomTracker
		);
	}
	
	return (baseScore_ + mva);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonHLTSeedMVAClassifierPhase2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonHLTSeedMVAClassifierPhase2);