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

#include "HLTrigger/MuonHLTSeedMVAClassifierPhase2/interface/SeedMvaEstimator2.h"


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

// ------ For CMSSW_12 ------
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
// ---------------------------

//
// class declaration
//

bool sortByMvaScore(const std::pair<unsigned, float> &A, const std::pair<unsigned, float> &B) {
	return (A.second > B.second);
};

class MuonHLTSeedMVAClassifierPhase2 : public edm::stream::EDProducer<> {
	public:
		explicit MuonHLTSeedMVAClassifierPhase2(const edm::ParameterSet&);
		~MuonHLTSeedMVAClassifierPhase2();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void produce(edm::Event&, const edm::EventSetup&) override;
		virtual void beginJob();
		virtual void endJob();

		//virtual void beginStream(edm::StreamID) override;
		//virtual void endStream() override;
		//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
		//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

		// ----------member data ---------------------------
		edm::EDGetTokenT<TrajectorySeedCollection>             t_Seed_;
		// edm::EDGetTokenT<l1t::MuonBxCollection>                t_L1Muon_;
		edm::EDGetTokenT<l1t::TrackerMuonCollection>           t_L1TkMu_;
		edm::EDGetTokenT<reco::RecoChargedCandidateCollection> t_L2Muon_;

		typedef std::vector< std::pair<SeedMvaEstimatorPhase2*, SeedMvaEstimatorPhase2*> > pairSeedMvaEstimator;
		pairSeedMvaEstimator mvaEstimator;

		edm::FileInPath mvaFile_B_0_;
		// edm::FileInPath mvaFile_B_1_;
		// edm::FileInPath mvaFile_B_2_;
		// edm::FileInPath mvaFile_B_3_;
		edm::FileInPath mvaFile_E_0_;
		// edm::FileInPath mvaFile_E_1_;
		// edm::FileInPath mvaFile_E_2_;
		// edm::FileInPath mvaFile_E_3_;

		std::vector<double> mvaScaleMean_B_;
		std::vector<double> mvaScaleStd_B_;
		std::vector<double> mvaScaleMean_E_;
		std::vector<double> mvaScaleStd_E_;

		const double mvaCut_B_;
		const double mvaCut_E_;

		const bool doSort_;
		const int nSeedsMax_B_;
		const int nSeedsMax_E_;

		// ------ For CMSSW_12 -------
		const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyESToken_;
		const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryESToken_;
		const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldESToken_;
		const edm::ESGetToken<GeometricDet, IdealGeometryRecord> geomDetESToken_;
		const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatorESToken_;
		// ---------------------------

		std::vector<float> getSeedMva(
			pairSeedMvaEstimator pairMvaEstimator,
			const TrajectorySeed& seed,
			GlobalVector global_p,
			GlobalPoint  global_x,
			edm::Handle<l1t::TrackerMuonCollection>& h_L1TkMu,
			edm::ESHandle<MagneticField>& magfieldH,
			const Propagator& propagatorAlong,
			GeometricSearchTracker* geomTracker,
			float offset
		);
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonHLTSeedMVAClassifierPhase2::MuonHLTSeedMVAClassifierPhase2(const edm::ParameterSet& iConfig):
	t_Seed_(  consumes<TrajectorySeedCollection>            (iConfig.getParameter<edm::InputTag>("src"))),
	// t_L1Muon_(consumes<l1t::MuonBxCollection>               (iConfig.getParameter<edm::InputTag>("L1Muon"))),
	t_L1TkMu_(consumes<l1t::TrackerMuonCollection>          (iConfig.getParameter<edm::InputTag>("L1TkMu"))),
	t_L2Muon_(consumes<reco::RecoChargedCandidateCollection>(iConfig.getParameter<edm::InputTag>("L2Muon"))),

	mvaFile_B_0_   (iConfig.getParameter<edm::FileInPath>("mvaFile_B_0")),
	// mvaFile_B_1_   (iConfig.getParameter<edm::FileInPath>("mvaFile_B_1")),
	// mvaFile_B_2_   (iConfig.getParameter<edm::FileInPath>("mvaFile_B_2")),
	// mvaFile_B_3_   (iConfig.getParameter<edm::FileInPath>("mvaFile_B_3")),
	mvaFile_E_0_   (iConfig.getParameter<edm::FileInPath>("mvaFile_E_0")),
	// mvaFile_E_1_   (iConfig.getParameter<edm::FileInPath>("mvaFile_E_1")),
	// mvaFile_E_2_   (iConfig.getParameter<edm::FileInPath>("mvaFile_E_2")),
	// mvaFile_E_3_   (iConfig.getParameter<edm::FileInPath>("mvaFile_E_3")),

	mvaScaleMean_B_(iConfig.getParameter<std::vector<double>>("mvaScaleMean_B")),
	mvaScaleStd_B_ (iConfig.getParameter<std::vector<double>>("mvaScaleStd_B")),
	mvaScaleMean_E_(iConfig.getParameter<std::vector<double>>("mvaScaleMean_E")),
	mvaScaleStd_E_ (iConfig.getParameter<std::vector<double>>("mvaScaleStd_E")),

	mvaCut_B_      (iConfig.getParameter<double>("mvaCut_B")),
	mvaCut_E_      (iConfig.getParameter<double>("mvaCut_E")),

	doSort_        (iConfig.getParameter<bool>("doSort")),
	nSeedsMax_B_   (iConfig.getParameter<int>("nSeedsMax_B")),
	nSeedsMax_E_   (iConfig.getParameter<int>("nSeedsMax_E")),

	trackerTopologyESToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
	trackerGeometryESToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
	magFieldESToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
	geomDetESToken_(esConsumes<GeometricDet, IdealGeometryRecord>()),
	propagatorESToken_(esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("", "PropagatorWithMaterialParabolicMf")))
{
	produces<TrajectorySeedCollection>();

	mvaEstimator = {
		make_pair( new SeedMvaEstimatorPhase2(mvaFile_B_0_, mvaScaleMean_B_, mvaScaleStd_B_),
		           new SeedMvaEstimatorPhase2(mvaFile_E_0_, mvaScaleMean_E_, mvaScaleStd_E_) )//,
		// make_pair( new SeedMvaEstimatorPhase2(mvaFile_B_1_, mvaScaleMean_B_, mvaScaleStd_B_),
		//            new SeedMvaEstimatorPhase2(mvaFile_E_1_, mvaScaleMean_E_, mvaScaleStd_E_) ),
		// make_pair( new SeedMvaEstimatorPhase2(mvaFile_B_2_, mvaScaleMean_B_, mvaScaleStd_B_),
		//            new SeedMvaEstimatorPhase2(mvaFile_E_2_, mvaScaleMean_E_, mvaScaleStd_E_) ),
		// make_pair( new SeedMvaEstimatorPhase2(mvaFile_B_3_, mvaScaleMean_B_, mvaScaleStd_B_),
		//            new SeedMvaEstimatorPhase2(mvaFile_E_3_, mvaScaleMean_E_, mvaScaleStd_E_) )
	};
}


MuonHLTSeedMVAClassifierPhase2::~MuonHLTSeedMVAClassifierPhase2()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called on each new Event  ------------
void MuonHLTSeedMVAClassifierPhase2::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	auto result = std::make_unique<TrajectorySeedCollection>();

	// edm::ESHandle<TrackerGeometry> trkGeom;
	// iSetup.get<TrackerDigiGeometryRecord>().get(trkGeom);
	edm::ESHandle<TrackerGeometry> trkGeom = iSetup.getHandle(trackerGeometryESToken_);

	// edm::ESHandle<GeometricDet> geomDet;
	// iSetup.get<IdealGeometryRecord>().get(geomDet);
	edm::ESHandle<GeometricDet> geomDet = iSetup.getHandle(geomDetESToken_);
	
	// edm::ESHandle<TrackerTopology> trkTopo;
	// iSetup.get<TrackerTopologyRcd>().get(trkTopo);
	edm::ESHandle<TrackerTopology> trkTopo = iSetup.getHandle(trackerTopologyESToken_);

	GeometricSearchTrackerBuilder builder;
	GeometricSearchTracker* geomTracker = builder.build(&(*geomDet), &(*trkGeom), &(*trkTopo));

	// edm::Handle<l1t::MuonBxCollection> h_L1Muon;
	// bool hasL1 = iEvent.getByToken( t_L1Muon_, h_L1Muon);

	edm::Handle<l1t::TrackerMuonCollection> h_L1TkMu;
	bool hasL1TkMu = iEvent.getByToken( t_L1TkMu_, h_L1TkMu);

	edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon;
	bool hasL2 = iEvent.getByToken( t_L2Muon_, h_L2Muon );

	edm::Handle<TrajectorySeedCollection> h_Seed;
	bool hasSeed = iEvent.getByToken( t_Seed_, h_Seed );

	std::cout << "hasSeed = " << hasSeed << std::endl;

	// edm::ESHandle<MagneticField> magfieldH;
	// iSetup.get<IdealMagneticFieldRecord>().get(magfieldH);
	edm::ESHandle<MagneticField> magfieldH = iSetup.getHandle(magFieldESToken_);

	// edm::ESHandle<Propagator> propagatorAlongH;
    // iSetup.get<TrackingComponentsRecord>().get("PropagatorWithMaterialParabolicMf", propagatorAlongH);
	edm::ESHandle<Propagator> propagatorAlongH = iSetup.getHandle(propagatorESToken_);
	std::unique_ptr<Propagator> propagatorAlong = SetPropagationDirection(*propagatorAlongH, alongMomentum);

	// if( !( hasL1 && hasL1TkMu && hasL2 && hasSeed ) ) {
	if( !( hasL1TkMu && hasL2 && hasSeed ) ) {
		std::cout << "MuonHLTSeedMVAClassifierPhase2::produce: !( hasL1 && hasL1TkMu && hasL2 && hasSeed )" << std::endl;
		return;
	}

	// -- sort seeds by MVA score and chooes top nSeedsMax_B_ / nSeedsMax_E_
	if( doSort_ ) {
		std::vector< std::pair<unsigned, float> > pairSeedIdxMvaScore_B = {};
		std::vector< std::pair<unsigned, float> > pairSeedIdxMvaScore_E = {};

		for( auto i=0U; i<h_Seed->size(); ++i ) {
			const auto& seed(h_Seed->at(i));

			GlobalVector global_p = trkGeom->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().momentum());
			GlobalPoint  global_x = trkGeom->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().position());

			// FIXME this should be configurable
			bool isB = ( std::abs( global_p.eta() ) < 0.9 );

			if( isB && nSeedsMax_B_ < 0 ) {
				result->emplace_back( seed );
				continue;
			}

			if( !isB && nSeedsMax_E_ < 0 ) {
				result->emplace_back( seed );
				continue;
			}

			std::vector<float> mvas = getSeedMva(
				mvaEstimator,
				seed,
				global_p,
				global_x,
				h_L1TkMu,
				magfieldH,
				*(propagatorAlong.get()),
				geomTracker,
				0.5
			);

			float logistic = 1 / ( 1 + std::exp(-mvas.at(0)) );
			
			if(isB)  pairSeedIdxMvaScore_B.push_back( make_pair( i, logistic ) );
			else     pairSeedIdxMvaScore_E.push_back( make_pair( i, logistic ) );
		}

		std::sort(pairSeedIdxMvaScore_B.begin(), pairSeedIdxMvaScore_B.end(), sortByMvaScore );
		std::sort(pairSeedIdxMvaScore_E.begin(), pairSeedIdxMvaScore_E.end(), sortByMvaScore );

		for( auto i=0U; i<pairSeedIdxMvaScore_B.size(); ++i ) {
			if((int)i == nSeedsMax_B_)  break;
			// std::cout << "B: " << i << " idx=" << pairSeedIdxMvaScore_B.at(i).first << " mva=" << pairSeedIdxMvaScore_B.at(i).second << std::endl;
			const auto& seed(h_Seed->at( pairSeedIdxMvaScore_B.at(i).first ));
			result->emplace_back( seed );
		}

		for( auto i=0U; i<pairSeedIdxMvaScore_E.size(); ++i ) {
			if((int)i == nSeedsMax_E_)  break;
			// std::cout << "E: " << i << " idx=" << pairSeedIdxMvaScore_E.at(i).first << " mva=" << pairSeedIdxMvaScore_E.at(i).second << std::endl;
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

			// FIXME this should be configurable
			bool isB = ( std::abs( global_p.eta() ) < 0.9 );

			if( isB && mvaCut_B_ <= 0. ) {
				result->emplace_back( seed );
				continue;
			}

			if( !isB && mvaCut_E_ <= 0. ) {
				result->emplace_back( seed );
				continue;
			}

			std::vector<float> mvas = getSeedMva(
				mvaEstimator,
				seed,
				global_p,
				global_x,
				h_L1TkMu,
				magfieldH,
				*(propagatorAlong.get()),
				geomTracker,
				0.5
			);

			float logistic = 1 / ( 1 + std::exp(-mvas.at(0)) );

			bool passMva = (
				(  isB && (logistic > mvaCut_B_) ) ||
				( !isB && (logistic > mvaCut_E_) )
			);

			if( passMva )  result->emplace_back( seed );
		}
	}

	iEvent.put(std::move(result));
}

std::vector<float> MuonHLTSeedMVAClassifierPhase2::getSeedMva(
	pairSeedMvaEstimator pairMvaEstimator,
	const TrajectorySeed& seed,
	GlobalVector global_p,
	GlobalPoint  global_x,
	edm::Handle<l1t::TrackerMuonCollection>& h_L1TkMu,
	edm::ESHandle<MagneticField>& magfieldH,
	const Propagator& propagatorAlong,
	GeometricSearchTracker* geomTracker,
	float offset = 0.5
) {
	std::vector<float> v_mva = {};

	for(auto ic=0U; ic<pairMvaEstimator.size(); ++ic) {
		// FIXME this should be configurable
		if( fabs( global_p.eta() ) < 0.9 ) {
			float mva = pairMvaEstimator.at(ic).first->computeMva(
				seed,
				global_p,
				global_x,
				h_L1TkMu,
				magfieldH,
				propagatorAlong,
				geomTracker
			);
			v_mva.push_back( (offset + mva) );
		}
		else {
			float mva = pairMvaEstimator.at(ic).second->computeMva(
				seed,
				global_p,
				global_x,
				h_L1TkMu,
				magfieldH,
				propagatorAlong,
				geomTracker
			);
			v_mva.push_back( (offset + mva) );
		}
	}
	// if( v_mva.size() != 4 ) {  // this should never happen
	// 	std::cout << "MuonHLTSeedMVAClassifierPhase2::getSeedMva: v_mva.size() != 4" << std::endl;
	// 	return { -99999., -99999., -99999., -99999. };
	// }
		if( v_mva.size() != 1 ) {  // this should never happen
		std::cout << "MuonHLTSeedMVAClassifierPhase2::getSeedMva: v_mva.size() != 1" << std::endl;
		return { -99999. };
	}

	return v_mva;
}

void MuonHLTSeedMVAClassifierPhase2::beginJob(){}

void MuonHLTSeedMVAClassifierPhase2::endJob(){
	// for( int i=0; i<4; ++i ) {
	for( int i=0; i<1; ++i ) {
		delete mvaEstimator.at(i).first;
		delete mvaEstimator.at(i).second;
	}
}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
/*
void MuonHLTSeedMVAClassifierPhase2::beginStream(edm::StreamID){}
*/

// ------------ method called once each stream after processing all runs, lumis and events  ------------
/*
void MuonHLTSeedMVAClassifierPhase2::endStream(){}
*/

// ------------ method called when starting to processes a run  ------------
/*
void MuonHLTSeedMVAClassifierPhase2::beginRun(edm::Run const&, edm::EventSetup const&){}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void MuonHLTSeedMVAClassifierPhase2::endRun(edm::Run const&, edm::EventSetup const&){}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MuonHLTSeedMVAClassifierPhase2::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MuonHLTSeedMVAClassifierPhase2::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonHLTSeedMVAClassifierPhase2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MuonHLTSeedMVAClassifierPhase2);
