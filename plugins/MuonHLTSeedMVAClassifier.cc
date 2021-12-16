// -*- C++ -*-
//
// Package:    HLTrigger/MuonHLTSeedMVAClassifier
// Class:      MuonHLTSeedMVAClassifier
// 
/**\class MuonHLTSeedMVAClassifier MuonHLTSeedMVAClassifier.cc HLTrigger/MuonHLTSeedMVAClassifier/plugins/MuonHLTSeedMVAClassifier.cc

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

#include "HLTrigger/MuonHLTSeedMVAClassifier/interface/SeedMvaEstimator.h"

//
// class declaration
//

bool sortByMvaScore(const std::pair<unsigned, float> &A, const std::pair<unsigned, float> &B) {
	return (A.second > B.second);
};

class MuonHLTSeedMVAClassifier : public edm::stream::EDProducer<> {
	public:
		explicit MuonHLTSeedMVAClassifier(const edm::ParameterSet&);
		~MuonHLTSeedMVAClassifier();

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
		edm::EDGetTokenT<l1t::TkMuonCollection>                t_L1TkMu_;
		edm::EDGetTokenT<reco::RecoChargedCandidateCollection> t_L2Muon_;

		typedef std::vector< std::pair<SeedMvaEstimator*, SeedMvaEstimator*> > pairSeedMvaEstimator;
		pairSeedMvaEstimator mvaEstimator;

		edm::FileInPath mvaFile_B_0_;
		edm::FileInPath mvaFile_B_1_;
		edm::FileInPath mvaFile_B_2_;
		edm::FileInPath mvaFile_B_3_;
		edm::FileInPath mvaFile_E_0_;
		edm::FileInPath mvaFile_E_1_;
		edm::FileInPath mvaFile_E_2_;
		edm::FileInPath mvaFile_E_3_;

		std::vector<double> mvaScaleMean_B_;
		std::vector<double> mvaScaleStd_B_;
		std::vector<double> mvaScaleMean_E_;
		std::vector<double> mvaScaleStd_E_;

		const double mvaCut_B_;
		const double mvaCut_E_;

		const bool doSort_;
		const int nSeedsMax_B_;
		const int nSeedsMax_E_;

		std::vector<float> getSeedMva(
			pairSeedMvaEstimator pairMvaEstimator,
			const TrajectorySeed& seed,
			GlobalVector global_p,
			GlobalPoint  global_x,
			// edm::Handle<l1t::MuonBxCollection>& h_L1Muon,
			edm::Handle<l1t::TkMuonCollection>& h_L1TkMu,
			edm::Handle<reco::RecoChargedCandidateCollection>& h_L2Muon,
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
MuonHLTSeedMVAClassifier::MuonHLTSeedMVAClassifier(const edm::ParameterSet& iConfig):
	t_Seed_(  consumes<TrajectorySeedCollection>            (iConfig.getParameter<edm::InputTag>("src"))),
	// t_L1Muon_(consumes<l1t::MuonBxCollection>               (iConfig.getParameter<edm::InputTag>("L1Muon"))),
	t_L1TkMu_(consumes<l1t::TkMuonCollection>               (iConfig.getParameter<edm::InputTag>("L1TkMu"))),
	t_L2Muon_(consumes<reco::RecoChargedCandidateCollection>(iConfig.getParameter<edm::InputTag>("L2Muon"))),

	mvaFile_B_0_   (iConfig.getParameter<edm::FileInPath>("mvaFile_B_0")),
	mvaFile_B_1_   (iConfig.getParameter<edm::FileInPath>("mvaFile_B_1")),
	mvaFile_B_2_   (iConfig.getParameter<edm::FileInPath>("mvaFile_B_2")),
	mvaFile_B_3_   (iConfig.getParameter<edm::FileInPath>("mvaFile_B_3")),
	mvaFile_E_0_   (iConfig.getParameter<edm::FileInPath>("mvaFile_E_0")),
	mvaFile_E_1_   (iConfig.getParameter<edm::FileInPath>("mvaFile_E_1")),
	mvaFile_E_2_   (iConfig.getParameter<edm::FileInPath>("mvaFile_E_2")),
	mvaFile_E_3_   (iConfig.getParameter<edm::FileInPath>("mvaFile_E_3")),

	mvaScaleMean_B_(iConfig.getParameter<std::vector<double>>("mvaScaleMean_B")),
	mvaScaleStd_B_ (iConfig.getParameter<std::vector<double>>("mvaScaleStd_B")),
	mvaScaleMean_E_(iConfig.getParameter<std::vector<double>>("mvaScaleMean_E")),
	mvaScaleStd_E_ (iConfig.getParameter<std::vector<double>>("mvaScaleStd_E")),

	mvaCut_B_      (iConfig.getParameter<double>("mvaCut_B")),
	mvaCut_E_      (iConfig.getParameter<double>("mvaCut_E")),

	doSort_        (iConfig.getParameter<bool>("doSort")),
	nSeedsMax_B_   (iConfig.getParameter<int>("nSeedsMax_B")),
	nSeedsMax_E_   (iConfig.getParameter<int>("nSeedsMax_E"))
{
	produces<TrajectorySeedCollection>();

	mvaEstimator = {
		make_pair( new SeedMvaEstimator(mvaFile_B_0_, mvaScaleMean_B_, mvaScaleStd_B_),
		           new SeedMvaEstimator(mvaFile_E_0_, mvaScaleMean_E_, mvaScaleStd_E_) ),
		make_pair( new SeedMvaEstimator(mvaFile_B_1_, mvaScaleMean_B_, mvaScaleStd_B_),
		           new SeedMvaEstimator(mvaFile_E_1_, mvaScaleMean_E_, mvaScaleStd_E_) ),
		make_pair( new SeedMvaEstimator(mvaFile_B_2_, mvaScaleMean_B_, mvaScaleStd_B_),
		           new SeedMvaEstimator(mvaFile_E_2_, mvaScaleMean_E_, mvaScaleStd_E_) ),
		make_pair( new SeedMvaEstimator(mvaFile_B_3_, mvaScaleMean_B_, mvaScaleStd_B_),
		           new SeedMvaEstimator(mvaFile_E_3_, mvaScaleMean_E_, mvaScaleStd_E_) )
	};
}


MuonHLTSeedMVAClassifier::~MuonHLTSeedMVAClassifier()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called on each new Event  ------------
void MuonHLTSeedMVAClassifier::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	auto result = std::make_unique<TrajectorySeedCollection>();

	edm::ESHandle<TrackerGeometry> trkGeom;
	iSetup.get<TrackerDigiGeometryRecord>().get(trkGeom);

	// edm::Handle<l1t::MuonBxCollection> h_L1Muon;
	// bool hasL1 = iEvent.getByToken( t_L1Muon_, h_L1Muon);

	edm::Handle<l1t::TkMuonCollection> h_L1TkMu;
	bool hasL1TkMu = iEvent.getByToken( t_L1TkMu_, h_L1TkMu);

	edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon;
	bool hasL2 = iEvent.getByToken( t_L2Muon_, h_L2Muon );

	edm::Handle< TrajectorySeedCollection > h_Seed;
	bool hasSeed = iEvent.getByToken( t_Seed_, h_Seed );

	// if( !( hasL1 && hasL1TkMu && hasL2 && hasSeed ) ) {
	if( !( hasL1TkMu && hasL2 && hasSeed ) ) {
		std::cout << "MuonHLTSeedMVAClassifier::produce: !( hasL1 && hasL1TkMu && hasL2 && hasSeed )" << std::endl;
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
				// h_L1Muon,
				h_L1TkMu,
				h_L2Muon,
				0.5
			);

			float softmax = std::exp(mvas.at(3)) / ( std::exp(mvas.at(0)) + std::exp(mvas.at(1)) + std::exp(mvas.at(2)) + std::exp(mvas.at(3)) );

			if(isB)  pairSeedIdxMvaScore_B.push_back( make_pair( i, softmax ) );
			else     pairSeedIdxMvaScore_E.push_back( make_pair( i, softmax ) );
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
				// h_L1Muon,
				h_L1TkMu,
				h_L2Muon,
				0.5
			);

			float softmax = std::exp(mvas.at(3)) / ( std::exp(mvas.at(0)) + std::exp(mvas.at(1)) + std::exp(mvas.at(2)) + std::exp(mvas.at(3)) );

			bool passMva = (
				(  isB && (softmax > mvaCut_B_) ) ||
				( !isB && (softmax > mvaCut_E_) )
			);

			if( passMva )  result->emplace_back( seed );
		}
	}

	iEvent.put(std::move(result));
}

std::vector<float> MuonHLTSeedMVAClassifier::getSeedMva(
	pairSeedMvaEstimator pairMvaEstimator,
	const TrajectorySeed& seed,
	GlobalVector global_p,
	GlobalPoint  global_x,
	// edm::Handle<l1t::MuonBxCollection>& h_L1Muon,
	edm::Handle<l1t::TkMuonCollection>& h_L1TkMu,
	edm::Handle<reco::RecoChargedCandidateCollection>& h_L2Muon,
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
				// h_L1Muon,
				h_L2Muon,
				h_L1TkMu
			);
			v_mva.push_back( (offset + mva) );
		}
		else {
			float mva = pairMvaEstimator.at(ic).second->computeMva(
				seed,
				global_p,
				global_x,
				// h_L1Muon,
				h_L2Muon,
				h_L1TkMu
			);
			v_mva.push_back( (offset + mva) );
		}
	}
	if( v_mva.size() != 4 ) {  // this should never happen
		std::cout << "MuonHLTSeedMVAClassifier::getSeedMva: v_mva.size() != 4" << std::endl;
		return { -99999., -99999., -99999., -99999. };
	}

	return v_mva;
}

void MuonHLTSeedMVAClassifier::beginJob(){}

void MuonHLTSeedMVAClassifier::endJob(){
	for( int i=0; i<4; ++i ) {
		delete mvaEstimator.at(i).first;
		delete mvaEstimator.at(i).second;
	}
}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
/*
void MuonHLTSeedMVAClassifier::beginStream(edm::StreamID){}
*/

// ------------ method called once each stream after processing all runs, lumis and events  ------------
/*
void MuonHLTSeedMVAClassifier::endStream(){}
*/

// ------------ method called when starting to processes a run  ------------
/*
void MuonHLTSeedMVAClassifier::beginRun(edm::Run const&, edm::EventSetup const&){}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void MuonHLTSeedMVAClassifier::endRun(edm::Run const&, edm::EventSetup const&){}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MuonHLTSeedMVAClassifier::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MuonHLTSeedMVAClassifier::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonHLTSeedMVAClassifier::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MuonHLTSeedMVAClassifier);
