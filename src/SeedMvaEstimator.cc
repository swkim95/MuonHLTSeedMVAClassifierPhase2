#include "HLTrigger/MuonHLTSeedMVAClassifier/interface/SeedMvaEstimator.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CommonTools/MVAUtils/interface/GBRForestTools.h"

#include "TMath.h"

using namespace std;

SeedMvaEstimator::SeedMvaEstimator(const edm::FileInPath& weightsfile, std::vector<double> scale_mean, std::vector<double> scale_std) {
  gbrForest_  = createGBRForest(weightsfile);
  scale_mean_ = scale_mean;
  scale_std_  = scale_std;
}

SeedMvaEstimator::~SeedMvaEstimator() {}

namespace {
  enum inputIndexes {
    kTsosErr0,         // 0
    // kTsosErr1,         // 1
    kTsosErr2,         // 1  2
    // kTsosErr3,         // 3
    // kTsosErr4,         // 4
    kTsosErr5,         // 2  5
    // kTsosErr6,         // 6
    // kTsosErr7,         // 7
    // kTsosErr8,         // 8
    kTsosErr9,         // 3  9
    // kTsosErr10,        // 10
    // kTsosErr11,        // 11
    // kTsosErr12,        // 12
    // kTsosErr13,        // 13
    kTsosErr14,        // 4  14
    kTsosDxdz,         // 5  15
    kTsosDydz,         // 6  16
    kTsosQbp,          // 7  17
    kTsosCharge,       // 8  18
    // kDRdRL1SeedP,      // 19
    // kDPhidRL1SeedP,    // 20
    // kDRdPhiL1SeedX,    // 21
    // kDPhidPhiL1SeedX,  // 22
    kDRdRL2SeedP,      // 9  23
    kDPhidRL2SeedP,    // 10  24
    // kDRdPhiL2SeedX,    // 25
    // kDPhidPhiL2SeedX,  // 26
    kDRL1TkMu,         // 11  27
    kDPhiL1TkMu,       // 12  28
    kHasL2,            // 13
    kLast              // 14  29
  };
}  // namespace

/*
void SeedMvaEstimator::getL1MuonVariables(
  const TrajectorySeed& seed,
  GlobalVector global_p,
  GlobalPoint  global_x,
  edm::Handle<l1t::MuonBxCollection> h_L1Muon,
  float& dRdRL1SeedP,
  float& dPhidRL1SeedP,
  float& dRdPhiL1SeedX,
  float& dPhidPhiL1SeedX ) const {

  for(int ibx = h_L1Muon->getFirstBX(); ibx<=h_L1Muon->getLastBX(); ++ibx)
  {
    if(ibx != 0) continue; // -- only take when ibx == 0 -- //
    for(auto it=h_L1Muon->begin(ibx); it!=h_L1Muon->end(ibx); it++)
    {
      l1t::MuonRef ref_L1Mu(h_L1Muon, distance(h_L1Muon->begin(h_L1Muon->getFirstBX()), it) );

      // FIXME: 7 should be configurable
      if(ref_L1Mu->hwQual() < 7)
        continue;

      float dR_L1SeedP_AtVtx   = reco::deltaR( ref_L1Mu->etaAtVtx(), ref_L1Mu->phiAtVtx(), global_p.eta(), global_p.phi());
      float dPhi_L1SeedP_AtVtx = reco::deltaPhi( ref_L1Mu->phiAtVtx(), global_p.phi());
      float dR_L1SeedX_AtVtx   = reco::deltaR( ref_L1Mu->etaAtVtx(), ref_L1Mu->phiAtVtx(), global_x.eta(), global_x.phi());
      float dPhi_L1SeedX_AtVtx = reco::deltaPhi( ref_L1Mu->phiAtVtx(), global_x.phi());

      if( dR_L1SeedP_AtVtx < dRdRL1SeedP ) {
        dRdRL1SeedP = dR_L1SeedP_AtVtx;
        dPhidRL1SeedP = dPhi_L1SeedP_AtVtx;
      }
      if( fabs(dPhi_L1SeedX_AtVtx) < fabs(dPhidPhiL1SeedX) ) {
        dRdPhiL1SeedX = dR_L1SeedX_AtVtx;
        dPhidPhiL1SeedX = dPhi_L1SeedX_AtVtx;
      }
    }
  }
}
*/

void SeedMvaEstimator::getL2MuonVariables(
  const TrajectorySeed& seed,
  GlobalVector global_p,
  GlobalPoint  global_x,
  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon,
  float& dRdRL2SeedP,
  float& dPhidRL2SeedP,
  float& dRdPhiL2SeedX,
  float& dPhidPhiL2SeedX ) const {

  for( unsigned int i_L2=0; i_L2<h_L2Muon->size(); i_L2++)
  {
    reco::RecoChargedCandidateRef ref_L2Mu(h_L2Muon, i_L2);

    float dR_L2SeedP   = reco::deltaR( *ref_L2Mu, global_p);
    float dPhi_L2SeedP = reco::deltaPhi( ref_L2Mu->phi(), global_p.phi());
    float dR_L2SeedX   = reco::deltaR( *ref_L2Mu, global_x);
    float dPhi_L2SeedX = reco::deltaPhi( ref_L2Mu->phi(), global_x.phi());

    if( dR_L2SeedP < dRdRL2SeedP ) {
      dRdRL2SeedP = dR_L2SeedP;
      dPhidRL2SeedP = dPhi_L2SeedP;
    }
    if( fabs(dPhi_L2SeedX) < fabs(dPhidPhiL2SeedX) ) {
      dRdPhiL2SeedX = dR_L2SeedX;
      dPhidPhiL2SeedX = dPhi_L2SeedX;
    }
  }
}

void SeedMvaEstimator::getL1TTVariables(
  const TrajectorySeed& seed,
  GlobalVector global_p,
  GlobalPoint  global_x,
  edm::Handle<l1t::TkMuonCollection> h_L1TkMu,
  float& DRL1TkMu,
  float& DPhiL1TkMu ) const {

  for(auto L1TkMu=h_L1TkMu->begin(); L1TkMu!=h_L1TkMu->end(); ++L1TkMu)
  {
    auto TkRef = L1TkMu->trkPtr();
    float DRL1TkMu_tmp   = reco::deltaR( *TkRef, global_p);
    float DPhiL1TkMu_tmp = reco::deltaPhi( TkRef->phi(), global_p.phi());
    if( DRL1TkMu_tmp < DRL1TkMu ) {
      DRL1TkMu   = DRL1TkMu_tmp;
      DPhiL1TkMu = DPhiL1TkMu_tmp;
    }
  }
}

float SeedMvaEstimator::computeMva( const TrajectorySeed& seed,
  GlobalVector global_p,
  GlobalPoint  global_x,
  //  edm::Handle<l1t::MuonBxCollection> h_L1Muon,
  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon,
  edm::Handle<l1t::TkMuonCollection> h_L1TkMu
) const {

  float var[kLast]{};

  var[kTsosErr0]   = seed.startingState().error(0);
  // var[kTsosErr1]   = seed.startingState().error(1);
  var[kTsosErr2]   = seed.startingState().error(2);
  // var[kTsosErr3]   = seed.startingState().error(3);
  // var[kTsosErr4]   = seed.startingState().error(4);
  var[kTsosErr5]   = seed.startingState().error(5);
  // var[kTsosErr6]   = seed.startingState().error(6);
  // var[kTsosErr7]   = seed.startingState().error(7);
  // var[kTsosErr8]   = seed.startingState().error(8);
  var[kTsosErr9]   = seed.startingState().error(9);
  // var[kTsosErr10]  = seed.startingState().error(10);
  // var[kTsosErr11]  = seed.startingState().error(11);
  // var[kTsosErr12]  = seed.startingState().error(12);
  // var[kTsosErr13]  = seed.startingState().error(13);
  var[kTsosErr14]  = seed.startingState().error(14);
  var[kTsosDxdz]   = seed.startingState().parameters().dxdz();
  var[kTsosDydz]   = seed.startingState().parameters().dydz();
  var[kTsosQbp]    = seed.startingState().parameters().qbp();
  var[kTsosCharge] = seed.startingState().parameters().charge();

  // FIXME: should be configurable
  float initDRdPhi = 99999.;

  // float dRdRL1SeedP = initDRdPhi;
  // float dPhidRL1SeedP = initDRdPhi;
  // float dRdPhiL1SeedX = initDRdPhi;
  // float dPhidPhiL1SeedX = initDRdPhi;
  // getL1MuonVariables( seed, global_p, global_x, h_L1Muon, dRdRL1SeedP, dPhidRL1SeedP, dRdPhiL1SeedX, dPhidPhiL1SeedX );

  float dRdRL2SeedP = initDRdPhi;
  float dPhidRL2SeedP = initDRdPhi;
  float dRdPhiL2SeedX = initDRdPhi;
  float dPhidPhiL2SeedX = initDRdPhi;
  getL2MuonVariables( seed, global_p, global_x, h_L2Muon, dRdRL2SeedP, dPhidRL2SeedP, dRdPhiL2SeedX, dPhidPhiL2SeedX );

  float hasL2 = float(dRdRL2SeedP < initDRdPhi);

  float DRL1TkMu = initDRdPhi;
  float DPhiL1TkMu = initDRdPhi;
  getL1TTVariables( seed, global_p, global_x, h_L1TkMu, DRL1TkMu, DPhiL1TkMu );

  // var[kDRdRL1SeedP]     = dRdRL1SeedP;
  // var[kDPhidRL1SeedP]   = dPhidRL1SeedP;
  // var[kDRdPhiL1SeedX]   = dRdPhiL1SeedX;
  // var[kDPhidPhiL1SeedX] = dPhidPhiL1SeedX;
  var[kDRdRL2SeedP]     = dRdRL2SeedP;
  var[kDPhidRL2SeedP]   = dPhidRL2SeedP;
  // var[kDRdPhiL2SeedX]   = dRdPhiL2SeedX;
  // var[kDPhidPhiL2SeedX] = dPhidPhiL2SeedX;
  var[kDRL1TkMu]        = DRL1TkMu;
  var[kDPhiL1TkMu]      = DPhiL1TkMu;
  var[kHasL2]           = hasL2;

  for(int iv=0; iv<kLast; ++iv) {
    var[iv] = (var[iv] - scale_mean_.at(iv)) / scale_std_.at(iv);
  }

  return gbrForest_->GetResponse( var );
}

/////////////////////////// 
//Phase2 SeedMvaEstimator//
///////////////////////////

// Phase2 enum namespace
namespace Phase2{
  enum inputIndexesPhase2 {
    kTsosErr0,          // 0
    kTsosErr2,          // 1
    kTsosErr5,          // 2
    kTsosErr9,          // 3
    kTsosErr14,         // 4
    kTsosDxdz,          // 5
    kTsosDydz,          // 6
    kTsosQbp,           // 7
    kTsosCharge,        // 8
    kDRdRL1TkMuSeedP,   // 9
    kDRdPhiL1TkMuSeedP, // 10
    kExpd2HitL1Tk1,     // 11
    kExpd2HitL1Tk2,     // 12
    kExpd2HitL1Tk3,     // 13
    kLast               // 14
  };
}

SeedMvaEstimatorPhase2::SeedMvaEstimatorPhase2(const edm::FileInPath& weightsfile, std::vector<double> scale_mean, std::vector<double> scale_std) {
  gbrForest_  = createGBRForest(weightsfile);
  scale_mean_ = scale_mean;
  scale_std_  = scale_std;
}

SeedMvaEstimatorPhase2::~SeedMvaEstimatorPhase2() {}

void SeedMvaEstimatorPhase2::getL1TTVariables( const TrajectorySeed& seed,
  GlobalVector global_p,
  GlobalPoint  global_x,
  edm::Handle<l1t::TkMuonCollection> h_L1TkMu,
  float& DRL1TkMu,
  float& DPhiL1TkMu ) const {

  for(auto L1TkMu=h_L1TkMu->begin(); L1TkMu!=h_L1TkMu->end(); ++L1TkMu)
  {
    auto TkRef = L1TkMu->trkPtr();
    float DRL1TkMu_tmp   = reco::deltaR( *TkRef, global_p);
    float DPhiL1TkMu_tmp = reco::deltaPhi( TkRef->phi(), global_p.phi());
    if( DRL1TkMu_tmp < DRL1TkMu ) {
      DRL1TkMu   = DRL1TkMu_tmp;
      DPhiL1TkMu = DPhiL1TkMu_tmp;
    }
  }
}

vector< LayerTSOS > SeedMvaEstimatorPhase2::getTsosOnPixels(
  TTTrack<Ref_Phase2TrackerDigi_> l1tk,
  edm::ESHandle<MagneticField>& magfieldH,
  const Propagator& propagatorAlong,
  GeometricSearchTracker* geomTracker
) const {
  vector< LayerTSOS > v_tsos = {};

  std::vector<BarrelDetLayer const*>  const&  bpix = geomTracker->pixelBarrelLayers();
  std::vector<ForwardDetLayer const*> const& nfpix = geomTracker->negPixelForwardLayers();
  std::vector<ForwardDetLayer const*> const& pfpix = geomTracker->posPixelForwardLayers();

  int chargeTk = l1tk.rInv() > 0. ? 1 : -1;
  GlobalPoint  gpos = l1tk.POCA();
  GlobalVector gmom = l1tk.momentum();

  FreeTrajectoryState fts = FreeTrajectoryState( gpos, gmom, chargeTk, magfieldH.product() );

  for(auto it = bpix.begin(); it != bpix.end(); ++it) {
    TrajectoryStateOnSurface tsos = propagatorAlong.propagate(fts, (**it).specificSurface());
    if( !tsos.isValid() )  continue;

    auto z0 = std::abs(tsos.globalPosition().z() - (**it).specificSurface().position().z());
    auto deltaZ = 0.5f * (**it).surface().bounds().length();
    deltaZ += 0.5f * (**it).surface().bounds().thickness() * std::abs(tsos.globalDirection().z()) / tsos.globalDirection().perp();
    bool is_compatible  = (z0 < deltaZ);

    if( is_compatible ) {
      v_tsos.push_back( make_pair( (const DetLayer*)(*it), tsos) );
    }
  }
  for(auto it = nfpix.begin(); it != nfpix.end(); ++it) {
    TrajectoryStateOnSurface tsos = propagatorAlong.propagate(fts, (**it).specificSurface());
    if( !tsos.isValid() )  continue;

    auto r2 = tsos.localPosition().perp2();
    float deltaR = 0.5f * (**it).surface().bounds().thickness() * tsos.localDirection().perp() / std::abs(tsos.localDirection().z());
    auto ri2 = std::max((**it).specificSurface().innerRadius() - deltaR, 0.f);
    ri2 *= ri2;
    auto ro2 = (**it).specificSurface().outerRadius() + deltaR;
    ro2 *= ro2;
    bool is_compatible  = ((r2 > ri2) && (r2 < ro2));

    if( is_compatible ) {
      v_tsos.push_back( make_pair( (const DetLayer*)(*it), tsos) );
    }
  }
  for(auto it = pfpix.begin(); it != pfpix.end(); ++it) {
    TrajectoryStateOnSurface tsos = propagatorAlong.propagate(fts, (**it).specificSurface());
    if( !tsos.isValid() )  continue;

    auto r2 = tsos.localPosition().perp2();
    float deltaR = 0.5f * (**it).surface().bounds().thickness() * tsos.localDirection().perp() / std::abs(tsos.localDirection().z());
    auto ri2 = std::max((**it).specificSurface().innerRadius() - deltaR, 0.f);
    ri2 *= ri2;
    auto ro2 = (**it).specificSurface().outerRadius() + deltaR;
    ro2 *= ro2;
    bool is_compatible  = ((r2 > ri2) && (r2 < ro2));

    if( is_compatible ) {
      v_tsos.push_back( make_pair( (const DetLayer*)(*it), tsos) );
    }
  }

  return v_tsos;
}

vector< pair<LayerHit, LayerTSOS> > SeedMvaEstimatorPhase2::getHitTsosPairs( const TrajectorySeed& seed,
  edm::Handle<l1t::TkMuonCollection> L1TkMuonHandle,
  edm::ESHandle<MagneticField>& magfieldH,
  const Propagator& propagatorAlong,
  GeometricSearchTracker* geomTracker
) const {
  vector< pair<LayerHit, LayerTSOS> > out = {};

  // FIXME: this is random choice
  float av_dr_min = 20.;

  // -- loop on L1TkMu
  for(auto L1TkMu=L1TkMuonHandle->begin(); L1TkMu!=L1TkMuonHandle->end(); ++L1TkMu) {

    const vector< LayerTSOS > v_tsos = getTsosOnPixels(
      *L1TkMu->trkPtr(),
      magfieldH,
      propagatorAlong,
      geomTracker
    );

    // -- loop on recHits
    vector<int> v_tsos_skip( v_tsos.size(), 0 );
    vector< pair<LayerHit, LayerTSOS> > hitTsosPair = {};
    int ihit = 0;
    for( auto hit = seed.recHits().first; hit!=seed.recHits().second; ++hit ) {
      // -- look for closest tsos by absolute distance
      // FIXME: this is random choice
      int the_tsos = -99999;
      float dr_min = 20.;
      for( auto i=0U; i<v_tsos.size(); ++i ) {
        if( v_tsos_skip.at(i) )  continue;
        float dr = ( v_tsos.at(i).second.globalPosition() - hit->globalPosition() ).mag();
        if( dr < dr_min ) {
          dr_min = dr;
          the_tsos = i;
          v_tsos_skip.at(i) = 1;
        }
      }

      if( the_tsos > -1 ) {
        const DetLayer* thelayer =  geomTracker->idToLayer( hit->geographicalId() );
        hitTsosPair.push_back( make_pair( make_pair( thelayer, &*hit), v_tsos.at(the_tsos) ) );
      }

      ihit++;
    } // loop on recHits

    // -- find tsos for all recHits?
    // FIXME: this is random choice
    if( (int)hitTsosPair.size() == ihit ) {
      float av_dr = 0.;
      for( auto it=hitTsosPair.begin(); it!=hitTsosPair.end(); ++it ) {
        auto hit  = it->first.second;
        auto tsos = it->second.second;
        av_dr += ( hit->globalPosition() - tsos.globalPosition() ).mag();
      }
      av_dr = av_dr > 0 ? av_dr / (float)hitTsosPair.size() : 1.e6;

      if( av_dr < av_dr_min ) {
        av_dr_min = av_dr;
        out = hitTsosPair;
      }
    }

  }  // loop on L1TkMu

  return out;
}

void SeedMvaEstimatorPhase2::getHitL1TkVatiables( const TrajectorySeed& seed,
  edm::Handle<l1t::TkMuonCollection> L1TkMuonHandle,
  edm::ESHandle<MagneticField>& magfieldH,
  const Propagator& propagatorAlong,
  GeometricSearchTracker* geomTracker,
  float& expd2HitL1Tk1,
  float& expd2HitL1Tk2,
  float& expd2HitL1Tk3 ) const {

    vector< pair<LayerHit, LayerTSOS> > hitTsosPair = getHitTsosPairs(
      seed,    
      L1TkMuonHandle,
      magfieldH,        
      propagatorAlong,
      geomTracker
    );

    bool found = (hitTsosPair.size() > 0);

    if(found) {
      float l1x1 = 99999.;
      float l1y1 = 99999.;
      float l1z1 = 99999.;
      float hitx1 = 99999.;
      float hity1 = 99999.;
      float hitz1 = 99999.;

      float l1x2 = 99999.;
      float l1y2 = 99999.;
      float l1z2 = 99999.;
      float hitx2 = 99999.;
      float hity2 = 99999.;
      float hitz2 = 99999.;

      float l1x3 = 99999.;
      float l1y3 = 99999.;
      float l1z3 = 99999.;
      float hitx3 = 99999.;
      float hity3 = 99999.;
      float hitz3 = 99999.;

      if (hitTsosPair.size() > 1) {
        auto hit1 = hitTsosPair.at(0).first.second;
        auto tsos1 = hitTsosPair.at(0).second.second;

        l1x1 = tsos1.globalPosition().x();
        l1y1 = tsos1.globalPosition().y();
        l1z1 = tsos1.globalPosition().z();
        hitx1 = hit1->globalPosition().x();
        hity1 = hit1->globalPosition().y();
        hitz1 = hit1->globalPosition().z();

        auto hit2 = hitTsosPair.at(1).first.second;
        auto tsos2 = hitTsosPair.at(1).second.second;

        l1x2 = tsos2.globalPosition().x();
        l1y2 = tsos2.globalPosition().y();
        l1z2 = tsos2.globalPosition().z();
        hitx2 = hit2->globalPosition().x();
        hity2 = hit2->globalPosition().y();
        hitz2 = hit2->globalPosition().z();
      }
      if (hitTsosPair.size() > 2) {
        auto hit3 = hitTsosPair.at(2).first.second;
        auto tsos3 = hitTsosPair.at(2).second.second;

        l1x3 = tsos3.globalPosition().x();
        l1y3 = tsos3.globalPosition().y();
        l1z3 = tsos3.globalPosition().z();
        hitx3 = hit3->globalPosition().x();
        hity3 = hit3->globalPosition().y();
        hitz3 = hit3->globalPosition().z();
      }

      // FIXME: (expd~ == NaN) was used in python training when (hit == 99999)
      // >> here we used 99999. , is it proper?
      if (hitx1 != 99999.){
        expd2HitL1Tk1 = TMath::Exp(TMath::Power((l1x1-hitx1),2) + TMath::Power((l1y1-hity1),2) + TMath::Power((l1z1-hitz1),2) );
      }
      if (hitx2 != 99999.){
        expd2HitL1Tk2 = TMath::Exp(TMath::Power((l1x2-hitx2),2) + TMath::Power((l1y2-hity2),2) + TMath::Power((l1z2-hitz2),2) );
      }
      if (hitx3 != 99999.){
        expd2HitL1Tk3 = TMath::Exp(TMath::Power((l1x3-hitx3),2) + TMath::Power((l1y3-hity3),2) + TMath::Power((l1z3-hitz3),2) );
      }
    }
  }

float SeedMvaEstimatorPhase2::computeMva( const TrajectorySeed& seed,
  GlobalVector global_p,
  GlobalPoint  global_x,
  edm::Handle<l1t::TkMuonCollection> h_L1TkMu,
  edm::ESHandle<MagneticField>& magfieldH,
  const Propagator& propagatorAlong,
  GeometricSearchTracker* geomTracker
) const {

  float Phase2var[Phase2::kLast]{};

  Phase2var[Phase2::kTsosErr0]   = seed.startingState().error(0);
  Phase2var[Phase2::kTsosErr2]   = seed.startingState().error(2);
  Phase2var[Phase2::kTsosErr5]   = seed.startingState().error(5);
  Phase2var[Phase2::kTsosErr9]   = seed.startingState().error(9);
  Phase2var[Phase2::kTsosErr14]  = seed.startingState().error(14);
  Phase2var[Phase2::kTsosDxdz]   = seed.startingState().parameters().dxdz();
  Phase2var[Phase2::kTsosDydz]   = seed.startingState().parameters().dydz();
  Phase2var[Phase2::kTsosQbp]    = seed.startingState().parameters().qbp();
  Phase2var[Phase2::kTsosCharge] = seed.startingState().parameters().charge();

  // FIXME: should be configurable
  float initDRdPhi = 99999.;

  float DRL1TkMuSeedP = initDRdPhi;
  float DPhiL1TkMuSeedP = initDRdPhi;
  getL1TTVariables( seed, global_p, global_x, h_L1TkMu, DRL1TkMuSeedP, DPhiL1TkMuSeedP );
  Phase2var[Phase2::kDRdRL1TkMuSeedP]        = DRL1TkMuSeedP;
  Phase2var[Phase2::kDRdPhiL1TkMuSeedP]      = DPhiL1TkMuSeedP;
  
  float expd2HitL1Tk1 = initDRdPhi;
  float expd2HitL1Tk2 = initDRdPhi;
  float expd2HitL1Tk3 = initDRdPhi;
  getHitL1TkVatiables( seed, h_L1TkMu, magfieldH, propagatorAlong, geomTracker, expd2HitL1Tk1, expd2HitL1Tk2, expd2HitL1Tk3);
  Phase2var[Phase2::kExpd2HitL1Tk1]   = expd2HitL1Tk1;
  Phase2var[Phase2::kExpd2HitL1Tk2]   = expd2HitL1Tk2;
  Phase2var[Phase2::kExpd2HitL1Tk3]   = expd2HitL1Tk3;

  for(int iv=0; iv<Phase2::kLast; ++iv) {
    Phase2var[iv] = (Phase2var[iv] - scale_mean_.at(iv)) / scale_std_.at(iv);
  }

  return gbrForest_->GetResponse( Phase2var );
}
