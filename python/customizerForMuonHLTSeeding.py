
import FWCore.ParameterSet.Config as cms
from HLTrigger.MuonHLTSeedMVAClassifierPhase2.mvaScale import *

def customizerFuncForMuonHLTSeeding(
    process, newProcessName = "MYHLT", WPName = "TEST",
    doSort = False, nSeedsMax_B = (-1), nSeedsMax_E = (-1),
    mvaCuts_B = (0), mvaCuts_E = (0), etaEdge = (1.2), baseScore = (0.5) ):

    print("doSort: ", doSort)
    print("nSeedsMax_B: ", nSeedsMax_B)
    print("nSeedsMax_E: ", nSeedsMax_E)
    print("mvaCuts_B: ", mvaCuts_B)
    print("mvaCuts_E: ", mvaCuts_E)
    print("etaEdge: ", etaEdge)
    print("baseScore:", baseScore)

    # -- Seed MVA Classifiers
    process.hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered = cms.EDProducer("MuonHLTSeedMVAClassifierPhase2",
        src    = cms.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeeds", "", newProcessName),
        L1TkMu = cms.InputTag("L1TkMuons", "", newProcessName),
        L2Muon = cms.InputTag("hltL2MuonFromL1TkMuonCandidates", "", newProcessName),


        mvaFile_B_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/xgb_Barrel_NThltIter2FromL1_0.xml"),
        mvaFile_E_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/xgb_Endcap_NThltIter2FromL1_0.xml"),

        mvaScaleMean_B = cms.vdouble(PU200_Barrel_NThltIter2FromL1_ScaleMean),
        mvaScaleStd_B  = cms.vdouble(PU200_Barrel_NThltIter2FromL1_ScaleStd),
        mvaScaleMean_E = cms.vdouble(PU200_Endcap_NThltIter2FromL1_ScaleMean),
        mvaScaleStd_E  = cms.vdouble(PU200_Endcap_NThltIter2FromL1_ScaleStd),

        doSort = cms.bool(doSort),
        nSeedsMax_B = cms.int32(nSeedsMax_B[0]),
        nSeedsMax_E = cms.int32(nSeedsMax_E[0]),

        etaEdge  = cms.double(etaEdge[0]),
        mvaCut_B = cms.double(mvaCuts_B[0]),
        mvaCut_E = cms.double(mvaCuts_E[0]),

        baseScore = cms.double(baseScore[0])
    )

    # -- Track Candidates
    process.hltIter2Phase2L3FromL1TkMuonCkfTrackCandidates.src = cms.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered", "", newProcessName)

    if hasattr(process, "dqmOutput"):
        process.dqmOutput.fileName = cms.untracked.string("DQMIO_%s.root" % WPName)

    if hasattr(process, "TFileService"):
        process.TFileService.fileName = cms.string("ntuple_%s.root" % WPName)

    if hasattr(process, "writeDataset"):
        process.writeDataset.fileName = cms.untracked.string("edmOutput_%s.root" % WPName)

    if hasattr(process, "ntupler"):
        process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered", "", newProcessName)

    if hasattr(process, "seedNtupler"):
        process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered", "", newProcessName)

    return process
