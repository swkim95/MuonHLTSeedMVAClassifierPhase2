
import FWCore.ParameterSet.Config as cms
from HLTrigger.MuonHLTSeedMVAClassifierPhase2.mvaScale import *

def customizerFuncForMuonHLTSeeding(
    process, newProcessName = "MYHLT", WPName = "TEST",
    doSort = False, nSeedsMax_B = (-1), nSeedsMax_E = (-1),
    mvaCuts_B = (0), mvaCuts_E = (0) ):

    print("doSort: ", doSort)
    print("nSeedsMax_B: ", nSeedsMax_B)
    print("nSeedsMax_E: ", nSeedsMax_E)
    print("mvaCuts_B: ", mvaCuts_B)
    print("mvaCuts_E: ", mvaCuts_E)

    # -- Seed MVA Classifiers
    process.hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered = cms.EDProducer("MuonHLTSeedMVAClassifierPhase2",
        src    = cms.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeeds", "", newProcessName),
        # src    = cms.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeeds", "", "HLT"),
        L1TkMu = cms.InputTag("L1TkMuons", "", newProcessName),
        # L1TkMu = cms.InputTag("hltL1TkMuons", "", newProcessName),
        # L1TkMu = cms.InputTag("L1TkMuons", "", "HLT"),
        L2Muon = cms.InputTag("hltL2MuonFromL1TkMuonCandidates", "", newProcessName),
        # L2Muon = cms.InputTag("hltL2MuonFromL1TkMuonCandidates", "", "HLT"),


        # mvaFile_B_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/DY_PU200_Barrel_Binary_NThltIter2FromL1_0.xml"),
        # mvaFile_B_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/CMSSW_13_0_9_Barrel_Binary_NThltIter2FromL1_0.xml"),
        # mvaFile_B_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/xgb_CMSSW_13_0_9_allData_updated_Barrel_NThltIter2FromL1_0.xml"),
        mvaFile_B_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/xgb_CMSSW_13_0_9_allData_Try_05_Barrel_NThltIter2FromL1_0.xml"),
        # mvaFile_E_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/DY_PU200_Endcap_Binary_NThltIter2FromL1_0.xml"),
        # mvaFile_E_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/CMSSW_13_0_9_Endcap_Binary_NThltIter2FromL1_0.xml"),
        # mvaFile_E_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/xgb_CMSSW_13_0_9_allData_updated_Endcap_NThltIter2FromL1_0.xml"),
        mvaFile_E_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/xgb_CMSSW_13_0_9_allData_Try_05_Endcap_NThltIter2FromL1_0.xml"),

        mvaScaleMean_B = cms.vdouble(PU200_Barrel_NThltIter2FromL1_ScaleMean),
        mvaScaleStd_B  = cms.vdouble(PU200_Barrel_NThltIter2FromL1_ScaleStd),
        mvaScaleMean_E = cms.vdouble(PU200_Endcap_NThltIter2FromL1_ScaleMean),
        mvaScaleStd_E  = cms.vdouble(PU200_Endcap_NThltIter2FromL1_ScaleStd),

        doSort = cms.bool(doSort),
        nSeedsMax_B = cms.int32(nSeedsMax_B[0]),
        nSeedsMax_E = cms.int32(nSeedsMax_E[0]),

        mvaCut_B = cms.double(mvaCuts_B[0]),
        mvaCut_E = cms.double(mvaCuts_E[0])
    )

    # -- Track Candidates
    process.hltIter2Phase2L3FromL1TkMuonCkfTrackCandidates.src = cms.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered", "", newProcessName)
    # process.hltIter2Phase2L3FromL1TkMuonCkfTrackCandidates.src = cms.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered", "", "HLT")

    # -- Sequences
    # process.HLTPhase2L3IOFromL1TkMuonTkCandidateSequence = cms.Sequence(
    #     process.hltPhase2L3MuonPixelTracksFilter+
    #     process.hltPhase2L3MuonPixelTracksFitter+
    #     process.hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions+
    #     process.hltPhase2L3FromL1TkMuonPixelLayerQuadruplets+
    #     process.hltPhase2L3FromL1TkMuonPixelTracksHitDoublets+
    #     process.hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets+
    #     process.hltPhase2L3FromL1TkMuonPixelTracks+
    #     process.hltPhase2L3FromL1TkMuonPixelVertices+
    #     process.hltPhase2L3FromL1TkMuonTrimmedPixelVertices+
    #     process.hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks+
    #     process.hltIter0Phase2L3FromL1TkMuonCkfTrackCandidates+
    #     process.hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks+
    #     process.hltIter0Phase2L3FromL1TkMuonTrackCutClassifier+
    #     process.hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity+
    #     process.hltIter2Phase2L3FromL1TkMuonClustersRefRemoval+
    #     process.hltIter2Phase2L3FromL1TkMuonMaskedMeasurementTrackerEvent+
    #     process.hltIter2Phase2L3FromL1TkMuonPixelLayerTriplets+
    #     process.hltIter2Phase2L3FromL1TkMuonPixelClusterCheck+
    #     process.hltIter2Phase2L3FromL1TkMuonPixelHitDoublets+
    #     process.hltIter2Phase2L3FromL1TkMuonPixelHitTriplets+
    #     process.hltIter2Phase2L3FromL1TkMuonPixelSeeds+
    #     process.hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered+  # HERE
    #     process.hltIter2Phase2L3FromL1TkMuonCkfTrackCandidates+
    #     process.hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracks+
    #     process.hltIter2Phase2L3FromL1TkMuonTrackCutClassifier+
    #     process.hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity+
    #     process.hltIter2Phase2L3FromL1TkMuonMerged
    # )

    # -- DQMOutput and Ntupler
    # if hasattr(process, "DQMOutput"):
    #     del process.DQMOutput

    if hasattr(process, "dqmOutput"):
        process.dqmOutput.fileName = cms.untracked.string("DQMIO_%s.root" % WPName)

    if hasattr(process, "TFileService"):
        process.TFileService.fileName = cms.string("ntuple_%s.root" % WPName)

    if hasattr(process, "writeDataset"):
        process.writeDataset.fileName = cms.untracked.string("edmOutput_%s.root" % WPName)

    if hasattr(process, "ntupler"):
        process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered", "", newProcessName)
        # process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered", "", "HLT")

    if hasattr(process, "seedNtupler"):
        process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered", "", newProcessName)
        # process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered", "", "HLT")

    return process
