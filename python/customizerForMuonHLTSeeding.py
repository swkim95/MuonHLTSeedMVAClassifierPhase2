
import FWCore.ParameterSet.Config as cms
import HLTrigger.MuonHLTSeedMVAClassifier.mvaScale as _mvaScale

def customizerFuncForMuonHLTSeeding(
    process, newProcessName = "MYHLT", version = "Run3v0", WPName = "TEST",
    doSort = False, nSeedsMax_B = (-1, -1, -1, -1, -1, -1, -1), nSeedsMax_E = (-1, -1, -1, -1, -1, -1, -1),
    mvaCuts_B = (0, 0, 0, 0, 0, 0, 0), mvaCuts_E = (0, 0, 0, 0, 0, 0, 0) ):

    print "\nCustomizing Seed MVA Classifier:"
    print "\tdoSort:      ", doSort
    print "\tnSeedsMax_B: ", nSeedsMax_B
    print "\tnSeedsMax_E: ", nSeedsMax_E
    print "\tmvaCuts_B:   ", mvaCuts_B
    print "\tmvaCuts_E:   ", mvaCuts_E

    # -- Seed MVA Classifiers
    process.hltIter2IterL3MuonPixelSeedsFiltered = cms.EDProducer("MuonHLTSeedMVAClassifier",
        src    = cms.InputTag("hltIter2IterL3MuonPixelSeeds", "", newProcessName),
        L1Muon = cms.InputTag("hltGtStage2Digis", "Muon", newProcessName),
        # L1Muon = cms.InputTag("hltGtStage2Digis", "Muon", "HLT"),
        # L1Muon = cms.InputTag("simGmtStage2Digis","",newProcessName),
        L2Muon = cms.InputTag("hltL2MuonCandidates", "", newProcessName),

        mvaFile_B_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter2_0.xml" % version),
        mvaFile_B_1 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter2_1.xml" % version),
        mvaFile_B_2 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter2_2.xml" % version),
        mvaFile_B_3 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter2_3.xml" % version),
        mvaFile_E_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter2_0.xml" % version),
        mvaFile_E_1 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter2_1.xml" % version),
        mvaFile_E_2 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter2_2.xml" % version),
        mvaFile_E_3 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter2_3.xml" % version),

        mvaScaleMean_B = cms.vdouble( getattr(_mvaScale, "%sBarrel_hltIter2_ScaleMean" % version) ),
        mvaScaleStd_B  = cms.vdouble( getattr(_mvaScale, "%sBarrel_hltIter2_ScaleStd" % version) ),
        mvaScaleMean_E = cms.vdouble( getattr(_mvaScale, "%sEndcap_hltIter2_ScaleMean" % version) ),
        mvaScaleStd_E  = cms.vdouble( getattr(_mvaScale, "%sEndcap_hltIter2_ScaleStd" % version) ),

        doSort = cms.bool(doSort),
        nSeedsMax_B = cms.int32(nSeedsMax_B[2]),
        nSeedsMax_E = cms.int32(nSeedsMax_E[2]),

        mvaCut_B = cms.double(mvaCuts_B[2]),
        mvaCut_E = cms.double(mvaCuts_E[2])
    )

    process.hltIter3IterL3MuonPixelSeedsFiltered = cms.EDProducer("MuonHLTSeedMVAClassifier",

        # Reject all seeds
        rejectAll = cms.bool(True),

        src    = cms.InputTag("hltIter3IterL3MuonPixelSeeds", "", newProcessName),
        L1Muon = cms.InputTag("hltGtStage2Digis", "Muon", newProcessName),
        # L1Muon = cms.InputTag("hltGtStage2Digis", "Muon", "HLT"),
        # L1Muon = cms.InputTag("simGmtStage2Digis","",newProcessName),
        L2Muon = cms.InputTag("hltL2MuonCandidates", "", newProcessName),

        mvaFile_B_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter3_0.xml" % version),
        mvaFile_B_1 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter3_1.xml" % version),
        mvaFile_B_2 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter3_2.xml" % version),
        mvaFile_B_3 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter3_3.xml" % version),
        mvaFile_E_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter3_0.xml" % version),
        mvaFile_E_1 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter3_1.xml" % version),
        mvaFile_E_2 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter3_2.xml" % version),
        mvaFile_E_3 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter3_3.xml" % version),

        mvaScaleMean_B = cms.vdouble( getattr(_mvaScale, "%sBarrel_hltIter3_ScaleMean" % version) ),
        mvaScaleStd_B  = cms.vdouble( getattr(_mvaScale, "%sBarrel_hltIter3_ScaleStd" % version) ),
        mvaScaleMean_E = cms.vdouble( getattr(_mvaScale, "%sEndcap_hltIter3_ScaleMean" % version) ),
        mvaScaleStd_E  = cms.vdouble( getattr(_mvaScale, "%sEndcap_hltIter3_ScaleStd" % version) ),

        doSort = cms.bool(doSort),
        nSeedsMax_B = cms.int32(nSeedsMax_B[3]),
        nSeedsMax_E = cms.int32(nSeedsMax_E[3]),

        mvaCut_B = cms.double(mvaCuts_B[3]),
        mvaCut_E = cms.double(mvaCuts_E[3])
    )

    process.hltIter2IterL3FromL1MuonPixelSeedsFiltered = cms.EDProducer("MuonHLTSeedMVAClassifier",
        src    = cms.InputTag("hltIter2IterL3FromL1MuonPixelSeeds", "", newProcessName),
        L1Muon = cms.InputTag("hltGtStage2Digis", "Muon", newProcessName),
        # L1Muon = cms.InputTag("hltGtStage2Digis", "Muon", "HLT"),
        # L1Muon = cms.InputTag("simGmtStage2Digis","",newProcessName),
        L2Muon = cms.InputTag("hltL2MuonCandidates", "", newProcessName),

        mvaFile_B_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter2FromL1_0.xml" % version),
        mvaFile_B_1 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter2FromL1_1.xml" % version),
        mvaFile_B_2 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter2FromL1_2.xml" % version),
        mvaFile_B_3 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter2FromL1_3.xml" % version),
        mvaFile_E_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter2FromL1_0.xml" % version),
        mvaFile_E_1 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter2FromL1_1.xml" % version),
        mvaFile_E_2 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter2FromL1_2.xml" % version),
        mvaFile_E_3 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter2FromL1_3.xml" % version),

        mvaScaleMean_B = cms.vdouble( getattr(_mvaScale, "%sBarrel_hltIter2FromL1_ScaleMean" % version) ),
        mvaScaleStd_B  = cms.vdouble( getattr(_mvaScale, "%sBarrel_hltIter2FromL1_ScaleStd" % version) ),
        mvaScaleMean_E = cms.vdouble( getattr(_mvaScale, "%sEndcap_hltIter2FromL1_ScaleMean" % version) ),
        mvaScaleStd_E  = cms.vdouble( getattr(_mvaScale, "%sEndcap_hltIter2FromL1_ScaleStd" % version) ),

        doSort = cms.bool(doSort),
        nSeedsMax_B = cms.int32(nSeedsMax_B[5]),
        nSeedsMax_E = cms.int32(nSeedsMax_E[5]),

        mvaCut_B = cms.double(mvaCuts_B[5]),
        mvaCut_E = cms.double(mvaCuts_E[5])
    )

    process.hltIter3IterL3FromL1MuonPixelSeedsFiltered = cms.EDProducer("MuonHLTSeedMVAClassifier",

        # Reject all seeds
        rejectAll = cms.bool(True),

        src    = cms.InputTag("hltIter3IterL3FromL1MuonPixelSeeds", "", newProcessName),
        L1Muon = cms.InputTag("hltGtStage2Digis", "Muon", newProcessName),
        # L1Muon = cms.InputTag("hltGtStage2Digis", "Muon", "HLT"),
        # L1Muon = cms.InputTag("simGmtStage2Digis","",newProcessName),
        L2Muon = cms.InputTag("hltL2MuonCandidates", "", newProcessName),

        mvaFile_B_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter3FromL1_0.xml" % version),
        mvaFile_B_1 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter3FromL1_1.xml" % version),
        mvaFile_B_2 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter3FromL1_2.xml" % version),
        mvaFile_B_3 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sBarrel_hltIter3FromL1_3.xml" % version),
        mvaFile_E_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter3FromL1_0.xml" % version),
        mvaFile_E_1 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter3FromL1_1.xml" % version),
        mvaFile_E_2 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter3FromL1_2.xml" % version),
        mvaFile_E_3 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifier/data/%sEndcap_hltIter3FromL1_3.xml" % version),

        mvaScaleMean_B = cms.vdouble( getattr(_mvaScale, "%sBarrel_hltIter3FromL1_ScaleMean" % version) ),
        mvaScaleStd_B  = cms.vdouble( getattr(_mvaScale, "%sBarrel_hltIter3FromL1_ScaleStd" % version) ),
        mvaScaleMean_E = cms.vdouble( getattr(_mvaScale, "%sEndcap_hltIter3FromL1_ScaleMean" % version) ),
        mvaScaleStd_E  = cms.vdouble( getattr(_mvaScale, "%sEndcap_hltIter3FromL1_ScaleStd" % version) ),

        doSort = cms.bool(doSort),
        nSeedsMax_B = cms.int32(nSeedsMax_B[6]),
        nSeedsMax_E = cms.int32(nSeedsMax_E[6]),

        mvaCut_B = cms.double(mvaCuts_B[6]),
        mvaCut_E = cms.double(mvaCuts_E[6])
    )

    # -- Track Candidates
    process.hltIter2IterL3MuonCkfTrackCandidates.src       = cms.InputTag("hltIter2IterL3MuonPixelSeedsFiltered", "", newProcessName)
    process.hltIter3IterL3MuonCkfTrackCandidates.src       = cms.InputTag("hltIter3IterL3MuonPixelSeedsFiltered", "", newProcessName)
    process.hltIter2IterL3FromL1MuonCkfTrackCandidates.src = cms.InputTag("hltIter2IterL3FromL1MuonPixelSeedsFiltered", "", newProcessName)
    process.hltIter3IterL3FromL1MuonCkfTrackCandidates.src = cms.InputTag("hltIter3IterL3FromL1MuonPixelSeedsFiltered", "", newProcessName)

    # -- Sequences
    process.HLTIterativeTrackingIteration2ForIterL3Muon = cms.Sequence(
        process.hltIter2IterL3MuonClustersRefRemoval+
        process.hltIter2IterL3MuonMaskedMeasurementTrackerEvent+
        process.hltIter2IterL3MuonPixelLayerTriplets+
        process.hltIter2IterL3MuonPixelClusterCheck+
        process.hltIter2IterL3MuonPixelHitDoublets+
        process.hltIter2IterL3MuonPixelHitTriplets+
        process.hltIter2IterL3MuonPixelSeeds+
        process.hltIter2IterL3MuonPixelSeedsFiltered+  # HERE
        process.hltIter2IterL3MuonCkfTrackCandidates+
        process.hltIter2IterL3MuonCtfWithMaterialTracks+
        process.hltIter2IterL3MuonTrackCutClassifier+
        process.hltIter2IterL3MuonTrackSelectionHighPurity
    )

    process.HLTIterativeTrackingIteration3ForIterL3Muon = cms.Sequence(
        process.hltIter3IterL3MuonClustersRefRemoval+
        process.hltIter3IterL3MuonMaskedMeasurementTrackerEvent+
        process.hltIter3IterL3MuonPixelLayerPairs+
        process.hltIter3IterL3MuonL2Candidates+
        process.hltIter3IterL3MuonTrackingRegions+
        process.hltIter3IterL3MuonPixelClusterCheck+
        process.hltIter3IterL3MuonPixelHitDoublets+
        process.hltIter3IterL3MuonPixelSeeds+
        process.hltIter3IterL3MuonPixelSeedsFiltered+  # HERE
        process.hltIter3IterL3MuonCkfTrackCandidates+
        process.hltIter3IterL3MuonCtfWithMaterialTracks+
        process.hltIter3IterL3MuonTrackCutClassifier+
        process.hltIter3IterL3MuonTrackSelectionHighPurity
    )

    process.HLTIterativeTrackingIteration2ForIterL3FromL1Muon = cms.Sequence(
        process.hltIter2IterL3FromL1MuonClustersRefRemoval+
        process.hltIter2IterL3FromL1MuonMaskedMeasurementTrackerEvent+
        process.hltIter2IterL3FromL1MuonPixelLayerTriplets+
        process.hltIter2IterL3FromL1MuonPixelClusterCheck+
        process.hltIter2IterL3FromL1MuonPixelHitDoublets+
        process.hltIter2IterL3FromL1MuonPixelHitTriplets+
        process.hltIter2IterL3FromL1MuonPixelSeeds+
        process.hltIter2IterL3FromL1MuonPixelSeedsFiltered+  # HERE
        process.hltIter2IterL3FromL1MuonCkfTrackCandidates+
        process.hltIter2IterL3FromL1MuonCtfWithMaterialTracks+
        process.hltIter2IterL3FromL1MuonTrackCutClassifier+
        process.hltIter2IterL3FromL1MuonTrackSelectionHighPurity
    )

    process.HLTIterativeTrackingIteration3ForIterL3FromL1Muon = cms.Sequence(
        process.hltIter3IterL3FromL1MuonClustersRefRemoval+
        process.hltIter3IterL3FromL1MuonMaskedMeasurementTrackerEvent+
        process.hltIter3IterL3FromL1MuonPixelLayerPairs+
        process.hltIter3IterL3FromL1MuonTrackingRegions+
        process.hltIter3IterL3FromL1MuonPixelClusterCheck+
        process.hltIter3IterL3FromL1MuonPixelHitDoublets+
        process.hltIter3IterL3FromL1MuonPixelSeeds+
        process.hltIter3IterL3FromL1MuonPixelSeedsFiltered+  # HERE
        process.hltIter3IterL3FromL1MuonCkfTrackCandidates+
        process.hltIter3IterL3FromL1MuonCtfWithMaterialTracks+
        process.hltIter3IterL3FromL1MuonTrackCutClassifier+
        process.hltIter3IterL3FromL1MuonTrackSelectionHighPurity
    )

    # -- DQMOutput and Ntupler
    # if hasattr(process, "DQMOutput"):
    #     del process.DQMOutput

    if hasattr(process, "dqmOutput"):
        process.dqmOutput.fileName = cms.untracked.string("DQMIO_%s_%s.root" % (version, WPName) )

    if hasattr(process, "TFileService"):
        process.TFileService.fileName = cms.string("ntuple_%s_%s.root" % (version, WPName) )

    if hasattr(process, "writeDataset"):
        process.writeDataset.fileName = cms.untracked.string("edmOutput_%s_%s.root" % (version, WPName) )

    if hasattr(process, "ntupler"):
        process.ntupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter2IterL3MuonPixelSeedsFiltered",       "", newProcessName)
        # process.ntupler.hltIter3IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter3IterL3MuonPixelSeedsFiltered",       "", newProcessName)
        process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter2IterL3FromL1MuonPixelSeedsFiltered", "", newProcessName)
        # process.ntupler.hltIter3IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter3IterL3FromL1MuonPixelSeedsFiltered", "", newProcessName)

    if hasattr(process, "seedNtupler"):
        process.seedNtupler.hltIter2IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter2IterL3MuonPixelSeedsFiltered",       "", newProcessName)
        # process.seedNtupler.hltIter3IterL3MuonPixelSeeds       = cms.untracked.InputTag("hltIter3IterL3MuonPixelSeedsFiltered",       "", newProcessName)
        process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter2IterL3FromL1MuonPixelSeedsFiltered", "", newProcessName)
        # process.seedNtupler.hltIter3IterL3FromL1MuonPixelSeeds = cms.untracked.InputTag("hltIter3IterL3FromL1MuonPixelSeedsFiltered", "", newProcessName)

    return process
