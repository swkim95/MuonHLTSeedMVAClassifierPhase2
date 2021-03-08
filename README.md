# CMS Run3 MuonHLTSeedMVAClassifier

## Setup
```shell
cmsrel CMSSW_11_0_0
cd CMSSW_11_0_0/src
cmsenv

git cms-init
git cms-addpkg HLTrigger/Muon
git clone https://github.com/wonpoint4/MuonHLTSeedMVAClassifier.git HLTrigger/MuonHLTSeedMVAClassifier

scram b -j 8
```

## HLT menu
### Get menu

```shell
# hltGetConfiguration /dev/CMSSW_11_0_0/GRun --globaltag 110X_mcRun3_2021_realistic_v6 \
--path HLTriggerFirstPath,HLT_Mu50_v*,HLTriggerFinalPath \
--input=root://xrootd-cms.infn.it//store/mc/Run3Winter20DRPremixMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8_HCAL/GEN-SIM-DIGI-RAW/110X_mcRun3_2021_realistic_v6-v2/280000/F2982585-CFEB-2345-BB4F-75447796B2F9.root  \
--process MYHLT --full --offline --mc --l1-emulator FullMC \
--timing --prescale none --max-events 100 --output none > HLT_Run3_DY.py
```
### Add this part in the last of menu
```shell
# -- Seer Classifier -- #                                                                            
from HLTrigger.MuonHLTSeedMVAClassifier.customizerForMuonHLTSeeding import *
doSort = False
nSeedMax_B = ( -1, -1, -1, 0, -1, -1, 0 )
nSeedMax_E = ( -1, -1, -1, 0, -1, -1, 0 )
mvaCuts_B = ( 0., 0., 0.00, 1e9, 0., 0.00, 1e9 )
mvaCuts_E = ( 0., 0., 0.00, 1e9, 0., 0.00, 1e9 )
process = customizerFuncForMuonHLTSeeding(process, "MYHLT", "Run3v6", "WP0p00", doSort, nSeedMax_B, nSeedMax_E, mvaCuts_B, mvaCuts_E )
# -- #
```
### Change by Working points
 * If you want to sort seeds by their scores
  - ex) limit the number of seeds at maximum 5 in Iter2, Iter2FromL1 Sequences
```shell
doSort = True
nSeedMax_B = ( -1, -1, 5, 0, -1, 5, 0 )
nSeedMax_E = ( -1, -1, 5, 0, -1, 5, 0 )
process = customizerFuncForMuonHLTSeeding(process, "MYHLT", "Run3v6", "N5", doSort, nSeedMax_B, nSeedMax_E, mvaCuts_B, mvaCuts_E )
```
 * If you want to reject seeds by score threshold
  - ex) limit the scores of seeds at minimum 0.04 in Iter2, Iter2FromL1 Sequences
```shell
doSort = False
mvaCuts_B = ( 0., 0., 0.04, 1e9, 0., 0.04, 1e9 )
mvaCuts_E = ( 0., 0., 0.04, 1e9, 0., 0.04, 1e9 )
process = customizerFuncForMuonHLTSeeding(process, "MYHLT", "Run3v6", "WP0p04", doSort, nSeedMax_B, nSeedMax_E, mvaCuts_B, mvaCuts_E )
```
