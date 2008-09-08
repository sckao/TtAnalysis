import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")

# The number of events to be processed.
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(4100) )
    
# For valgrind studies
# process.ProfilerService = cms.Service("ProfilerService",
#    lastEvent = cms.untracked.int32(13),
#    firstEvent = cms.untracked.int32(3),
#    paths = cms.untracked.vstring('p1')
#)

# Include the RandomNumberGeneratorService definition
process.load("FastSimulation/Configuration/RandomServiceInitialization_cff")

process.load("Configuration.GenProduction.MadGraph_XQCUT20_10TeV_WmEle2J_cff")
# Generate ttbar events
#  process.load("FastSimulation/Configuration/ttbar_cfi")
# Generate multijet events with different ptHAT bins
#  process.load("FastSimulation/Configuration/QCDpt80-120_cfi")
#  process.load("FastSimulation/Configuration/QCDpt600-800_cfi")
# Generate Minimum Bias Events
#  process.load("FastSimulation/Configuration/MinBiasEvents_cfi")
# Generate muons with a flat pT particle gun, and with pT=10.
# process.load("FastSimulation/Configuration/FlatPtMuonGun_cfi")
# replace FlatRandomPtGunSource.PGunParameters.PartID={130}
# Generate di-electrons with pT=35 GeV
# process.load("FastSimulation/Configuration/DiElectrons_cfi")

# Famos sequences (Frontier conditions)
process.load("FastSimulation/Configuration/CommonInputsFake_cff")
process.load("FastSimulation/Configuration/FamosSequences_cff")

# Parametrized magnetic field (new mapping, 4.0 and 3.8T)
#process.load("Configuration.StandardSequences.MagneticField_40T_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = True

# If you want to turn on/off pile-up
process.famosPileUp.PileUpSimulator.averageNumber = 5.0    
# You may not want to simulate everything for your study
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True
# process.famosSimHits.SimulateMuons = False

# Produce Tracks and Clusters
#process.p1 = cms.Path(process.famosWithTracksAndEcalClusters)
process.p1 = cms.Path(process.famosWithEverything)

# To write out events (not need: FastSimulation _is_ fast!)
process.o1 = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string("/data/top/sckao/WmEle2J/WmEle2J_Exclusive_Fast.root"),
    outputCommands = cms.untracked.vstring("keep *",
                                           "drop *_mix_*_*")
    )

process.outpath = cms.EndPath(process.o1)

# Keep the logging output to a nice level #
# process.Timing =  cms.Service("Timing")
# process.load("FWCore/MessageService/MessageLogger_cfi")
# process.MessageLogger.destinations = cms.untracked.vstring("pyDetailedInfo.txt")

