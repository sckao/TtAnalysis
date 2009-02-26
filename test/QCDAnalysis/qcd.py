import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
#process.load("MagneticField.Engine.volumeBasedMagneticField_1103l_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'IDEAL_V9::All'

#from TrackingTools.TrackAssociator.default_cfi import * 

process.source = cms.Source("PoolSource",
    debugFlag = cms.untracked.bool(False),
    debugVebosity = cms.untracked.uint32(10),
    skipEvents = cms.untracked.uint32(SKIPNUM),
    fileNames = cms.untracked.vstring(
#'file:/data/top/sckao/MuEnrichedQCD/QCD_PAT1_1.root',
#'file:/data/top/sckao/MuEnrichedQCD/QCD_PAT1_2.root',
#'file:/data/top/sckao/MuEnrichedQCD/QCD_PAT1_3.root',
'file:/data/top/sckao/MuEnrichedQCDSkim/FuckPAT.root'
#'file:/data/top/sckao/Tt210Full/RelValTt214_full4.root'
    )
)

# replace the source files from a file list
#import PhysicsTools.TtAnalysis.qcdlist_cff as fileList
import PhysicsTools.TtAnalysis.qcdlist_skim1 as fileList
process.source.fileNames = fileList.fileNames

process.maxEvents = cms.untracked.PSet(
    # 1 jet filter ; 2841508 for 100 /pb , every 56830 => 2 /pb
    input = cms.untracked.int32(56830)
)
process.MessageLogger = cms.Service("MessageLogger")

process.ttAna = cms.EDFilter("TtAnalysis",

    #TrackAssociatorParameterBlock,
    #TrackAssociatorParameters, 
    debug    = cms.untracked.bool(False),
    btag     = cms.untracked.bool(False),
    needTree = cms.untracked.bool(False),
    trigOn   = cms.untracked.bool(False),
    rootFileName = cms.untracked.string('/data/top/sckao/Fall08QCDAna/ANAFILE'),
    genParticles = cms.InputTag("genParticles"),
    genJetSource = cms.InputTag("iterativeCone5GenJets"),
    electronSource = cms.InputTag("selectedLayer1Electrons"),
    photonSource   = cms.InputTag("selectedLayer1Photons"),
    jetSource      = cms.InputTag("selectedLayer1Jets"),
    metSource      = cms.InputTag("selectedLayer1METs"),
    muonSource     = cms.InputTag("selectedLayer1Muons"),
    caloSource     = cms.InputTag("towerMaker"),
    triggerSource  = cms.InputTag("TriggerResults","","HLT"),
    #triggerSource  = cms.InputTag("TriggerResults","","PAT"),
    recoMuons      = cms.untracked.string('muons'),
    leptonFlavour  = cms.string('muon')
)

process.jetAna = cms.EDFilter("JetAnalysis",

    debug    = cms.untracked.bool(False),
    rootFileName = cms.untracked.string('/data/top/sckao/Fall08QCDAna/JANAFILE'),
    genParticles = cms.InputTag("genParticles"),
    genJetSource = cms.InputTag("iterativeCone5GenJets"),
    electronSource = cms.InputTag("selectedLayer1Electrons"),
    jetSource      = cms.InputTag("selectedLayer1Jets"),
    metSource      = cms.InputTag("selectedLayer1METs"),
    muonSource     = cms.InputTag("selectedLayer1Muons"),
    caloSource     = cms.InputTag("towerMaker"),
    recoMuons      = cms.untracked.string('muons'),
)

process.muAna = cms.EDFilter("MuonAnalysis",
    
    debug    = cms.untracked.bool(False),
    rootFileName = cms.untracked.string('/data/top/sckao/Fall08QCDAna/MANAFILE'),
    genParticles = cms.InputTag("genParticles"),
    genJetSource = cms.InputTag("iterativeCone5GenJets"),
    electronSource = cms.InputTag("selectedLayer1Electrons"),
    jetSource      = cms.InputTag("selectedLayer1Jets"),
    muonSource     = cms.InputTag("selectedLayer1Muons"),
    metSource      = cms.InputTag("selectedLayer1METs"),
    caloSource     = cms.InputTag("towerMaker"),
    recoMuons      = cms.untracked.string('muons'),
)   


process.p = cms.Path(process.ttAna + process.jetAna + proecess.muAna )
#process.ttAna.TrackAssociatorParameters.useEcal = False
#process.ttAna.TrackAssociatorParameters.useHcal = False
#process.ttAna.TrackAssociatorParameters.useCalo = True
#process.ttAna.TrackAssociatorParameters.useHO = False
#process.ttAna.TrackAssociatorParameters.useMuon = False

