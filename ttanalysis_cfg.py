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
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
'file:/data/top/sckao/MGTtJets/TtJets_PAT1_0.root'
    )
)

# replace the source files from a file list
#import PhysicsTools.TtAnalysis.ttjetslist_cff as fileList
import PhysicsTools.TtAnalysis.ttjetslist_skim1 as fileList
process.source.fileNames = fileList.fileNames

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(8353)
)
process.MessageLogger = cms.Service("MessageLogger")

process.ttAna = cms.EDAnalyzer("TtAnalysis",

    #TrackAssociatorParameterBlock,
    #TrackAssociatorParameters, 
    debug    = cms.untracked.bool(False),
    btag     = cms.untracked.bool(False),
    needTree = cms.untracked.bool(False),
    trigOn   = cms.untracked.bool(False),
    rootFileName = cms.untracked.string('ttj_1Jskim_NoBTag1.root'),
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
    # recoAlgo : zero, beta, dmMin, ptMin
    recoAlgo       = cms.untracked.string('zero'),
    leptonFlavour  = cms.string('muon')
)

process.jetAna = cms.EDAnalyzer("JetAnalysis",

    debug    = cms.untracked.bool(False),
    rootFileName = cms.untracked.string('ttj_JetEtAnalysis.root'),
    genParticles = cms.InputTag("genParticles"),
    genJetSource = cms.InputTag("iterativeCone5GenJets"),
    electronSource = cms.InputTag("selectedLayer1Electrons"),
    jetSource      = cms.InputTag("selectedLayer1Jets"),
    muonSource     = cms.InputTag("selectedLayer1Muons"),
    caloSource     = cms.InputTag("towerMaker"),
    recoMuons      = cms.untracked.string('muons')
)   
    
process.muAna = cms.EDAnalyzer("MuonAnalysis",
    
    debug    = cms.untracked.bool(False),
    rootFileName = cms.untracked.string('ttj_IsoMuAnalysis.root'),
    genParticles = cms.InputTag("genParticles"),
    genJetSource = cms.InputTag("iterativeCone5GenJets"),
    electronSource = cms.InputTag("selectedLayer1Electrons"),
    jetSource      = cms.InputTag("selectedLayer1Jets"),
    muonSource     = cms.InputTag("selectedLayer1Muons"),
    metSource      = cms.InputTag("selectedLayer1METs"),
    caloSource     = cms.InputTag("towerMaker"),
    recoMuons      = cms.untracked.string('muons')
)   
    

process.p = cms.Path( process.ttAna + process.jetAna + process.muAna )
#process.ttAna.TrackAssociatorParameters.useEcal = False
#process.ttAna.TrackAssociatorParameters.useHcal = False
#process.ttAna.TrackAssociatorParameters.useCalo = True
#process.ttAna.TrackAssociatorParameters.useHO = False
#process.ttAna.TrackAssociatorParameters.useMuon = False

