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
process.GlobalTag.globaltag = 'IDEAL_V5::All'

from TrackingTools.TrackAssociator.default_cfi import * 

process.source = cms.Source("PoolSource",
    debugFlag = cms.untracked.bool(False),
    debugVebosity = cms.untracked.uint32(10),
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
'file:/data/top/sckao/Tt210Full/RelValTt214_full1.root',
'file:/data/top/sckao/Tt210Full/RelValTt214_full2.root',
'file:/data/top/sckao/Tt210Full/RelValTt214_full3.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.MessageLogger = cms.Service("MessageLogger")

process.ttAna = cms.EDFilter("TtAnalysis",

    TrackAssociatorParameterBlock,
    TrackAssociatorParameters, 
    debug    = cms.untracked.bool(False),
    needTree = cms.untracked.bool(False),
    rootFileName = cms.untracked.string('tt_test_full210.root'),
    genParticles = cms.InputTag("genParticles"),
    genJetSource = cms.InputTag("iterativeCone5GenJets"),
    electronSource = cms.InputTag("selectedLayer1Electrons"),
    photonSource   = cms.InputTag("selectedLayer1Photons"),
    jetSource      = cms.InputTag("selectedLayer1Jets"),
    metSource      = cms.InputTag("selectedLayer1METs"),
    muonSource     = cms.InputTag("selectedLayer1Muons"),
    caloSource     = cms.InputTag("towerMaker"),
    recoMuons      = cms.untracked.string('muons'),
    leptonFlavour  = cms.string('muon')
)

process.p = cms.Path(process.ttAna)
process.ttAna.TrackAssociatorParameters.useEcal = False
process.ttAna.TrackAssociatorParameters.useHcal = False
process.ttAna.TrackAssociatorParameters.useCalo = True
process.ttAna.TrackAssociatorParameters.useHO = False
process.ttAna.TrackAssociatorParameters.useMuon = False

