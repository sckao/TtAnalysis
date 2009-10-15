import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'IDEAL_V12::All'

#from TrackingTools.TrackAssociator.default_cfi import * 

process.source = cms.Source("PoolSource",
    debugFlag = cms.untracked.bool(False),
    debugVebosity = cms.untracked.uint32(10),
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
'file:/home/cms/sckao/Data/FullPat2213/ttPatFast171_1.root'

    ), duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

)

# replace the source files from a file list
#import PhysicsTools.TtAnalysis.wjetslist_cff as fileList
#import PhysicsTools.TtAnalysis.wjetslist_skim as fileList
#rocess.source.fileNames = fileList.fileNames

process.maxEvents = cms.untracked.PSet(
    # for 1J Skim 100 pb-1 => 9067
    input = cms.untracked.int32(1200)
)
process.MessageLogger = cms.Service("MessageLogger")

process.ttNtp = cms.EDFilter("TtNtupleProd",

    debug    = cms.untracked.bool(False),
    btag     = cms.untracked.bool(False),
    bTagCut  = cms.untracked.double(5),
    nbtagged = cms.untracked.int32(0),
    JEScale  = cms.untracked.double(1.),
    bTagAlgo = cms.untracked.string('trackCountingHighEffBJetTags'),
    needTree = cms.untracked.bool(True),
    trigOn   = cms.untracked.bool(True),
    rootFileName = cms.untracked.string('tt171_test2.root'),
    genParticles = cms.InputTag("genParticles"),
    genJetSource = cms.InputTag("sisCone5GenJets"),
    electronSource = cms.InputTag("cleanLayer1Electrons"),
    photonSource   = cms.InputTag("cleanLayer1Photons"),
    jetSource      = cms.InputTag("cleanLayer1JetsSC5"),
    jptSource      = cms.InputTag("ZSPJetCorJetIcone5"),
    caloSource     = cms.InputTag("towerMaker"),
    tcMetSource    = cms.InputTag("tcMet"),
    metSource      = cms.InputTag("layer1METsSC5"),
    genmetSource   = cms.InputTag("genMet"),
    muonSource     = cms.InputTag("cleanLayer1Muons"),
    triggerSource  = cms.InputTag("TriggerResults","","HLT"),
    recoMuons      = cms.untracked.string('muons'),
    #recoAlgo       = cms.untracked.string('zero'),
    recoAlgo       = cms.untracked.string('kConstrain'),
)


process.p = cms.Path( process.ttNtp )
#process.ttAna.TrackAssociatorParameters.useEcal = False
#process.ttAna.TrackAssociatorParameters.useHcal = False
#process.ttAna.TrackAssociatorParameters.useCalo = True
#process.ttAna.TrackAssociatorParameters.useHO = False
#process.ttAna.TrackAssociatorParameters.useMuon = False

