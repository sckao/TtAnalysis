import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'IDEAL_V11::All'

from TrackingTools.TrackAssociator.default_cfi import * 

process.source = cms.Source("PoolSource",
    debugFlag = cms.untracked.bool(False),
    debugVebosity = cms.untracked.uint32(10),
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
#'file:/data/top/sckao/Fall08HLTWJets/WJets_PATSkim1_3.root'
'file:/home/cms/sckao/Data/FastPat2213/WJets_PAT_Test.root',
'file:/home/cms/sckao/Data/FastPat2213/WJets_PAT_Test2.root'

    )
)

# replace the source files from a file list
#import PhysicsTools.TtAnalysis.wjetslist_cff as fileList
#import PhysicsTools.TtAnalysis.wjetslist_skim as fileList
#rocess.source.fileNames = fileList.fileNames

process.maxEvents = cms.untracked.PSet(
    # for 1J Skim 100 pb-1 => 44000
    input = cms.untracked.int32(44000)
    #input = cms.untracked.int32(-1)
)
process.MessageLogger = cms.Service("MessageLogger")

process.ttAna = cms.EDFilter("TtAnalysis",

    TrackAssociatorParameterBlock,
    TrackAssociatorParameters, 
    debug    = cms.untracked.bool(False),
    btag     = cms.untracked.bool(True),
    bTagCut  = cms.untracked.double(2),
    #bTagAlgo = cms.untracked.string('softMuonBJetTags'),
    #bTagAlgo = cms.untracked.string('jetProbabilityBJetTags'),
    bTagAlgo = cms.untracked.string("trackCountingHighEffBJetTags"),
    needTree = cms.untracked.bool(False),
    trigOn   = cms.untracked.bool(True),
    rootFileName = cms.untracked.string('wjets_pat2_SCJ_TC_btag_tt.root'),
    genParticles = cms.InputTag("genParticles"),
    #genJetSource = cms.InputTag("iterativeCone5GenJets"),
    genJetSource = cms.InputTag("sisCone5GenJets"),
    electronSource = cms.InputTag("cleanLayer1Electrons"),
    photonSource   = cms.InputTag("cleanLayer1Photons"),
    jetSource      = cms.InputTag("cleanLayer1Jets"),
    jptSource      = cms.InputTag("ZSPJetCorJetIcone5"),
    tcMetSource    = cms.InputTag("tcMet"),
    metSource      = cms.InputTag("layer1METsTC"),
    genmetSource   = cms.InputTag("genMet"),
    muonSource     = cms.InputTag("cleanLayer1Muons"),
    caloSource     = cms.InputTag("towerMaker"),
    triggerSource  = cms.InputTag("TriggerResults","","HLT"),
    #triggerSource  = cms.InputTag("TriggerResults","","PAT"),
    #recoMuons      = cms.untracked.string('paramMuons'),
    recoMuons      = cms.untracked.string('muons'),
    recoAlgo       = cms.untracked.string('zero'),
)

process.ttAna1 = cms.EDFilter("TtAnalysis",

    TrackAssociatorParameterBlock,
    TrackAssociatorParameters, 
    debug    = cms.untracked.bool(False),
    btag     = cms.untracked.bool(False),
    bTagCut  = cms.untracked.double(2),
    bTagAlgo = cms.untracked.string("trackCountingHighEffBJetTags"),
    needTree = cms.untracked.bool(False),
    trigOn   = cms.untracked.bool(True),
    rootFileName = cms.untracked.string('wjets_pat2_SCJ_TC_tt.root'),
    genParticles = cms.InputTag("genParticles"),
    #genJetSource = cms.InputTag("iterativeCone5GenJets"),
    genJetSource = cms.InputTag("sisCone5GenJets"),
    electronSource = cms.InputTag("cleanLayer1Electrons"),
    photonSource   = cms.InputTag("cleanLayer1Photons"),
    jetSource      = cms.InputTag("cleanLayer1Jets"),
    jptSource      = cms.InputTag("ZSPJetCorJetIcone5"),
    tcMetSource    = cms.InputTag("tcMet"),
    metSource      = cms.InputTag("layer1METsTC"),
    genmetSource   = cms.InputTag("genMet"),
    muonSource     = cms.InputTag("cleanLayer1Muons"),
    caloSource     = cms.InputTag("towerMaker"),
    triggerSource  = cms.InputTag("TriggerResults","","HLT"),
    recoMuons      = cms.untracked.string('muons'),
    recoAlgo       = cms.untracked.string('zero'),
)

process.jetAna = cms.EDFilter("JetAnalysis",

    debug    = cms.untracked.bool(False),
    bTagCut  = cms.untracked.double(2),
    #bTagAlgo = cms.untracked.string('softMuonBJetTags'),
    #bTagAlgo = cms.untracked.string('jetProbabilityBJetTags'),
    bTagAlgo = cms.untracked.string("trackCountingHighEffBJetTags"),
    rootFileName = cms.untracked.string('wjets_pat2_SCJ_TC_btag_jm.root'),
    genParticles = cms.InputTag("genParticles"),
    #genJetSource = cms.InputTag("iterativeCone5GenJets"),
    genJetSource = cms.InputTag("sisCone5GenJets"),
    electronSource = cms.InputTag("cleanLayer1Electrons"),
    jetSource      = cms.InputTag("cleanLayer1Jets"),
    metSource      = cms.InputTag("layer1METsTC"),
    genmetSource   = cms.InputTag("genMet"),
    jptSource      = cms.InputTag("ZSPJetCorJetIcone5"),
    tcMetSource    = cms.InputTag("tcMet"),
    muonSource     = cms.InputTag("cleanLayer1Muons"),
    caloSource     = cms.InputTag("towerMaker"),
    recoMuons      = cms.untracked.string('muons')

)

process.muAna = cms.EDFilter("MuonAnalysis",

    debug    = cms.untracked.bool(False),
    rootFileName = cms.untracked.string('/data/top/sckao/wjets_IsoMuAnalysis.root'),
    genParticles = cms.InputTag("genParticles"),
    genJetSource = cms.InputTag("iterativeCone5GenJets"),
    electronSource = cms.InputTag("selectedLayer1Electrons"),
    jetSource      = cms.InputTag("selectedLayer1Jets"),
    muonSource     = cms.InputTag("selectedLayer1Muons"),
    metSource      = cms.InputTag("selectedLayer1METs"),
    caloSource     = cms.InputTag("towerMaker"),
    recoMuons      = cms.untracked.string('muons')
)

process.p = cms.Path( process.ttAna + process.ttAna1 +process.jetAna )
#process.ttAna.TrackAssociatorParameters.useEcal = False
#process.ttAna.TrackAssociatorParameters.useHcal = False
#process.ttAna.TrackAssociatorParameters.useCalo = True
#process.ttAna.TrackAssociatorParameters.useHO = False
#process.ttAna.TrackAssociatorParameters.useMuon = False

