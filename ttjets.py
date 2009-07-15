import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'IDEAL_V12::All'

# for match muon in a jet by using trackAssociator
#from TrackingTools.TrackAssociator.default_cfi import * 

process.source = cms.Source("PoolSource",
    debugFlag = cms.untracked.bool(False),
    debugVebosity = cms.untracked.uint32(10),
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(

'file:/home/cms/sckao/Top/CMSSW_2_2_13/src/TopPhysics/TtAnalysis/test/TMass/patTest_SCJ_TC.root'
#'file:/home/cms/sckao/Top/CMSSW_2_2_13/src/TopPhysics/TtAnalysis/test/TMass/patTest_ICJ_TC.root'
#'file:/home/cms/sckao/Top/CMSSW_2_2_13/src/TopPhysics/TtAnalysis/test/TMass/patTest_oldMET_skim.root'
#'file:/home/cms/sckao/Data/FastPat2213/WJets_PAT_Test.root'

    )
)

# replace the source files from a file list
#import TopPhysics.TtAnalysis.ttjetslist_skim2 as fileList
#process.source.fileNames = fileList.fileNames

process.maxEvents = cms.untracked.PSet(
    # 100/pb => 9191 
    input = cms.untracked.int32(9191)
    #input = cms.untracked.int32(24000)
)
process.MessageLogger = cms.Service("MessageLogger")

#from PhysicsTools.PatAlgos.tools.trigTools import switchOffTriggerMatchingOld
#switchOffTriggerMatchingOld( process )

process.ttAna = cms.EDAnalyzer("TtAnalysis",

    # for matching Muon in a Jet
    #TrackAssociatorParameterBlock,
    #TrackAssociatorParameters, 
    debug    = cms.untracked.bool(False),
    btag     = cms.untracked.bool(True),
    bTagCut  = cms.untracked.double(2),
    #bTagAlgo = cms.untracked.string('softMuonBJetTags'),
    bTagAlgo = cms.untracked.string("trackCountingHighEffBJetTags"),
    needTree = cms.untracked.bool(False),
    trigOn   = cms.untracked.bool(False),
    rootFileName = cms.untracked.string('pat2_SCJ_TC_btag_tt.root'),
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

process.ttAna1 = cms.EDAnalyzer("TtAnalysis",

    # for matching Muon in a Jet
    #TrackAssociatorParameterBlock,
    #TrackAssociatorParameters, 
    debug    = cms.untracked.bool(False),
    btag     = cms.untracked.bool(False),
    bTagCut  = cms.untracked.double(2),
    bTagAlgo = cms.untracked.string("trackCountingHighEffBJetTags"),
    needTree = cms.untracked.bool(False),
    trigOn   = cms.untracked.bool(False),
    rootFileName = cms.untracked.string('pat2_SCJ_TC_tt.root'),
    genParticles = cms.InputTag("genParticles"),
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

process.jetAna = cms.EDAnalyzer("JetAnalysis",

    debug    = cms.untracked.bool(False),
    bTagCut  = cms.untracked.double(2),
    #bTagAlgo = cms.untracked.string('softMuonBJetTags'),
    bTagAlgo = cms.untracked.string("trackCountingHighEffBJetTags"),
    #bTagAlgo = cms.untracked.string('jetProbabilityBJetTags'),
    rootFileName = cms.untracked.string('pat2_SCJ_TC_btag_jm.root'),
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
    #recoMuons      = cms.untracked.string('paramMuons')
    recoMuons      = cms.untracked.string('muons')
)

process.muAna = cms.EDAnalyzer("MuonAnalysis",

    debug    = cms.untracked.bool(False),
    rootFileName = cms.untracked.string('ttj_IsoMuAnalysis.root'),
    genParticles = cms.InputTag("genParticles"),
    genJetSource = cms.InputTag("iterativeCone5GenJets"),
    electronSource = cms.InputTag("cleanLayer1Electrons"),
    jetSource      = cms.InputTag("cleanLayer1Jets"),
    muonSource     = cms.InputTag("cleanLayer1Muons"),
    metSource      = cms.InputTag("layer1METs"),
    caloSource     = cms.InputTag("towerMaker"),
    recoMuons      = cms.untracked.string('paramMuons')
)


#process.p = cms.Path( process.ttAna +  process.jetAna + process.muAna )
#process.p = cms.Path( process.ttAna + process.ttAna1 + process.ttAna2 + process.ttAna3 )
process.p = cms.Path( process.ttAna + process.ttAna1 + process.jetAna )
#process.p = cms.Path( process.ttAna )
#process.ttAna.TrackAssociatorParameters.useEcal = False
#process.ttAna.TrackAssociatorParameters.useHcal = False
#process.ttAna.TrackAssociatorParameters.useCalo = True
#process.ttAna.TrackAssociatorParameters.useHO = False
#process.ttAna.TrackAssociatorParameters.useMuon = False

