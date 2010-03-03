import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'IDEAL_V12::All'

#from TrackingTools.TrackAssociator.default_cfi import * 

process.source = cms.Source("PoolSource",
    debugFlag = cms.untracked.bool(False),
    debugVebosity = cms.untracked.uint32(10),
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
'dcache:/pnfs/cms/WAX/resilient/sckao/PAT336/wj_pat336.root'
#'file:/uscms_data/d2/sckao/Coll2360/col2360_124120pat.root'

    ), duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

)

# replace the source files from a file list
import TopPhysics.TtAnalysis.wjetPATfile_list as fileList
process.source.fileNames = fileList.fileNames

process.maxEvents = cms.untracked.PSet(
    # for 1J Skim 100 pb-1 => 9067
    input = cms.untracked.int32(1000)
)
process.MessageLogger = cms.Service("MessageLogger")

process.ttNtp = cms.EDAnalyzer("TtNtupleProd",

    # General Setup
    debug    = cms.untracked.bool(False),
    needTree = cms.untracked.bool(True),
    trigOn   = cms.untracked.bool(True),
    isData   = cms.untracked.bool(True),
    eventId  = cms.untracked.int32(0),
    rootFileName = cms.untracked.string('wj_336test.root'),
    eventId  = cms.untracked.int32(0),
    # Jet/MET Setup
    btag     = cms.untracked.bool(False),
    bTagCut  = cms.untracked.double(5),
    bTagAlgo = cms.untracked.string('trackCountingHighEffBJetTags'),
    nbtagged = cms.untracked.int32(0),
    ##                          ( pT, eta, JES )
    jetSetup       = cms.vdouble( 30, 2.6, 1.0 ),
    jetSource      = cms.InputTag("cleanLayer1Jets"),
    metSource      = cms.InputTag("layer1METs"),
    caloSource     = cms.InputTag("towerMaker"),
    recoMetSource  = cms.InputTag("tcMet"),
    genJetSource   = cms.InputTag("ak5GenJets"),
    # Muon Setup
    ##                          ( pT, eta, Iso )
    muSetup        = cms.vdouble( 15, 2.1, 0.1 ),
    muonSource     = cms.InputTag("cleanLayer1Muons"),
    recoMuons      = cms.untracked.string('muons'),
    # e/gamma Setup
    ##                          ( pT, eta, Iso, H/E,  E/P )
    eleSetup        = cms.vdouble(15, 2.1, 0.9, 0.02, 0.9 ),
    electronSource = cms.InputTag("cleanLayer1Electrons"),
    photonSource   = cms.InputTag("cleanLayer1Photons"),

    # other setup
    genParticles = cms.InputTag("genParticles"),
    triggerSource  = cms.InputTag("TriggerResults","","HLT"),
    recoAlgo       = cms.untracked.string('kConstrain')
    #recoAlgo       = cms.untracked.string('zero'),
)

process.p = cms.Path( process.ttNtp )
#process.ttAna.TrackAssociatorParameters.useEcal = False
#process.ttAna.TrackAssociatorParameters.useHcal = False
#process.ttAna.TrackAssociatorParameters.useCalo = True
#process.ttAna.TrackAssociatorParameters.useHO = False
#process.ttAna.TrackAssociatorParameters.useMuon = False

