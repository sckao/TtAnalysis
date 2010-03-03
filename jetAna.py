import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'IDEAL_V12::All'

# for match muon in a jet by using trackAssociator
#from TrackingTools.TrackAssociator.default_cfi import * 

process.source = cms.Source("PoolSource",
    debugFlag = cms.untracked.bool(False),
    debugVebosity = cms.untracked.uint32(10),
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(

'file:/uscms_data/d2/sckao/PAT336/tt171_pat336.root'
#'dcache:/pnfs/cms/WAX/resilient/sckao/QCDPAT336/qcd_pat336.root'
    )
)

# replace the source files from a file list
#import TopPhysics.TtAnalysis.qcdPATfile_list as fileList
#process.source.fileNames = fileList.fileNames


process.maxEvents = cms.untracked.PSet(
    # 100/pb => 9191 
    input = cms.untracked.int32(20000)
    #input = cms.untracked.int32(24000)
)
process.MessageLogger = cms.Service("MessageLogger")

#from PhysicsTools.PatAlgos.tools.trigTools import switchOffTriggerMatchingOld
#switchOffTriggerMatchingOld( process )

process.jetAna = cms.EDAnalyzer("JetAnalysis",

    debug    = cms.untracked.bool(False),
    isData   = cms.untracked.bool(False),
    rootFileName = cms.untracked.string('tt171_jetInfo.root'),
    ##                          ( Et, eta, JES )
    jetSetup       = cms.vdouble( 20, 2.6, 1.0 ),
    bTagCut  = cms.untracked.double(5),
    bTagAlgo = cms.untracked.string("trackCountingHighEffBJetTags"),
    #bTagAlgo = cms.untracked.string('softMuonBJetTags'),
    #bTagAlgo = cms.untracked.string('jetProbabilityBJetTags'),

    jetSource      = cms.InputTag("cleanLayer1Jets"),
    metSource      = cms.InputTag("layer1METs"),
    recoMetSource  = cms.InputTag("tcMet"),
    muSetup        = cms.vdouble( 5, 2.4, 0.3 ),
    muonSource     = cms.InputTag("cleanLayer1Muons"),
    eleSetup        = cms.vdouble( 5, 2.4, 0.2, 0.02, 0.9 ),
    electronSource = cms.InputTag("cleanLayer1Electrons"),
    genParticles = cms.InputTag("genParticles"),
    genJetSource = cms.InputTag("ak5GenJets"),
    caloSource     = cms.InputTag("towerMaker")
)

process.lepAna = cms.EDAnalyzer("MuonAnalysis",

    debug    = cms.untracked.bool(False),
    isData   = cms.untracked.bool(False),
    rootFileName = cms.untracked.string('tt171_MuInfo.root'),

    ## Muon Setup
    ##                          ( pT, eta, Iso )
    muSetup        = cms.vdouble( 15, 2.1, 0.1 ),
    recoMuons      = cms.untracked.string('muons'),
    muonSource     = cms.InputTag("cleanLayer1Muons"),
    ## Electron Setup
    ##                          ( pT, eta, Iso, H/E,  E/P )
    eleSetup        = cms.vdouble( 5, 2.4, 0.9, 0.02, 0.9 ),
    electronSource = cms.InputTag("cleanLayer1Electrons"),
    genParticles = cms.InputTag("genParticles")
)


#process.p = cms.Path( process.lepAna )
process.p = cms.Path( process.jetAna + process.lepAna )

#process.ttAna.TrackAssociatorParameters.useEcal = False
#process.ttAna.TrackAssociatorParameters.useHcal = False
#process.ttAna.TrackAssociatorParameters.useCalo = True
#process.ttAna.TrackAssociatorParameters.useHO = False
#process.ttAna.TrackAssociatorParameters.useMuon = False

