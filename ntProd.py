import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load("Configuration.StandardSequences.MagneticField_cff")

#process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")


process.source = cms.Source("PoolSource",
    #debugFlag = cms.untracked.bool(False),
    #debugVebosity = cms.untracked.uint32(10),
    #skipEvents = cms.untracked.uint32(1683),
    fileNames = cms.untracked.vstring(
#'dcache:/pnfs/cms/WAX/resilient/sckao/PAT361/tt_pat356_1_1.root'
'file:/uscms_data/d2/sckao/tt_pat387D6.root'

    ), duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

)

# replace the source files from a file list
#import TopPhysics.TtAnalysis.ttPATfile_list as fileList
#process.source.fileNames = fileList.fileNames

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.MessageLogger = cms.Service("MessageLogger")

process.ttNtp = cms.EDAnalyzer("TtNtupleProd",

    # General Setup
    debug    = cms.untracked.bool(False),
    isData   = cms.untracked.bool(False),
    eventId  = cms.untracked.int32(0),
    #rootFileName = cms.untracked.string('/uscms_data/d2/sckao/Ntp361/tt_361_0jpfA.root'),
    rootFileName = cms.untracked.string('tt_387_0j.root'),
    numberOfJets = cms.untracked.int32(0),
    # primary vertex selection
    pvSource = cms.InputTag('offlinePrimaryVertices'),
    pvNDOF   = cms.untracked.double(4),
    pvMaxZ   = cms.untracked.double(15.),
    pvMaxRho = cms.untracked.double(2.),
    beamSpotSource = cms.InputTag('offlineBeamSpot'),
    # trigger Source
    trigSource = cms.InputTag('patTriggerEvent'),
    trigTag    = cms.untracked.string('HLT_Mu9'),

    # Jet/MET Setup
    #btag     = cms.untracked.bool(False),
    bTagCut  = cms.untracked.double(5),
    bTagAlgo = cms.untracked.string('trackCountingHighEffBJetTags'),
    ##                          ( pT, eta, JES  emF  )
    jetSetup       = cms.vdouble( 20, 2.4, 1.00, 0.01 ),
    jetSource      = cms.InputTag("selectedPatJetsPF"),
    metSource      = cms.InputTag("patMETsPF"),
    #jetSource      = cms.InputTag("selectedPatJets"),
    #metSource      = cms.InputTag("patMETs"),
    #caloSource     = cms.InputTag("towerMaker"),
    recoMetSource  = cms.InputTag("tcMet"),
    #genJetSource   = cms.InputTag("ak5GenJets"),
    # Muon Setup
    ##                          ( pT, eta, Iso, nHits, chi2,  d0(Bsp)  )
    muSetup        = cms.vdouble( 15, 2.1, 0.15,   11,   10.,  0.02 ),
    muonSource     = cms.InputTag("selectedPatMuonsPF"),
    recoMuons      = cms.untracked.string('muons'),
    # e/gamma veto Setup
    ##                          ( pT, eta, Iso, H/E,  E/P )
    eleSetup        = cms.vdouble(15, 2.5, 0.2, 0.02, 0.8 ),
    #eleSetup        = cms.vdouble(15, 2.5, 0.2, 0.02, 0.8 ),
    electronSource = cms.InputTag("selectedPatElectronsPF"),
    photonSource   = cms.InputTag("selectedPatPhotons"),

    # other setup
    genParticles = cms.InputTag("genParticles"),
)

process.p = cms.Path( process.ttNtp )
#process.ttAna.TrackAssociatorParameters.useEcal = False
#process.ttAna.TrackAssociatorParameters.useHcal = False
#process.ttAna.TrackAssociatorParameters.useCalo = True
#process.ttAna.TrackAssociatorParameters.useHO = False
#process.ttAna.TrackAssociatorParameters.useMuon = False

