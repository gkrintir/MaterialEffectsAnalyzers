import FWCore.ParameterSet.Config as cms

EnergyLossValidation = cms.EDAnalyzer('EnergyLossValidation', 
     SimHitTags = cms.VInputTag(
        cms.InputTag('famosSimHits', 'TrackerHits')
        #cms.InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"),
        #cms.InputTag("g4SimHits", "TrackerHitsPixelBarrelHighTof"),
        #cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof"), 
        #cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapHighTof"),
        #cms.InputTag("g4SimHits", "TrackerHitsTECLowTof"), 
        #cms.InputTag("g4SimHits", "TrackerHitsTECHighTof"),
        #cms.InputTag("g4SimHits", "TrackerHitsTIBLowTof"), 
        #cms.InputTag("g4SimHits", "TrackerHitsTIBHighTof"),
        #cms.InputTag("g4SimHits", "TrackerHitsTIDLowTof"), 
        #cms.InputTag("g4SimHits", "TrackerHitsTIDHighTof"),
        #cms.InputTag("g4SimHits", "TrackerHitsTOBLowTof"), 
        #cms.InputTag("g4SimHits", "TrackerHitsTOBHighTof")
     ),
     particleSelection = cms.PSet(
        # Particles with the given PDG code(s) will be considered
        particleTypes = cms.untracked.vint32(211),
        # Sign sensitive; should anti-particles be analyzed, please set 'True'
        selectAntiParticle = cms.untracked.bool(True)
     ),
     bins_p =  cms.untracked.vdouble(.4, .8, 1.2, 2.1, 3.4, 5.2),
     #bins_E =  cms.untracked.vdouble( .4, 1.2, 1.5, 1.9, 2.1),
     particleFilter = cms.PSet(
     # Particles with eta < etaMim are not selected 
     etaMin = cms.untracked.double(-2.5),
     # Particles with eta > etaMax are not selected
     etaMax = cms.untracked.double(2.5),
     # Particles with trans. momentum < pTMin (GeV/c) are not selected
     pTMin = cms.untracked.double(0.05),
     # Particles with trans. momentum > pTMax (GeV/c) are not selected
     pTMax = cms.untracked.double(100.0),
     # Particles with momentum < pMin (GeV/c) are not selected
     pMin = cms.untracked.double(0.0),
     # Particles with momentum > pMax (GeV/c) are not selected
     pMax = cms.untracked.double(0.0),
     # Particles with energy < EMin (GeV) are not selected
     EMin = cms.untracked.double(0.0),
     # Particles with energy > EMax (GeV) are not selected
     EMax = cms.untracked.double(0.0),
     # Particles with the given PDG code(s) will be rejected
     pdgIdsToFilter = cms.untracked.vint32(2212),
     # Sign sensitive; should anti-particles be rejected, please set 'True'
     filterAntiParticle = cms.untracked.bool(True)
     ),
     histo1DFormating = cms.VPSet (
        cms.PSet(
         title = cms.untracked.string('test'),
         name = cms.untracked.string('default'),
         labelx = cms.untracked.string('dEdX'),
         labely = cms.untracked.string(''),
         rangex = cms.untracked.vdouble(200, 0., 1e-2)
         ),
        cms.PSet(
         title = cms.untracked.string('test'),
         name = cms.untracked.string('default'),
         labelx = cms.untracked.string('dEdX'),
         labely = cms.untracked.string(''),
         rangex = cms.untracked.vdouble(200, 0., 1e-2)
         ),
     ),
     histo2DFormating = cms.VPSet (
        cms.PSet(
         title = cms.untracked.string('test'),
         name = cms.untracked.string('default'),
         labelx = cms.untracked.string('p'),
         labely = cms.untracked.string('dEdX'),
         rangex = cms.untracked.vdouble(60., 0., 100.),
         rangey = cms.untracked.vdouble(200., 0., 1e-2),
         )
     )
     
 )                          

    



