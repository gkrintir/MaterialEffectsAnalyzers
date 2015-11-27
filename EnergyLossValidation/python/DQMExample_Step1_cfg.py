import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

import os

options = VarParsing ('analysis')

options.register(
    'isFSim',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "is fullsim"
)

options.parseArguments()

process = cms.Process('RECODQM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = 'PHYS14_50_V1'

# load DQM
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")

# my analyzer
process.load("Analyzer.MaterialEffectsAnalyzer.EnergyLossValidation_cfi")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# Input source
if options.isFSim:
  FSimfile = open('ppions_FULLSIM_GENSIM_05_10.txt', 'r')
  FSimlist = cms.untracked.vstring()
  FSimlist.extend( [line.strip() for line in FSimfile.read().splitlines()] )
  process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = FSimlist
    )
else:
  process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
      'file:/afs/cern.ch/user/g/gkrintir/github/GenSim/CMSSW_7_3_0/src/Generation_Output/FastSim/GENSIM/mygun_ppions_FASTSIM_SameFULLSIMConditions_noflatpT.root'
       )
)

#For debugging
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('nalysis_mygyn_muons_FASTSIM_SameFULLSIMConditions_noflatpT.root')
                                   )


process.dqmSaver.workflow = "/CMSSW_3_1_1/RelVal/TrigVal_Hgg"

# Path and EndPath definitions
process.dqmoffline_step = cms.Path(process.DQMExample_Step1)
process.dqmsave_step = cms.Path(process.dqmSaver)
#process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
process.schedule = cms.Schedule(
    process.dqmoffline_step,
    #process.DQMoutput_step,
    process.dqmsave_step
    )

# Keep the logging output to a nice level
process.MessageLogger.destinations = ['newE_.txt']
