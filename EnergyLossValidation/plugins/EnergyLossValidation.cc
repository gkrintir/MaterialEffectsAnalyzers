// -*- C++ -*-
//
// Package:    Analyzer/DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Analyzer/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Georgios Krintiras
//         Created:  Fri, 13 Mar 2015 11:02:38 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "HepPDT/defs.h"
#include "HepPDT/TableBuilder.hh"
#include "HepPDT/ParticleDataTable.hh"
#include "SimGeneral/HepPDTRecord/interface/PdtEntry.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH2.h"
#include "TMath.h"
#include "TString.h"
#include <algorithm>
#include <list>
#include <map>
#include <vector>
#include <set>
#include <sstream>
#include <stdlib.h> 
#include <string>

//
// class declaration
// migrate to:https://twiki.cern.ch/twiki/bin/viewauth/CMS/ThreadedDQM
//https://github.com/cirkovic/my-cmssw/blob/master/DQMOffline/RecoB/plugins/BTagPerformanceAnalyzerOnData.cc
//https://github.com/cms-tau-pog/HighPtTau_539/blob/master/PhysicsTools/PatAlgos/plugins/PATGenCandsFromSimTracksProducer.cc


class DemoAnalyzer : public DQMEDAnalyzer {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      struct histo1DFormating
      {
	  std::string title;
	  std::string name;
	  std::string labelx;
	  std::string labely;
	  std::vector<double>  rangex;
      };
  
      struct histo2DFormating
      {
          std::string title;
          std::string name;
          std::string labelx;
          std::string labely;
          std::vector<double> rangex;
          std::vector<double> rangey;
      };


   private:
  
      virtual void dqmBeginRun(edm::Run const& ,edm::EventSetup const& ) override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  
      // ----------member data ---------------------------
      template <typename T> 
      void bookHistosPerParticle(std::vector<T> & pdg, DQMStore::IBooker & ibooker);

  
      void bookEnergyLossesRelatedInfo2D( std::vector<MonitorElement*>& histos1D, DQMStore::IBooker & ibooker,
					  const std::string &subDet, const std::string &pdgName, int iHisto );
  
      void bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>& histos2D, DQMStore::IBooker & ibooker, 
					  const std::string &subDet, const std::string &pdgName, int iHisto );
  

      typedef std::list<std::vector<MonitorElement*> > matEffHistos;
      typedef std::map<std::string,  matEffHistos > matEffHistosToPdgIdsToSelect;
      matEffHistosToPdgIdsToSelect histos1DToPdgIdsToSelect_;
      matEffHistosToPdgIdsToSelect histos2DToPdgIdsToSelect_;
  

      void defaultHisto1DFormating(int iHisto, bool publishedParticleFilterWarning);
      void defaultHisto2DFormating(int iHisto, bool publishedParticleFilterWarnings);
      void fillHisto1DInfo(matEffHistosToPdgIdsToSelect::const_iterator iterator, int subdetPosBegin, int subdetPosEnd, int bin, float dE, float dx);
      void fillHisto2DInfo(matEffHistosToPdgIdsToSelect::const_iterator iterator, int subdetPosBegin, int subdetPosEnd, float momentumAtEntry, float dE, float dx);
      std::vector<int> findSubDet (int detId, int subdetId, const TrackerTopology* const tTopo);
      bool passesFilterSelection(edm::SimTrackContainer::const_iterator simTrack);
      void parseConfiguration(const edm::ParameterSet& iConfig);
      void parseHisto1DFormating(const std::vector<edm::ParameterSet>& histo1DFormatingVPSets);
      void parseHisto2DFormating(const std::vector<edm::ParameterSet>& histo2DFormatingVPSets);
      void publishParticleFilterWarnings();
      void publishHisto1DFormatingWarnings(const std::vector<edm::ParameterSet>& histo1DFormatingVPSets);
      void publishHisto2DFormatingWarnings(const std::vector<edm::ParameterSet>& histo2DFormatingVPSets);

      //SimHits
      std::vector<edm::InputTag> simHitsTag_;
  
      //SimTracks
      std::vector<edm::InputTag> trackerContainers_;

      //
      edm::ParameterSet particleSelection_;
      std::vector<PdtEntry> pdtEntries_;   
      std::vector<int> pdgIdsToSelect_;
      bool selectAntiParticle_;

      //
      std::vector<double> binsToProject_;
      bool usePBins; 
      bool useEBins;    
      edm::ESHandle<HepPDT::ParticleDataTable> pdgTable_;

      //
      edm::ParameterSet particleFilter_;
      double etaMin_;
      double etaMax_;
      double pTMin_;
      double pTMax_;
      double pMin_;
      double pMax_;
      double EMin_;
      double EMax_;
      std::vector<int> pdgIdsToFilter_;
      bool filterAntiParticle_;

      // 
      std::vector<histo1DFormating> configuredHistos1D_;
      std::vector<histo2DFormating> configuredHistos2D_;

      // Detectors
      std::map<std::string, int> subDetsToBinsToProject_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig):
    simHitsTag_(iConfig.getParameter<std::vector<edm::InputTag> >("SimHitTags"))
 
{
    parseConfiguration(iConfig);
}


DemoAnalyzer::~DemoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------


void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
#ifdef rrDEBUG
  std::cout << "Famos analysis" << std::endl;
#endif
  // get event and run number
#ifdef rrDEBUG
  int t_Run   = iEvent.id().run();
  int t_Event = iEvent.id().event();
  std::cout
    << " #################################### Run " << t_Run 
    << " Event "                                    << t_Event << " #################################### " 
    << std::endl;
#endif

    using namespace edm;

    //Retrieve tracker topology from geometry
    ESHandle< TrackerTopology > tTopoHandle;
    iSetup.get<IdealGeometryRecord>().get(tTopoHandle);
    const TrackerTopology* const tTopo = tTopoHandle.product();

    //Retrieve the Monte Carlo truth Particles(SimTracks)
    Handle< SimTrackContainer > simTracksHandle;
    iEvent.getByLabel(trackerContainers_[0].label(),simTracksHandle); 
    const SimTrackContainer simTracks = *(simTracksHandle.product());

    std::vector<double>::iterator bin;
    typedef std::map<std::string,  std::list<std::vector<MonitorElement*> > >::iterator histosToPdgIdsToSelect;

    for ( SimTrackContainer::const_iterator  simTrack = simTracks.begin(); simTrack != simTracks.end(); ++simTrack) 
    {
	for (unsigned i=0; i<trackerContainers_.size(); ++i) 
	{
           //Retrieve the low level info(SimHits)
	   Handle<std::vector<PSimHit> > simHitCollection;
	   iEvent.getByLabel(trackerContainers_[i], simHitCollection);
	   const std::vector<PSimHit>& simHits = *simHitCollection.product();
	   if(simHitCollection.isValid()) 
	   {
	       for( std::vector<PSimHit>::const_iterator Hits=simHits.begin(); Hits!=simHits.end(); ++Hits) 
	       {
		   int simHitID = (*Hits).trackId();
		   int simTrackID = (*simTrack).trackId();
		   int simHitparticleID = (*Hits).particleType();
		   if(simHitID == simTrackID && passesFilterSelection(simTrack))
		   {
		       std::set<unsigned int> detIds;
		       unsigned int detId = (*Hits).detUnitId();
		       unsigned int subdetId = DetId(detId).subdetId();
		       //Calculate the transversed thickness (in cm)
		       float dx = TMath::Power(
					       TMath::Power( (*Hits).entryPoint().x() - (*Hits).exitPoint().x(), 2) + 
					       TMath::Power( (*Hits).entryPoint().y() - (*Hits).exitPoint().y(), 2) + 
					       TMath::Power( (*Hits).entryPoint().z() - (*Hits).exitPoint().z(), 2), 
					       1/2.);
		       //Retrieve the experienced energy loss (in GeV)
		       float dE = (*Hits).energyLoss();
		       //Retrieve momentum/energy (in GeV, c=1)
		       float momentumAtEntry = 0;
		       if (usePBins)
		       {
			   momentumAtEntry = (*Hits).momentumAtEntry().mag();
			   binsToProject_.push_back(momentumAtEntry);
		       }
		       else
		       {
			   momentumAtEntry = (*Hits).momentumAtEntry().mag();
			   const HepPDT::ParticleData* PData = pdgTable_->particle(HepPDT::ParticleID(abs(simHitparticleID))) ;
			   float mass = PData->mass().value() ;
			   momentumAtEntry+= sqrt(pow(momentumAtEntry,2) + pow(mass,2));
			   binsToProject_.push_back(momentumAtEntry);
		       }
		       
		       sort(binsToProject_.begin(),binsToProject_.end());	
		       bin=find(binsToProject_.begin(),binsToProject_.end(), momentumAtEntry);

		       int binToProject = std::distance(binsToProject_.begin(), bin);
		       for( histosToPdgIdsToSelect iterator = histos2DToPdgIdsToSelect_.begin(); iterator != histos2DToPdgIdsToSelect_.end(); ++iterator )
		       {
			   std::vector<int> subdetPosBeginAndEnd = findSubDet(detId, subdetId, tTopo );
			   fillHisto2DInfo(iterator, subdetPosBeginAndEnd[0], subdetPosBeginAndEnd[1], momentumAtEntry, dE, dx);
		       }
 
		       for( histosToPdgIdsToSelect iterator = histos1DToPdgIdsToSelect_.begin(); iterator != histos1DToPdgIdsToSelect_.end(); ++iterator )
		       {
		           if (binToProject!=0 && binToProject<int(binsToProject_.size()-1))
			   {
			       std::vector<int> subdetPosBeginAndEnd = findSubDet(detId, subdetId, tTopo );
			       fillHisto1DInfo(iterator, subdetPosBeginAndEnd[0], subdetPosBeginAndEnd[1], binToProject, dE, dx);
			   }
			   binsToProject_.erase(std::remove(binsToProject_.begin(), binsToProject_.end(), momentumAtEntry), binsToProject_.end());
		       }
		   }
	       }
	   }
	}
    }
    
#ifdef THIS_IS_AN_EVENT_EXAMPLE
    Handle<ExampleData> pIn;
    iEvent.getByLabel("example",pIn);
#endif
    
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
#endif
}


void 
DemoAnalyzer::dqmBeginRun(edm::Run const& r, edm::EventSetup const& es)
{
  es.getData( pdgTable_ ) ;
  
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}


template <typename T> 
void 
DemoAnalyzer::bookHistosPerParticle(std::vector<T> & pdgNames, DQMStore::IBooker & ibooker)
{
    for (unsigned int iPart=0; iPart<pdgNames.size(); ++iPart) 
    {
        for ( auto& x: subDetsToBinsToProject_)
	{
	    std::vector<MonitorElement*> histos1D;
	    std::vector<MonitorElement*> histos2D;

	    std::string subDetName; 
	    subDetName.append(x.first);

	    int nHistos;
	    nHistos = x.second;

	    for(unsigned int iHisto = 0; iHisto < nHistos; iHisto++) 
	    {   
	        bookEnergyLossesRelatedInfo1D( histos1D, ibooker, subDetName, pdgNames[iPart], iHisto );
		histos1DToPdgIdsToSelect_[pdgNames[iPart]].push_back(histos1D);
		if (iHisto == 0)
		{
		    bookEnergyLossesRelatedInfo2D( histos2D, ibooker, subDetName, pdgNames[iPart], iHisto );
		    histos2DToPdgIdsToSelect_[pdgNames[iPart]].push_back(histos2D);
		}
	    }
	}
    }

}


void 
DemoAnalyzer::bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>& histos1D, DQMStore::IBooker & ibooker, 
					     const std::string &subDet, const std::string &pdgName, int iHisto )
{
    std::string name(configuredHistos1D_[iHisto].name.data());
    std::string title(configuredHistos1D_[iHisto].title.data());
    std::string labelx(configuredHistos1D_[iHisto].labelx.data());
    std::string labely(configuredHistos1D_[iHisto].labely.data());
    int nBins = configuredHistos1D_[iHisto].rangex[0];
    double lowBinX = configuredHistos1D_[iHisto].rangex[1];
    double highBinX = configuredHistos1D_[iHisto].rangex[2];
    std::string gevOverC = usePBins ? "GeV/c" : "GeV/c^{2}";
    if (name.compare("default")==0 && title.compare("default")==0) 
    {
         histos1D.push_back( ibooker.book1D( Form( "SimHit_%s_%s_%s_%u", labelx.data(), subDet.data(), pdgName.data(), iHisto+1 ) ,
					     Form( "(%.2f,%.2f) %s;%s;%s", binsToProject_[iHisto] , binsToProject_[iHisto+1], gevOverC.data(),
						   labelx.data(), labely.data() ) , 
					     nBins , lowBinX, highBinX ) );
    }
    else
    {
	 histos1D.push_back( ibooker.book1D( Form( "SimHit_%s_%s_%s_%u", name.data(), subDet.data(),  pdgName.data(), iHisto+1), 
					     Form( "(%.2f,%.2f) %s_%s;%s;%s", binsToProject_[iHisto] , binsToProject_[iHisto+1], gevOverC.data(),
						   title.data(), labelx.data(), labely.data() ) ,
					     nBins , lowBinX, highBinX ) );
    }

}


void 
DemoAnalyzer::bookEnergyLossesRelatedInfo2D( std::vector<MonitorElement*>& histos2D, DQMStore::IBooker & ibooker,
                                             const std::string &subDet, const std::string &pdgName, int iHisto)

{
    std::string name(configuredHistos2D_[iHisto].name.data());
    std::string title(configuredHistos2D_[iHisto].title.data());
    std::string labelx(configuredHistos2D_[iHisto].labelx.data());
    std::string labely(configuredHistos2D_[iHisto].labely.data());
    int nBinsX = configuredHistos2D_[iHisto].rangex[0];
    double lowBinX = configuredHistos2D_[iHisto].rangex[1];
    double highBinX = configuredHistos2D_[iHisto].rangex[2];
    int nBinsY = configuredHistos2D_[iHisto].rangey[0];
    double lowBinY = configuredHistos2D_[iHisto].rangey[1];
    double highBinY = configuredHistos2D_[iHisto].rangey[2];
    std::string gevOverC = usePBins ? "GeV/c" : "GeV/c^{2}";
    if (name.compare("default")==0 && title.compare("default")==0)
    {
        histos2D.push_back( ibooker.book2D( Form( "SimHit_%s_vs_%s_%s_%s_%u", labelx.data(), labely.data(), subDet.data(), pdgName.data(), iHisto+1) ,
					    Form( ";%s (%s);%s", labelx.data(), gevOverC.data(), labely.data() ) ,
					    nBinsX , lowBinX, highBinX,  nBinsY , lowBinY, highBinY) );
    }
    else
    {
        histos2D.push_back( ibooker.book2D( Form( "SimHit_%s_%s_%s", name.data(), subDet.data(),  pdgName.data() ),
					    Form( "%s;%s;%s", title.data(), labelx.data(), labely.data() ) ,
					    nBinsX , lowBinX, highBinX, nBinsY , lowBinY, highBinY ) );
    }

}


void DemoAnalyzer::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const & iRun,
				  edm::EventSetup const & iSetup)
{
     ibooker.cd();
     ibooker.setCurrentFolder("testMaterialEffects");
     
     std::vector<std::string> pdgNames; 
     if (!pdtEntries_.empty()) 
     {
         for (std::vector<PdtEntry>::iterator itp = pdtEntries_.begin(), edp = pdtEntries_.end(); itp != edp; ++itp) 
	 {
	     itp->setup(iSetup); // decode string->pdgId and vice-versa
	     pdgNames.push_back(itp->name());
	     pdgIdsToSelect_.push_back(std::abs(itp->pdgId()));
	 }
     }
     else
     {
         pdgNames.push_back("allParticles");
     }
     
     sort(pdgIdsToSelect_.begin(),pdgIdsToSelect_.end());
     publishParticleFilterWarnings();
     bookHistosPerParticle(pdgNames, ibooker);
	  
}


void
DemoAnalyzer::defaultHisto1DFormating(int iHisto, bool publishedParticleFilterWarnings=false)
{
    //========
    //dynamic title and name (will be substistuted in bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>& , DQMStore::IBooker & , const std::string &, const std::string &, int ) ) 
    //========
    if(!publishedParticleFilterWarnings)
    {
	configuredHistos1D_[iHisto].title = "default"; 
	configuredHistos1D_[iHisto].name = "default";
	configuredHistos1D_[iHisto].labelx = "dEdx";
	configuredHistos1D_[iHisto].labely = "";
    }
    configuredHistos1D_[iHisto].rangex.push_back(200); //nBins
    configuredHistos1D_[iHisto].rangex.push_back(0.); //xMin
    configuredHistos1D_[iHisto].rangex.push_back(0.01); //xMax
        
}


void
DemoAnalyzer::defaultHisto2DFormating(int iHisto, bool publishedParticleFilterWarnings=false)
{
    //========
    //dynamic title and name (will be substistuted in bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>& , DQMStore::IBooker  & , const std::string &, const std::string &, int ) )
    //========
    if(!publishedParticleFilterWarnings)
    {
        configuredHistos2D_[iHisto].title= "default";
	configuredHistos2D_[iHisto].name = "default";
	configuredHistos2D_[iHisto].labelx = "p";
	configuredHistos2D_[iHisto].labely = "dEdx";
    }
    configuredHistos2D_[iHisto].rangex.push_back(50); //nBinsX
    configuredHistos2D_[iHisto].rangex.push_back(0.); //xMin
    configuredHistos2D_[iHisto].rangex.push_back(60); //xMax
    configuredHistos2D_[iHisto].rangey.push_back(200); //nBinsY
    configuredHistos2D_[iHisto].rangey.push_back(0.); //yMin
    configuredHistos2D_[iHisto].rangey.push_back(0.01); //yMax
}


void 
DemoAnalyzer::fillHisto1DInfo(matEffHistosToPdgIdsToSelect::const_iterator iterator,  int subdetPosBegin, int subdetPosEnd, 
			      int binToProject, float dE, float dx)
{
    matEffHistos histos = iterator->second;
    int detct = 0;
    for (std::list<std::vector<MonitorElement*> >::iterator it1 = histos.begin(); it1 != histos.end(); ++it1)
    {
      int count =0 ;
      std::vector<MonitorElement*>::iterator it2;
      for (it2 = (*it1).begin(); it2 != (*it1).end(); ++ it2)
      {
	  if ((binToProject-1) == count && (detct/int(binsToProject_.size()-2))>=subdetPosBegin && (detct/int(binsToProject_.size()-2))<subdetPosEnd)
	  {
	      (*it2)->Fill(dE/dx);
	  }
	  count++;
      }
      detct++;
    }
}


void
DemoAnalyzer::fillHisto2DInfo(matEffHistosToPdgIdsToSelect::const_iterator iterator,  int subdetPosBegin, int subdetPosEnd, 
			      float momentumAtEntry, float dE, float dx)
{

  matEffHistos histos = iterator->second;
  int detct = 0;
  for (std::list<std::vector<MonitorElement*> >::iterator it1 = histos.begin(); it1 != histos.end(); ++it1)
  {
      int count =0 ;
      std::vector<MonitorElement*>::iterator it2;
      for (it2 = (*it1).begin(); it2 != (*it1).end(); ++ it2)
      {
          if (detct/1>=subdetPosBegin && detct/1<subdetPosEnd)
	  {
	      (*it2)->Fill(momentumAtEntry, dE/dx);
	  }
          count++;
      }
      detct++;
  }
}


std::vector<int> 
DemoAnalyzer::findSubDet (int detId, int subdetId, const TrackerTopology* const tTopo)
{
    std::vector<int> subdetPosBeginAndEnd; 
    if (subdetId == static_cast<int>(PixelSubdetector::PixelBarrel) || subdetId == static_cast<int>(PixelSubdetector::PixelEndcap)) 
    {
	subdetPosBeginAndEnd.push_back(0);
	subdetPosBeginAndEnd.push_back(1);
    }
    else if (subdetId == StripSubdetector::TIB || subdetId == StripSubdetector::TOB || 
	     subdetId == StripSubdetector::TID || subdetId == StripSubdetector::TEC )
    {
	if (subdetId == StripSubdetector::TEC && tTopo->tecRing(detId)>=5) 
	{
	    subdetPosBeginAndEnd.push_back(2);
	    subdetPosBeginAndEnd.push_back(3);
	}
	else 
	{
	    subdetPosBeginAndEnd.push_back(1);
	    subdetPosBeginAndEnd.push_back(2);	    
	}
    }
    
    return subdetPosBeginAndEnd;
}    


bool
DemoAnalyzer::passesFilterSelection(edm::SimTrackContainer::const_iterator simTrack)
{

   if ((*simTrack).momentum().Eta() < etaMin_)
   {
       return false;
   }
   if ((*simTrack).momentum().Eta() > etaMax_)
   {
       return false;
   }
   if ((*simTrack).momentum().Pt() < pTMin_)
   {
       return false;
   }
   if ((*simTrack).momentum().Pt() > pTMax_)
   {
       return false;
   }
   if ((*simTrack).momentum().P() < pMin_)
   {
       return false;
   }
   if ((*simTrack).momentum().E() < EMin_)
   {
       return false;
   }
   int particleIDToFilter = filterAntiParticle_ ? (*simTrack).type() : fabs((*simTrack).type());
   if (std::find(pdgIdsToFilter_.begin(), pdgIdsToFilter_.end(), particleIDToFilter) != pdgIdsToFilter_.end())
   {
       return false;
   }
   if (!pdtEntries_.empty()) 
   {
       int particleIDToSelect = selectAntiParticle_ ? (*simTrack).type() : fabs((*simTrack).type());
       if (std::find(pdgIdsToSelect_.begin(), pdgIdsToSelect_.end(), particleIDToSelect) != pdgIdsToSelect_.end())
       {
	   return true; 
       }
       else
       {
	   return false;
       }
   }

   return true;

}


void
DemoAnalyzer::parseConfiguration(const edm::ParameterSet& iConfig)
{
    for (unsigned int i_subDet=0; i_subDet<simHitsTag_.size(); ++i_subDet)
    {
	trackerContainers_.push_back(simHitsTag_[i_subDet]);
    }

    // Possibly allow a list of particle types
    if (iConfig.exists("particleSelection"))
    {
    	particleSelection_ = iConfig.getParameter<edm::ParameterSet>("particleSelection");
	pdtEntries_ = particleSelection_.getUntrackedParameter<std::vector<PdtEntry> >("particleTypes");
	selectAntiParticle_ =  particleSelection_.getUntrackedParameter<bool>("selectAntiParticle");
    }
    
    if (iConfig.exists("bins_p") && iConfig.exists("bins_E"))
    {	
        edm::LogWarning("TooManyBins") << "Only momentum bins will be used";
	binsToProject_ = iConfig.getUntrackedParameter <std::vector <double> >("bins_p");
	usePBins = true;useEBins = false;
    }
    else if (iConfig.exists("bins_p") && !iConfig.exists("bins_E"))
    {
	binsToProject_ = iConfig.getUntrackedParameter <std::vector <double> >("bins_p");
	usePBins = true;useEBins = false;
    }
    else if (!iConfig.exists("bins_p") && iConfig.exists("bins_E"))
    {
        binsToProject_ = iConfig.getUntrackedParameter <std::vector <double> >("bins_E");
	usePBins = false;useEBins = true;
    }


    sort(binsToProject_.begin(),binsToProject_.end());
    
    if (iConfig.exists("particleFilter"))
    {
    	particleFilter_ = iConfig.getParameter<edm::ParameterSet>("particleFilter");
	etaMin_ = particleFilter_.getUntrackedParameter<double>("etaMin");
        etaMax_ = particleFilter_.getUntrackedParameter<double>("etaMax");
	pTMin_ = particleFilter_.getUntrackedParameter<double>("pTMin");
	pTMax_ = particleFilter_.getUntrackedParameter<double>("pTMax");
	pMin_ = particleFilter_.getUntrackedParameter<double>("pMin");
	pMax_ = particleFilter_.getUntrackedParameter<double>("pMax");
	EMin_ = particleFilter_.getUntrackedParameter<double>("EMin");
	EMax_ = particleFilter_.getUntrackedParameter<double>("EMax");
	pdgIdsToFilter_ = particleFilter_.getUntrackedParameter<std::vector<int> >("pdgIdsToFilter");
	filterAntiParticle_ =  particleFilter_.getUntrackedParameter<bool>("filterAntiParticle");
	sort(pdgIdsToFilter_.begin(),pdgIdsToFilter_.end());
    }

    
    subDetsToBinsToProject_.insert(std::pair<std::string,int>("Pixels",binsToProject_.size()-1));
    subDetsToBinsToProject_.insert(std::pair<std::string,int>("StripsNoTEC5to7",binsToProject_.size()-1));
    subDetsToBinsToProject_.insert(std::pair<std::string,int>("StripsOnlyTEC5to7",binsToProject_.size()-1));
    

    if (iConfig.exists("histo1DFormating"))
    {
        const std::vector<edm::ParameterSet>& histo1DFormatingVPSets = iConfig.getParameter<std::vector<edm::ParameterSet>>("histo1DFormating");
	parseHisto1DFormating(histo1DFormatingVPSets);
    }
    if (iConfig.exists("histo2DFormating"))
    {
        const std::vector<edm::ParameterSet>& histo2DFormatingVPSets = iConfig.getParameter<std::vector<edm::ParameterSet>>("histo2DFormating");
	parseHisto2DFormating(histo2DFormatingVPSets);
    }

}


void
DemoAnalyzer::parseHisto1DFormating(const std::vector<edm::ParameterSet>& histo1DFormatingVPSets)
{
    for (unsigned int iHisto=0; iHisto< binsToProject_.size()-1; ++iHisto)
    {   
        configuredHistos1D_.push_back(histo1DFormating());
        if(iHisto<histo1DFormatingVPSets.size()) //Read parameters given in the configuration file 
	{
	    configuredHistos1D_[iHisto].title = histo1DFormatingVPSets[iHisto].getParameter<std::string>("title");
	    configuredHistos1D_[iHisto].name = histo1DFormatingVPSets[iHisto].getParameter<std::string>("name");
	    configuredHistos1D_[iHisto].labelx = histo1DFormatingVPSets[iHisto].getUntrackedParameter<std::string>("labelx");
	    configuredHistos1D_[iHisto].labely = histo1DFormatingVPSets[iHisto].getUntrackedParameter<std::string>("labely");
	    configuredHistos1D_[iHisto].rangex = histo1DFormatingVPSets[iHisto].getUntrackedParameter<std::vector<double> >("rangex");
	}
	else //Default parameters 
	{
	    defaultHisto1DFormating(iHisto);
	}
	    
    }
    publishHisto1DFormatingWarnings(histo1DFormatingVPSets);
    
}


void
DemoAnalyzer::parseHisto2DFormating(const std::vector<edm::ParameterSet>& histo2DFormatingVPSets)
{
    for (unsigned int iHisto=0; iHisto< subDetsToBinsToProject_.size(); ++iHisto)
    {
        configuredHistos2D_.push_back(histo2DFormating());
        if(iHisto<histo2DFormatingVPSets.size()) //Read parameters given in the configuration file
        {
	    configuredHistos2D_[iHisto].title = histo2DFormatingVPSets[iHisto].getParameter<std::string>("title");
	    configuredHistos2D_[iHisto].name = histo2DFormatingVPSets[iHisto].getParameter<std::string>("name");
	    configuredHistos2D_[iHisto].labelx = histo2DFormatingVPSets[iHisto].getUntrackedParameter<std::string>("labelx");
	    configuredHistos2D_[iHisto].labely = histo2DFormatingVPSets[iHisto].getUntrackedParameter<std::string>("labely");
	    configuredHistos2D_[iHisto].rangex = histo2DFormatingVPSets[iHisto].getUntrackedParameter<std::vector<double> >("rangex");
	    configuredHistos2D_[iHisto].rangey = histo2DFormatingVPSets[iHisto].getUntrackedParameter<std::vector<double> >("rangey");
	}
	else //Default parameters
	{
	    defaultHisto2DFormating(iHisto);
	}
    }
    publishHisto2DFormatingWarnings(histo2DFormatingVPSets);
    
}


void
DemoAnalyzer::publishHisto1DFormatingWarnings(const std::vector<edm::ParameterSet>&  histo1DFormatingVPSets)
{
    if(histo1DFormatingVPSets.size()>binsToProject_.size()+1)
    {
	edm::LogWarning("TooMany1DHistograms") << "Number of reserved 1D Histos will be reduced to number of binsToProject minus one";
    }
    for (unsigned int iHisto=0; iHisto< binsToProject_.size()-1; ++iHisto)
    {  
        if(iHisto<histo1DFormatingVPSets.size())
	{
	    if (configuredHistos1D_[iHisto].rangex.size()!=3)
	    {
		edm::LogWarning("NotProperBinConfiguration1DHistograms") << "Wrong bin configuration in the reserved 1D Histos; will be reduced to the default ones";
		//Then, use default parameters only for the bins
		defaultHisto1DFormating(iHisto, true);
	    }
	}
    }
}


void
DemoAnalyzer::publishHisto2DFormatingWarnings(const std::vector<edm::ParameterSet>&  histo2DFormatingVPSets)
{
    if(histo2DFormatingVPSets.size()>subDetsToBinsToProject_.size())   
    {
        edm::LogWarning("TooMany2DHistograms") << "Number of reserved 2D Histos will be reduced to number of reserved sub-detectors";
    }
    for (unsigned int iHisto=0; iHisto< subDetsToBinsToProject_.size(); ++iHisto)
    {
        if(iHisto<histo2DFormatingVPSets.size()) 
        {
	    if (configuredHistos2D_[iHisto].rangex.size()!=3 || configuredHistos2D_[iHisto].rangey.size()!=3)
	    {
		edm::LogWarning("NotProperBinConfiguration2DHistograms") << "Wrong bin configuration in the reserved 2D Histos; both x- and y-axes will be reduced to the default ones";
		//Then, use default parameters only for the bins
		defaultHisto2DFormating(iHisto, true);
	    }
	}
    }
}


void
DemoAnalyzer::publishParticleFilterWarnings()
{
    if (etaMin_== etaMax_)
    {
        edm::LogWarning("IdenticalFilterValues") << "Unique (min=max) pseudo-rapidity filter value will be used";
    }
    if (pTMin_== pTMax_)
    {
	if (pTMin_< std::numeric_limits<double>::epsilon())
	{
	    edm::LogWarning("DisabledFilterValues") << "No pT filter value will be applied";
	    pTMin_ = -1.;
	    pTMax_ = std::numeric_limits<double>::max();
	}
	else
	{
	    edm::LogWarning("IdenticalFilterValues") << "Unique (min=max) pT filter value will be used with machine tolerance";
	    pTMin_ = pTMin_ - std::numeric_limits<double>::epsilon();
	    pTMax_ = pTMax_ + std::numeric_limits<double>::epsilon();
	}
    }
    if (pMin_== pMax_)
    {
        if (pMin_< std::numeric_limits<double>::epsilon())
	{
	    edm::LogWarning("DisabledFilterValues") << "No p (momentum) filter value will be applied";
	    pMin_ = -1.;
            pMax_ = std::numeric_limits<double>::max();
	}
        else
	{
	    edm::LogWarning("IdenticalFilterValues") << "Unique p (momentum) filter value will be used";
	    pMin_ = pMin_ - std::numeric_limits<double>::epsilon();
            pMax_ = pMax_ + std::numeric_limits<double>::epsilon();
	}
    }
    if (EMin_== EMax_)
    {
        if (EMin_< std::numeric_limits<double>::epsilon())
        {
	    edm::LogWarning("DisabledFilterValues") << "No E (energy) filter value will be applied";
	    EMin_ = -1.;
            EMax_ = std::numeric_limits<double>::max();
	}
        else
        {
	    edm::LogWarning("IdenticalFilterValues") << "Unique E (energy) filter value will be used";
	    EMin_ = pMin_ - std::numeric_limits<double>::epsilon();
            EMax_ = pMax_ + std::numeric_limits<double>::epsilon();
	}
    }
    
    std::set<int> intersect;
    std::set<int>::iterator it;
    std::set_intersection(pdgIdsToSelect_.begin(), pdgIdsToSelect_.end(), pdgIdsToFilter_.begin(), pdgIdsToFilter_.end(), std::inserter(intersect,intersect.begin()));
    if (intersect.size()!=0)
    {
	edm::LogWarning("IdenticalToFilterAndToSelectPDGValues") << "No PDGfilter can be applied for identical PDGselect values!";
	for (it=intersect.begin(); it!=intersect.end(); ++it)
	{
	    pdgIdsToFilter_.erase(std::remove(pdgIdsToFilter_.begin(), pdgIdsToFilter_.end(), *it), pdgIdsToFilter_.end());
	}
    }
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
