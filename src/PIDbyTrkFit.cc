#include "PIDbyTrkFit.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCObject.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCParameters.h>
#include <marlinutil/GeometryUtil.h>
#include <UTIL/LCRelationNavigator.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

using namespace lcio ;
using namespace marlin ;

MyPIDbyTrkFit aMyPIDbyTrkFit;

MyPIDbyTrkFit::MyPIDbyTrkFit() :
	Processor("PIDbyTrkFit"),
	m_nRun(0),
	m_nEvt(0),
	m_nRunSum(0),
	m_nEvtSum(0),

	m_trk_pdg{},
	m_trk_charge{},
	m_trk_genStat{},
	m_highest_weight{},
	m_ndf{},
	m_chi2{},
	m_chi2_ndf{},
	m_ndf_electron{},
	m_chi2_electron{},
	m_chi2_ndf_electron{},
	m_ndf_muon{},
	m_chi2_muon{},
	m_chi2_ndf_muon{},
	m_ndf_pion{},
	m_chi2_pion{},
	m_chi2_ndf_pion{},
	m_ndf_kaon{},
	m_chi2_kaon{},
	m_chi2_ndf_kaon{},
	m_ndf_proton{},
	m_chi2_proton{},
	m_chi2_ndf_proton{},

	m_pTFile(NULL),
	m_pTTree(NULL)

{
	// modify processor description
        _description = "MyPIDbyTrkFit finds Track parameters regarding different mass hypothesis" ;


        // register steering parameters: name, description, class-variable, default value
            registerInputCollection(LCIO::MCPARTICLE,
     				       "MCParticleCollection" ,
     				       "Name of the MCParticle collection"  ,
     				       m_mcParticleCollection,
     				       std::string("MCParticle")
     				       );

            registerInputCollection(LCIO::LCRELATION,
     				       "MarlinTrkTracksMCTruthLink" ,
     				       "Name of the Track-MCParticle Relations collection for input tracks"  ,
     				       m_TrackMCParticleRelCol ,
     				       std::string("MarlinTrkTracksMCTruthLink")
     				       );

            registerInputCollection(LCIO::LCRELATION,
     				       "MCTruthMarlinTrkTracksLink" ,
     				       "Name of the MCParticle-Track Relations collection for input MCParticles"  ,
     				       m_MCParticleTrackRelCol ,
     				       std::string("MCTruthMarlinTrkTracksLink")
     				       );

            registerOutputCollection(LCIO::MCPARTICLE,
     				       "MCParticleCollectionVector" ,
     				       "Name of the MCParticle collection for Kaons"  ,
     				       m_col_mcp,
     				       std::string("MCParticleCollectionVector")
     				       );

            registerInputCollection( LCIO::TRACK,
     				       "MarlinTrkTracksCollection" ,
     				       "Name of the MarlinTrkTracks collection"  ,
     				       m_MarlinTrkTracks,
     				       std::string("MarlinTrkTracks")
     				       );

            registerInputCollection( LCIO::TRACK,
     				       "MarlinTrkTracksCollectionPion" ,
     				       "Name of the MarlinTrkTracks collection"  ,
     				       m_MarlinTrkTracksPION,
     				       std::string("MarlinTrkTracksPion")
     				       );

            registerInputCollection( LCIO::TRACK,
     				       "MarlinTrkTracksCollectionMuon" ,
     				       "Name of the MarlinTrkTracks collection"  ,
     				       m_MarlinTrkTracksMUON,
     				       std::string("MarlinTrkTracksMuon")
     				       );

            registerInputCollection( LCIO::TRACK,
     				       "MarlinTrkTracksCollectionKaon" ,
     				       "Name of the MarlinTrkTracks collection"  ,
     				       m_MarlinTrkTracksKAON,
     				       std::string("MarlinTrkTracksKaon")
     				       );

            registerInputCollection( LCIO::TRACK,
     				       "MarlinTrkTracksCollectionProton" ,
     				       "Name of the MarlinTrkTracks collection"  ,
     				       m_MarlinTrkTracksPROTON,
     				       std::string("MarlinTrkTracksProton")
     				       );

            registerInputCollection( LCIO::TRACK,
     				       "MarlinTrkTracksCollectionElectron" ,
     				       "Name of the MarlinTrkTracks collection"  ,
     				       m_MarlinTrkTracksELECTRON,
     				       std::string("MarlinTrkTracksElectron")
     				       );

            registerProcessorParameter("RootFile",
     				       "Name of the output root file",
     				       m_rootFile,
     				       std::string("TrkFitGOF.root")
     				       );
}

void MyPIDbyTrkFit::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    // usually a good idea to
    printParameters() ;

	m_nRun = 0 ;
	m_nEvt = 0 ;
	m_nRunSum = 0;
	m_nEvtSum = 0;
	this->Clear();

	m_pTFile = new TFile(m_rootFile.c_str(), "recreate");
	m_pTTree = new TTree("TrackParameters", "TrackParameters");
	m_pTTree->SetDirectory(0);
	m_pTTree->Branch("run", &m_nRun, "run/I");
	m_pTTree->Branch("event", &m_nEvt, "event/I");

	m_pTTree->Branch("track_PDG", &m_trk_pdg);
	m_pTTree->Branch("trk_GeneratorStatus", &m_trk_genStat);
	m_pTTree->Branch("highest_weight", &m_highest_weight);
	m_pTTree->Branch("track_charge", &m_trk_charge);
	m_pTTree->Branch("NDF", &m_ndf);
	m_pTTree->Branch("CHI2", &m_chi2);
	m_pTTree->Branch("CHI2_NDF", &m_chi2_ndf);
	m_pTTree->Branch("NDF_electron", &m_ndf_electron);
	m_pTTree->Branch("CHI2_electron", &m_chi2_electron);
	m_pTTree->Branch("CHI2_NDF_electron", &m_chi2_ndf_electron);
	m_pTTree->Branch("NDF_muon", &m_ndf_muon);
	m_pTTree->Branch("CHI2_muon", &m_chi2_muon);
	m_pTTree->Branch("CHI2_NDF_muon", &m_chi2_ndf_muon);
	m_pTTree->Branch("NDF_pion", &m_ndf_pion);
	m_pTTree->Branch("CHI2_pion", &m_chi2_pion);
	m_pTTree->Branch("CHI2_NDF_pion", &m_chi2_ndf_pion);
	m_pTTree->Branch("NDF_kaon", &m_ndf_kaon);
	m_pTTree->Branch("CHI2_kaon", &m_chi2_kaon);
	m_pTTree->Branch("CHI2_NDF_kaon", &m_chi2_ndf_kaon);
	m_pTTree->Branch("NDF_proton", &m_ndf_proton);
	m_pTTree->Branch("CHI2_proton", &m_chi2_proton);
	m_pTTree->Branch("CHI2_NDF_proton", &m_chi2_ndf_proton);

}
void MyPIDbyTrkFit::processRunHeader( LCRunHeader *pLCRunHeader)
{
	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;
}



void MyPIDbyTrkFit::processEvent( LCEvent *pLCEvent )
{

	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();

	++m_nEvtSum;

	if ((m_nEvtSum % 100) == 0)
         std::cout << " processed events: " << m_nEvtSum << std::endl;

    // New reco particle collection for the semi-leptonic decays

	this->Clear();
	this->ExtractCollections(pLCEvent);
	this->CheckGOF(pLCEvent);
	m_pTTree->Fill();
}



void MyPIDbyTrkFit::check( LCEvent *pLCEvent )
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void MyPIDbyTrkFit::end()
{

	m_pTFile->cd();
	m_pTTree->Write();

	m_pTFile->Close();
	delete m_pTFile;

    //   std::cout << "MySLDecayFinder::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;

}

void MyPIDbyTrkFit::Clear()
{

	m_trk_pdg.clear();
	m_trk_charge.clear();
	m_trk_genStat.clear();
	m_highest_weight.clear();
	m_ndf.clear();
	m_chi2.clear();
	m_chi2_ndf.clear();
	m_ndf_electron.clear();
	m_chi2_electron.clear();
	m_chi2_ndf_electron.clear();
	m_ndf_muon.clear();
	m_chi2_muon.clear();
	m_chi2_ndf_muon.clear();
	m_ndf_pion.clear();
	m_chi2_pion.clear();
	m_chi2_ndf_pion.clear();
	m_ndf_kaon.clear();
	m_chi2_kaon.clear();
	m_chi2_ndf_kaon.clear();
	m_ndf_proton.clear();
	m_chi2_proton.clear();
	m_chi2_ndf_proton.clear();

}

void MyPIDbyTrkFit::ExtractCollections(EVENT::LCEvent *pLCEvent)
{
	try
     {
     	const EVENT::LCCollection *pLCCollectionPion = pLCEvent->getCollection(m_MarlinTrkTracksPION);

     	for (unsigned int i = 0, nElements = pLCCollectionPion->getNumberOfElements(); i < nElements; ++i)
		{
			Track *pionTrack = dynamic_cast<EVENT::Track*>(pLCCollectionPion->getElementAt(i));

			if (NULL == pionTrack)
			throw EVENT::Exception("Collection type mismatch");

//			if (!pMCParticle->getParents().empty())
//				continue;

		}
	}
     catch (...)
     {
		streamlog_out(WARNING) << "Could not extract track collection for using pion mass " << std::endl;
     }
}


void MyPIDbyTrkFit::CheckGOF(EVENT::LCEvent *pLCEvent)
{
	try
	{
		const EVENT::LCCollection *mcpCollection = pLCEvent->getCollection(m_mcParticleCollection);
		const EVENT::LCCollection *trkCollection = pLCEvent->getCollection(m_MarlinTrkTracks);
		const EVENT::LCCollection *trkCollectionPion = pLCEvent->getCollection(m_MarlinTrkTracksPION);
		const EVENT::LCCollection *trkCollectionKaon = pLCEvent->getCollection(m_MarlinTrkTracksKAON);
		const EVENT::LCCollection *trkCollectionMuon = pLCEvent->getCollection(m_MarlinTrkTracksMUON);
		const EVENT::LCCollection *trkCollectionElectron = pLCEvent->getCollection(m_MarlinTrkTracksELECTRON);
		const EVENT::LCCollection *trkCollectionProton = pLCEvent->getCollection(m_MarlinTrkTracksPROTON);
		const EVENT::LCCollection *relCollection = pLCEvent->getCollection(m_TrackMCParticleRelCol);

		LCRelationNavigator track2mcNav(pLCEvent->getCollection( m_TrackMCParticleRelCol ));
//		streamlog_out(DEBUG) << " got track2mcNav from " << track2mcNav.getFromType() << " to " << track2mcNav.getToType() << std::endl;

		unsigned int ntrks = trkCollection->getNumberOfElements();
		unsigned int ntrksPion = trkCollectionPion->getNumberOfElements();
		unsigned int ntrksKaon = trkCollectionKaon->getNumberOfElements();
		unsigned int ntrksMuon = trkCollectionMuon->getNumberOfElements();
		unsigned int ntrksElectron = trkCollectionElectron->getNumberOfElements();
		unsigned int ntrksProton = trkCollectionProton->getNumberOfElements();

		if (ntrksPion != ntrks)
		{
			streamlog_out(DEBUG) << " <<<<<=====  Number of Pion tracks mismatch with the number of original track  =====>>>>>" << std::endl;
			streamlog_out(DEBUG) << " Number of original track = " << ntrks << std::endl;
			streamlog_out(DEBUG) << " Number of Pion track = " << ntrksPion << std::endl;
		}
		if (ntrksKaon != ntrks)
		{
			streamlog_out(DEBUG) << " <<<<<=====  Number of Kaon tracks mismatch with the number of original track  =====>>>>>" << std::endl;
			streamlog_out(DEBUG) << " Number of original track = " << ntrks << std::endl;
			streamlog_out(DEBUG) << " Number of Kaon track = " << ntrksKaon << std::endl;
		}
		if (ntrksMuon != ntrks)
		{
			streamlog_out(DEBUG) << " <<<<<=====  Number of Muon tracks mismatch with the number of original track  =====>>>>>" << std::endl;
			streamlog_out(DEBUG) << " Number of original track = " << ntrks << std::endl;
			streamlog_out(DEBUG) << " Number of Muon track = " << ntrksMuon << std::endl;
		}
		if (ntrksElectron != ntrks)
		{
			streamlog_out(DEBUG) << " <<<<<=====  Number of Electron tracks mismatch with the number of original track  =====>>>>>" << std::endl;
			streamlog_out(DEBUG) << " Number of original track = " << ntrks << std::endl;
			streamlog_out(DEBUG) << " Number of Electron track = " << ntrksElectron << std::endl;
		}
		if (ntrksProton != ntrks)
		{
			streamlog_out(DEBUG) << " <<<<<=====  Number of Proton tracks mismatch with the number of original track  =====>>>>>" << std::endl;
			streamlog_out(DEBUG) << " Number of original track = " << ntrks << std::endl;
			streamlog_out(DEBUG) << " Number of Proton track = " << ntrksProton << std::endl;
		}

		for (unsigned int i = 0; i < ntrks; ++i)
		{
			auto origTrack = dynamic_cast<EVENT::Track*>(trkCollection->getElementAt(i));
			auto electronTrack = dynamic_cast<EVENT::Track*>(trkCollectionElectron->getElementAt(i));
			auto muonTrack = dynamic_cast<EVENT::Track*>(trkCollectionMuon->getElementAt(i));
			auto pionTrack = dynamic_cast<EVENT::Track*>(trkCollectionPion->getElementAt(i));
			auto kaonTrack = dynamic_cast<EVENT::Track*>(trkCollectionKaon->getElementAt(i));
			auto protonTrack = dynamic_cast<EVENT::Track*>(trkCollectionProton->getElementAt(i));

			const EVENT::LCObjectVec& mcpvec = track2mcNav.getRelatedToObjects(origTrack);
			const EVENT::FloatVec&  trackweightvec = track2mcNav.getRelatedToWeights(origTrack);

			double maxweight = 0.;
			int imcpmax = 0;
			for ( unsigned int imcp = 0; imcp < mcpvec.size(); imcp++ )
			{
				if ( trackweightvec.at(imcp) > maxweight )
				{
					maxweight = trackweightvec.at(imcp);
					imcpmax = imcp;
				}
			}
			m_highest_weight.push_back(maxweight);
			MCParticle *mcpLinked = (MCParticle *) mcpvec.at(imcpmax);
			m_trk_genStat.push_back(mcpLinked->getGeneratorStatus());
			int trk_PDG = mcpLinked->getPDG();
			m_trk_pdg.push_back(trk_PDG);
			float trk_charge = mcpLinked->getCharge();
			m_trk_charge.push_back(trk_charge);

			m_chi2.push_back(origTrack->getChi2());
			m_ndf.push_back(origTrack->getNdf());
			if (origTrack->getNdf()!=0)
			{
				m_chi2_ndf.push_back(origTrack->getChi2()/origTrack->getNdf());
			}

			m_chi2_electron.push_back(electronTrack->getChi2());
			m_ndf_electron.push_back(electronTrack->getNdf());
			if (electronTrack->getNdf()!=0)
			{
				m_chi2_ndf_electron.push_back(electronTrack->getChi2()/electronTrack->getNdf());
			}

			m_chi2_muon.push_back(muonTrack->getChi2());
			m_ndf_muon.push_back(muonTrack->getNdf());
			if (muonTrack->getNdf()!=0)
			{
				m_chi2_ndf_muon.push_back(muonTrack->getChi2()/muonTrack->getNdf());
			}

			m_chi2_pion.push_back(pionTrack->getChi2());
			m_ndf_pion.push_back(pionTrack->getNdf());
			if (pionTrack->getNdf()!=0)
			{
				m_chi2_ndf_pion.push_back(pionTrack->getChi2()/pionTrack->getNdf());
			}

			m_chi2_kaon.push_back(kaonTrack->getChi2());
			m_ndf_kaon.push_back(kaonTrack->getNdf());
			if (kaonTrack->getNdf()!=0)
			{
				m_chi2_ndf_kaon.push_back(kaonTrack->getChi2()/kaonTrack->getNdf());
			}

			m_chi2_proton.push_back(protonTrack->getChi2());
			m_ndf_proton.push_back(protonTrack->getNdf());
			if (protonTrack->getNdf()!=0)
			{
				m_chi2_ndf_proton.push_back(protonTrack->getChi2()/protonTrack->getNdf());
			}
		}
	}
     catch (...)
     {
		streamlog_out(WARNING) << "Could not extract track collection for using pion mass " << std::endl;
     }
}
