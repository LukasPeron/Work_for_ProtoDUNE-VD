// THIS CODE IS A HEADER IN C LANGUAGE //

#ifndef IMPORT_H_INCLUDED
#define IMPORT_H_INCLUDED

// LArSoft Libraries
R__ADD_INCLUDE_PATH("lardataobj/RawData/RawDigit.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Wire.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Hit.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Shower.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Track.h")
R__ADD_INCLUDE_PATH("gallery/Event.h")

// Standard libraries
#include <vector>
#include <iostream>
#include <string>

using namespace art;

//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//IMPORT RECO

//Import a vector of 4 lists (X,Y,Z,Q) carrying respectively the x,y,z coordinates and the charge q deposited on the collection plane for all reconstructed energy depositions

std::vector<std::vector<double>> GetTrack(string path, string label, string ext, int EvN=0, string ang=""){
    string rootfilelabel = path + label + ext;
    string trcktag = "pandoraTrack";
    art::InputTag trck_tag(trcktag);
    vector<string> filename(1,rootfilelabel);
    gallery::Event ev(filename);
    for (int i=0; i<EvN; i++)
    {
        ev.next();
    }
    auto Tracks = ev.getValidHandle<vector<recob::Track>>(trck_tag);
    float ftrackLen;

    ofstream Flux((path+label+"_track_size_"+ang+".txt").c_str(), ios::app);
    ofstream Flux2((path+label+"_track_lenght_"+ang+".txt").c_str(), ios::app);

    std::vector<std::vector<double>> Track_info(4);
    Flux << EvN << ";" << Tracks->size() << endl;
    for (unsigned itrk = 0; itrk < Tracks->size(); ++itrk) {
        const recob::Track& track = Tracks->at(itrk);
        ftrackLen = track.Length();
        Track_info.at(3).push_back(track.Chi2());
        Flux2 << EvN << ";" << ftrackLen << ";" << track.Chi2() << endl;
        for(int i(0);i<track.NumberTrajectoryPoints();i++)
        {
            if(track.LocationAtPoint(i).X()>-999)
            {
                Track_info.at(0).push_back(track.LocationAtPoint(i).X());
                Track_info.at(1).push_back(track.LocationAtPoint(i).Y());
                Track_info.at(2).push_back(track.LocationAtPoint(i).Z());
            }
        }
    }

    return Track_info;
}

std::vector<std::vector<double>> XYZQ(
    string path,
    string label, // "NbOfEvents TypeOfParticle _ Energy" -> ex : "10elec_1GeV"
    string ext,
    int EvN = 0 // Event number
    ){ 
    
    //Path to the Reco root file associated to the label : 
    
    string rootfilelabel = path + label + ext;
    //module for SpacePoints : reco3d
    string recotag = "reco3d"; //pandora
    art::InputTag reco_tag(recotag);
    vector<string> filename(1,rootfilelabel);
    


    //Event selection
    gallery::Event ev(filename);
    for (int i=0; i<EvN; i++)
    {
        ev.next();
    }
    // auto const depolist = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(rootfilelabel);
    // cout << endl << "Nombre de depots : " << depolist->size() << endl;
    auto SPs = ev.getValidHandle<vector<recob::SpacePoint>>(reco_tag); //Get the list of reconstructed SpacePoints (SPs)
    art::FindManyP<recob::Hit> fmHits(SPs, ev, reco_tag); //Get the list of Hits from each plane associated to each SP
    
    //Vairiables Initialisation
    std::vector<std::vector<double>> Pts(4);
    double Q(0);
    
    for (size_t j = 0, sz = SPs->size(); j != sz; ++j)
    {
        const std::vector<art::Ptr<recob::Hit>>& Hits = fmHits.at(j);
        for (size_t i = 0, sz = Hits.size(); i != sz; ++i)
        {
                if (Hits[i]->WireID().Plane == 2) 
                //Only Hits integral from the collection plane (ID = 2) carry real charge information
                {
                    
                	//Condition to avoid twice counted SPs
                    //double newQ = Hits.at(i)->Integral();
                    double newQ = Hits.at(i)->SummedADC();
                    if (newQ!=Q)
                	{
                        //Charge
                        if(auto it = std::find(begin(Pts.at(3)), end(Pts.at(3)), newQ); it != std::end(Pts.at(3))){ cout << i << " " << j << " " << newQ  << " " << *it << endl;}
                        Pts.at(3).push_back(newQ);
				
                        //SpacePoints 
                        double const X = SPs->at(j).XYZ()[0];
                        Pts.at(0).push_back(X);
                        double const Y = SPs->at(j).XYZ()[1];
                        Pts.at(1).push_back(Y);
                        double const Z = SPs->at(j).XYZ()[2];
                        Pts.at(2).push_back(Z);

                		Q = newQ;
                        // cout << i << " " << j << " " << Q << endl;
                		// cout << Q << " (X="<<X << " Y=" << Y << " Z="<< Z<< ")"<< endl;
                	}
                }
        }
        
    }
    
    double qtot(0);
    for (int i(0); i<Pts.at(3).size(); i++)
            {
                qtot+=Pts.at(3).at(i);
            }
            cout << "event_sim_reco qtot = " << qtot << endl;

    return Pts;
}

std::vector<std::vector<double>> XYZQ_pandora(
    string path, 
    string label, // "NbOfEvents TypeOfParticle _ Energy" -> ex : "10elec_1GeV"
    string ext,
    int EvN = 0 // Event number
    ){ 
    
    //Path to the Reco root file associated to the label : 
    string rootfilelabel = path + label + ext;

    //module for SpacePoints : reco3d
    string recotag = "pandora"; //pandora
    art::InputTag reco_tag(recotag);
    vector<string> filename(1,rootfilelabel);
    
    //Event selection
    gallery::Event ev(filename);
    for (int i=0; i<EvN; i++)
    {
        ev.next();
    }
    
    // auto const depolist = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(rootfilelabel);
    // cout << endl << "Nombre de depots : " << depolist->size() << endl;
    auto SPs = ev.getValidHandle<vector<recob::SpacePoint>>(reco_tag); //Get the list of reconstructed SpacePoints (SPs)
    art::FindManyP<recob::Hit> fmHits(SPs, ev, reco_tag); //Get the list of Hits from each plane associated to each SP
    
    //Vairiables Initialisation
    std::vector<std::vector<double>> Pts(4);
    double Q(0);
    
    for (size_t j = 0, sz = SPs->size(); j != sz; ++j)
    {
        const std::vector<art::Ptr<recob::Hit>>& Hits = fmHits.at(j);
        for (size_t i = 0, sz = Hits.size(); i != sz; ++i)
        {
                if (Hits[i]->WireID().Plane == 2) 
                //Only Hits integral from the collection plane (ID = 2) carry real charge information
                {
                	//Condition to avoid twice counted SPs
                    //double newQ = Hits.at(i)->Integral();
                    double newQ = Hits.at(i)->SummedADC();
                    if (newQ!=Q)
                	{
                        //Charge
                        if(auto it = std::find(begin(Pts.at(3)), end(Pts.at(3)), newQ); it != std::end(Pts.at(3))){ cout << i << " " << j << " " << newQ  << " " << *it << endl;}
                        Pts.at(3).push_back(newQ);
				
                        //SpacePoints 
                        double const X = SPs->at(j).XYZ()[0];
                        Pts.at(0).push_back(X);
                        double const Y = SPs->at(j).XYZ()[1];
                        Pts.at(1).push_back(Y);
                        double const Z = SPs->at(j).XYZ()[2];
                        Pts.at(2).push_back(Z);
				
                		Q = newQ;
                        // cout << i << " " << j << " " << Q << endl;
                		// cout << Q << " (X="<<X << " Y=" << Y << " Z="<< Z<< ")"<< endl;
                	}
                }
        }
        
    }

    double qtot(0);
    for (int i(0); i<Pts.at(3).size(); i++)
            {
                qtot+=Pts.at(3).at(i);
            }
            cout << "event_sim_reco qtot = " << qtot << endl;

    return Pts;
}
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//IMPORT SIM G4

//Import a vector of 4 lists (X,Y,Z,E) carrying respectively the x,y,z coordinates and the energy deposited in the LArTPC for all energy depositions simulated by G4
std::vector<std::vector<double>> XYZE(
    string path,
    string label, // "NbOfEvents TypeOfParticle _ Energy" -> ex : "10elec_1GeV"
    string ext,
    int EvN = 0 // Event number
    ){
    
    //Path to the Sim root file associated to the label :
    string rootfilelabel = path + label + ext;
    
    //module for the G4 sim infos :
    string simtag = "largeant:LArG4DetectorServicevolTPCActive"; // "largeant:LArG4DetectorServicevolTPCActive"; IonAndScint
    art::InputTag sim_tag(simtag);
    vector<string> filename(1, rootfilelabel);
    
    //Event selection
    gallery::Event ev(filename); 
    for (int i=0; i<EvN; i++)
    {
        ev.next();
    }    
    
    auto const depolist = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(sim_tag);
    
    //Initialisation
    size_t ndepos = depolist->size();
    vector<vector<double>> Pts(4);
    
    for (size_t depo_i = 0; depo_i < ndepos; depo_i++)
    {
        const sim::SimEnergyDeposit& depo = depolist->at(depo_i);
        geo::Point_t xyz = depo.MidPoint(); // MidPoint contain the coordinates of the center of the deposition
        
        std::vector<double> Pt; // Vector to fill with x,y,z coordinates
        Pts.at(0).push_back(xyz.X());
        Pts.at(1).push_back(xyz.Y());
        Pts.at(2).push_back(xyz.Z());
        
        // double e = depo.Energy(); //Energy of the Deposition
        double ne = depo.NumElectrons();
        Pts.at(3).push_back(ne);
    }
    
    return Pts;
}


std::vector<std::vector<double>> GetOriginSimDepo(string path, 
                                     string label, 
                                     string ext, 
                                     int EvN=0)
{
    string rootfilelabel = path + label + ext;
    
    //module for the G4 sim infos :
    string simtag = "IonAndScint"; // "largeant:LArG4DetectorServicevolTPCActive";
    art::InputTag sim_tag(simtag);
    vector<string> filename(1, rootfilelabel);
    
    //Event selection
    gallery::Event ev(filename); 
    for (int i=0; i<EvN; i++)
    {
        ev.next();
    }    
    
    auto const depolist = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(sim_tag);
    
    //Initialisation
    size_t ndepos = depolist->size();
    std::vector<std::vector<double>> Pts(3);
    
    for (size_t depo_i = 0; depo_i < ndepos; depo_i++)
    {
        const sim::SimEnergyDeposit& depo = depolist->at(depo_i);
        if(depo.OrigTrackID()==0 && depo.PdgCode() == 11) //bonne origine + c'est un electron
        {
            geo::Point_t xyz = depo.MidPoint(); // MidPoint contain the coordinates of the center of the deposition
            // cout << "Find one original electron !"<< endl;
            Pts.at(0).push_back(xyz.X());
            Pts.at(1).push_back(xyz.Y());
            Pts.at(2).push_back(xyz.Z());
        }
        
    }
    
    return Pts;
}


//Import the (x,y,z) coordinates of the origin of the beam choosen for the simulation
std::vector<double> Origin(string label, int EvN = 0)
{
    std::vector<double> O(3);
    /*REFERENCE 
        O[0] = 94.8;
        O[1] = -142.6;
        O[2] = 299.3;
    */
    //Currently an operationnal function but bound to become an automatic way to access simulation parameters through TBranches
    O[0] = 10;
    O[1] = 10;
    O[2] = 100;
    return O;
}

#endif
