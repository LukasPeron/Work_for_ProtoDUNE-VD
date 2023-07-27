// THIS CODE IS IN C LANGUAGE //

#include "import.h"
#include "tools.h"
#include <vector>
#include <cmath>
#include <iostream>
#include "TGraph.h"
#include <string>

using namespace art;

//Prototype for final "Stat file generator"

//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//Data ordering :

void ex_stat_builder(string label)
{
    int EVmin = 1;
    int EVmax = 10;
    for (int evn = EVmin-1; evn<EVmax; evn++)
    {
        string evlabel = label + "_" + to_string(evn+1);
        cout << endl << endl << "-----------------  EVENT #" << evn+1 << "  ------------------------"<< endl;
        std::vector<std::vector<double>> Pts = XYZQ(label, evn);
        std::vector<std::vector<double>> SimPts = XYZE(label, evn);
        std::vector<double> O = Origin(label,evn);
        
    }
}

vector<double> UnitVectShowerAxe(vector<vector<double>> Pts,vector<double> O)
{
    std::vector<double> B = Barycenter(Pts,Pts.at(3));
    vector<double> v = Vect2Pts(O,B);
    double vNorm = sqrt(scalar(v,v));
    v[0]/=vNorm;
    v[1]/=vNorm;
    v[2]/=vNorm;    
    return v;
}

vector<double> DistancesBarycentres(string label)
{
    int EVmin = 1;
    int EVmax = 10;
    vector<double> DBs;
    for (int evn = EVmin-1; evn<EVmax; evn++)
    {
        string evlabel = label + "_" + to_string(evn+1);
        cout << endl << endl << "-----------------  EVENT #" << evn+1 << "  ------------------------"<< endl;
        std::vector<std::vector<double>> Pts = XYZQ(label, evn);
        std::vector<std::vector<double>> SimPts = XYZE(label, evn);
        
        //Barycenters
        std::vector<double> Charges = Pts.at(3);
        std::vector<double> Energies = SimPts.at(3);
        std::vector<double> O = Origin(label,evn);
        Pts = ShiftXOpti_Pts(Pts,label,evn);
        std::vector<double> B = Barycenter(Pts, Charges);
        std::vector<double> Bs = Barycenter_simu(SimPts, Energies);
        std::vector<double> axis = UnitVectShowerAxe(SimPts, O);
        
        DBs.push_back(Distance_algebrique(B,Bs,axis));
    }
    return DBs;
}

void EnergyRecoStat(string label,
                    double BinSzX, double BinSzY, double BinSzZ, 
    double xmin=0,double xmax=400, 
    double ymin = -400 , double ymax = 400, 
    double zmin = -400, double zmax = 400)
{
    //Infos à stater :
    //Nb de bins communs, charge tot et énergie tot, rapport Etot/Qtot
    //SharedFilledBins
    //Moyenne 
    //etc.
    
    int EVmin = 1;
    int EVmax = 10;
    // vector<double> DBs; //not useful here
    for (int evn = EVmin-1; evn<EVmax; evn++)
    {
        string evlabel = label + "_" + to_string(evn+1);
        cout << endl << endl << "-----------------  EVENT #" << evn+1 << "  ------------------------"<< endl;
        std::vector<std::vector<double>> Pts = XYZQ(label, evn);
        std::vector<std::vector<double>> SimPts = XYZE(label, evn);
    
        //Useful const
        int N = Pts.at(0).size();
        int Nsim = SimPts.at(0).size();
        int NbinX = (xmax - xmin) / BinSzX;
        int NbinY = (ymax - ymin) / BinSzY;
        int NbinZ = (zmax - zmin) / BinSzZ;
        
        cout << endl << "***" << endl;
        cout << "Nb of Simulated Depositions : " << Nsim << endl;
        cout << "Nb of Reconstructed Depositions : " << N << endl;
        
        cout << endl << "***" << endl << "Bins Info" << endl;
        cout << "X Direction : " << NbinX << "bins of size " << BinSzX << " cm (" << xmin << "cm < x < " << xmax << "cm)" << endl;
        cout << "Y Direction : " << NbinY << "bins of size " << BinSzY << " cm (" << ymin << "cm < y < " << ymax << "cm)" << endl;
        cout << "Z Direction : " << NbinZ << "bins of size " << BinSzZ << " cm (" << zmin << "cm < z < " << zmax << "cm)" << endl;
        
        cout << endl << "***" << endl << "Normalization process : " << endl;
        //Normalization
        double Qtot(0);
        double Etot(0);
        for (int i = 0; i < N; i++) 
        {
            Qtot+=Pts.at(3).at(i);
        }
        
            for (int i = 0; i < Nsim; i++) 
        {
            Etot+=SimPts.at(3).at(i);
        }
        
        cout << "Integrated Charge (Reco) : " << Qtot << " ADC" << endl;
        cout << "Integrated Energy (Sim) : " << Etot << " MeV" << endl;
    }
}

//----------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//Stat Histogram
void plot_Hist(vector<double> List, string hlabel)
{
    TCanvas *c = new TCanvas(hlabel.c_str(),"Histogram");
    c->cd();
    
    TH1D *h = new TH1D(hlabel.c_str(),"Ecart Barycentres;Distance between SIM and RECO Barycenters (cm);NbEvent",40,-5,50);
    for (int i=0, sz = List.size();i<sz;i++)
    {
        h->Fill(List.at(i));
    }
    h->Draw("HIST");
    
//     auto legend = new TLegend();
//     legend->SetHeader("Legend","C");
//     legend->AddEntry(h,"")
//     legend->Draw()
}

//----------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//Stat Infos


void Stats()
{
    string label = "10muon_10GeV";
    plot_Hist(DistancesBarycentres(label),label);
}
