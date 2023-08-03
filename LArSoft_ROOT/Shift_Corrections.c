// THIS CODE IS IN C LANGUAGE //

#include "import.h"
#include "tools.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "TGraph.h"
#include <string>
#include "TVector3.h"

using namespace art;

//Very beginning of the prototype for final "Optimization file infos generator"
// Fusion with main script if not expanded enough

//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//Shift Calculation Tools

double SharedFilledBins(
    std::vector<std::vector<double>> Pts, 
    std::vector<std::vector<double>> SimPts,
    string plotlabel,
    double BinSzX, double BinSzY, double BinSzZ,
    string hlabel,
    double xmin=-400,double xmax=400, 
    double ymin = -400 , double ymax = 400, 
    double zmin = 0, double zmax = 400)
{
    int N = Pts.at(0).size();
    int Nsim = SimPts.at(0).size();
    int NbinX = (xmax - xmin) / BinSzX;
    int NbinY = (ymax - ymin) / BinSzX;
    int NbinZ = (zmax - zmin) / BinSzZ;
    
    string hl = hlabel + "_" + plotlabel;
    TH3D *h = new TH3D(hl.c_str(), "Energy;X(cm);Y(cm);Z(cm)", NbinX, xmin, xmax, NbinY, ymin,ymax , NbinZ ,zmin , zmax);
    vector<int> QFilledBins;
    vector<int> EFilledBins;
    vector<int> SharedFilledBins;

    //Filling
    for (int i = 0; i < N; i++) 
    {
        double x = Pts.at(0).at(i);
        double y = Pts.at(1).at(i);
        double z = Pts.at(2).at(i);
        int bin = h->FindBin(x,y,z);
        bool New = true;
        for (int k = 0, sz = QFilledBins.size(); k < sz; k++)
        {
            int i = QFilledBins.at(k);
            if (i==bin)
            {
                New = false;
            }
        }
        if (New)
        {
            QFilledBins.push_back(bin);
        }
    }
    for (int i = 0; i < Nsim; i++) 
    {
        double x = SimPts.at(0).at(i);
        double y = SimPts.at(1).at(i);
        double z = SimPts.at(2).at(i);
        int bin = h->FindBin(x,y,z);
        bool New = true;
        for (int k = 0, sz = EFilledBins.size(); k < sz; k++)
        {
            int i = EFilledBins.at(k);
            if (i==bin)
            {
                New = false;
            }
        }
        if (New)
        {
            EFilledBins.push_back(bin);
        }
    }
    for (int k = 0, sz = QFilledBins.size(); k < sz; k++)
    {
        int i = QFilledBins.at(k);
        for (int l = 0, sze = EFilledBins.size(); l < sze; l++)
        {
            int j = EFilledBins.at(l);
            if (i==j)
            {
                SharedFilledBins.push_back(i);
            }
        }
    }
    return SharedFilledBins.size();
}

//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//Optimisation along 1 given axis

void plot_SFB(std::vector<std::vector<double>> Pts, 
    std::vector<std::vector<double>> SimPts,
    string plotlabel,
    double XShiftMin, double XShiftMax, int NXShift,
    double BinSzX, double BinSzY, double BinSzZ, 
    double xmin=-400,double xmax=400, 
    double ymin = -400 , double ymax = 400, 
    double zmin = 0, double zmax = 400)
{
    double XShift[NXShift], SFBX[NXShift];
    for (int i=0; i<NXShift ; i++)
    {
        double xshift = XShiftMin + (i*(XShiftMax-XShiftMin)/(NXShift-1));
        XShift[i] = xshift;
        SFBX[i] = SharedFilledBins(Shift_Pts(Pts,xshift,0,0,0) ,SimPts,plotlabel,BinSzX,BinSzY,BinSzZ,to_string(i) + "_x_",xmin,xmax,ymin,ymax,zmin,zmax);
    }
//     double YShift[NXShift], SFBY[NXShift];
//     for (int i=0; i<NXShift ; i++)
//     {
//         double yshift = XShiftMin + (i*(XShiftMax-XShiftMin)/(NXShift-1));
//         YShift[i] = yshift;
//         SFBY[i] = SharedFilledBins(Shift_Pts(Pts,0,yshift,0,0) ,SimPts,plotlabel,BinSzX,BinSzY,BinSzZ,to_string(i) + "_y_",xmin,xmax,ymin,ymax,zmin,zmax);
//     }
//     double ZShift[NXShift], SFBZ[NXShift];
//     for (int i=0; i<NXShift ; i++)
//     {
//         double zshift = XShiftMin + (i*(XShiftMax-XShiftMin)/(NXShift-1));
//         ZShift[i] = zshift;
//         SFBZ[i] = SharedFilledBins(Shift_Pts(Pts,0,0,zshift,0) ,SimPts,plotlabel,BinSzX,BinSzY,BinSzZ,to_string(i) + "_z_",xmin,xmax,ymin,ymax,zmin,zmax);
//     }
    string clabel = "cSFB_" + plotlabel;
    TCanvas* cSFB = new TCanvas(clabel.c_str(),"Sim/Reco xshift");
//     cSFB->Divide(3);
    
    cSFB->cd(1);
    TGraph* gX = new TGraph(NXShift,XShift,SFBX);
    gX->SetTitle("Nb of Bins Filled in both Simulation and Reconstruction vs shift on x axis;X Shift (cm);Nb of SharedFilledBins");
    gX->Draw("AC*");
    
//     cSFB->cd(2);    
//     TGraph* gY = new TGraph(NXShift,YShift,SFBY);
//     gY->SetTitle("Nb of Bins Filled in both Simulation and Reconstruction vs shift on y axis;Y Shift (cm);Nb of SharedFilledBins");
//     gY->Draw("AC*");
//     
//     cSFB->cd(3);    
//     TGraph* gZ = new TGraph(NXShift,ZShift,SFBZ);
//     gZ->SetTitle("Nb of Bins Filled in both Simulation and Reconstruction vs shift on z axis;Z Shift (cm);Nb of SharedFilledBins");
//     gZ->Draw("AC*");
}

double Opti_xshift(std::vector<std::vector<double>> Pts, 
    std::vector<std::vector<double>> SimPts,
    string plotlabel,
    double BinSzX, double BinSzY, double BinSzZ, 
    double eps = 0.1,
    double xmin=-400,double xmax=400, 
    double ymin = -400 , double ymax = 400, 
    double zmin = 0, double zmax = 400)

{
    double x(0);
    double sfb = SharedFilledBins(Pts,SimPts,plotlabel,BinSzX,BinSzY,BinSzZ,"",xmin,xmax,ymin,ymax,zmin,zmax);
    bool MaxNotFound = true;
    int i(0);
    double old_sfb(0);
    while (old_sfb<sfb)
    {
        x-=1;
        double new_sfb = SharedFilledBins(Shift_Pts(Pts,x,0,0,0), SimPts,plotlabel,BinSzX,BinSzY,BinSzZ,to_string(i),xmin,xmax,ymin,ymax,zmin,zmax);
        old_sfb=sfb;
        sfb = new_sfb;
        //cout << " x=" << x << " ; old : " << old_sfb << "  ;  new : " << sfb << endl; 
        i++;
    }
    old_sfb=sfb-1;
    while (old_sfb<sfb)
    {
        x+=eps;
        double new_sfb = SharedFilledBins(Shift_Pts(Pts,x,0,0,0), SimPts,plotlabel,BinSzX,BinSzY,BinSzZ,to_string(i)+"_",xmin,xmax,ymin,ymax,zmin,zmax);
        old_sfb=sfb;
        sfb = new_sfb;
        //cout << " x=" << x << " ; old : " << old_sfb << "  ;  new : " << sfb << endl; 
        i++;
    }
    return x;
}

//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//Prototype : Optimisation with Gradient Technique

vector<double> Grad_SFB(vector<double> X,double eps,
    std::vector<std::vector<double>> Pts, 
    std::vector<std::vector<double>> SimPts,
    string plotlabel,
    double BinSzX, double BinSzY, double BinSzZ,
    string hlg,
    double xmin=-400,double xmax=400, 
    double ymin = -400 , double ymax = 400, 
    double zmin = 0, double zmax = 400)
{
    cout << endl<< "*";
    int d = X.size();
    cout << "*";
    vector<double> G(d);
    cout << "*";
    for (int i=0; i<d; i++)
    {
        cout << endl<< "*";
        vector<double> Xm;
        vector<double> Xp;
        for (int j=0; j<d; j++)
        {
            if (j==i)
            {
                Xm[j] = X.at(j) - eps;
                Xp[j] = X.at(j) + eps;
            }
            else
            {
                Xm[j] = X.at(j);
                Xp[j] = X.at(j);
            }
        }
        cout << "*";
        string hl = hlg + to_string(i);
        double sfbp = SharedFilledBins(Shift_Pts(Pts,Xp.at(0),Xp.at(1),Xp.at(2)), SimPts,plotlabel,BinSzX,BinSzY,BinSzZ,hl+"p",xmin,xmax,ymin,ymax,zmin,zmax);
        cout<< "*";
        double sfbm = SharedFilledBins(Shift_Pts(Pts,Xm.at(0),Xm.at(1),Xm.at(2)), SimPts,plotlabel,BinSzX,BinSzY,BinSzZ,hl+"m",xmin,xmax,ymin,ymax,zmin,zmax);
        cout<< "*";
        G[i]=(sfbp-sfbm)/eps;
    }
    return G;
}

vector<double> MinSFB(std::vector<std::vector<double>> Pts, 
                        std::vector<std::vector<double>> SimPts,
                        string plotlabel,
                        double BinSzX, double BinSzY, double BinSzZ, 
                        int Niter =20,double tol=0.1,double eps = 0.5,
                        double xmin=-400,double xmax=400, 
                        double ymin = -400 , double ymax = 400, 
                        double zmin = 0, double zmax = 400)
{
    vector<double> X(3,0);
    double step(1);
    int n(0);
    bool MaxNotFound = true;
    while (MaxNotFound && n<Niter)
    {
        cout << "Step "<< (3*n) <<" OK" << endl;
        vector<double> g = Grad_SFB(X,eps,Pts,SimPts,plotlabel,BinSzX,BinSzY,BinSzZ,to_string(n));
        cout << "Step "<< (3*n) + 1 <<" OK" << endl;

        double normg = sqrt((g.at(0)*g.at(0)) + (g.at(1)*g.at(1)) + (g.at(2)*g.at(2)));
        X[0] -= step * g.at(0)/normg;
        X[1] -= step * g.at(1)/normg;
        X[2] -= step * g.at(2)/normg;
        if (normg<tol)
        {
            MaxNotFound = false;
        }
        step *= 0.8;
        n++;
        cout << "Step "<< (3*n) + 2 <<" OK" << endl;
    }
    return X;
}

//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//Terminal Requests

void Shift_tests()
{
    string label = "10elec_1GeV";
    int EVmin = 1;
    int EVmax = 1;
    double X[EVmax-EVmin+1];
    for (int evn = EVmin-1; evn<EVmax; evn++)
    {
        string evlabel = label + "_" + to_string(evn+1);
        cout << endl << endl << "-----------------  EVENT #" << evn+1 << "  ------------------------"<< endl;
        std::vector<std::vector<double>> Pts = XYZQ(label, evn);
        std::vector<std::vector<double>> SimPts = XYZE(label, evn);
//         plot_SFB(Pts,SimPts,evlabel,-0.5,0.5,25,1,1,1);
//         X[evn] = Opti_xshift(Pts,SimPts,evlabel,0.5,0.5,0.5);
//         plot_SFB(Pts,SimPts,evlabel+"_precize",-5,-1,50,1,1,1);
        vector<double> X = MinSFB(Pts,SimPts,evlabel,0.5,0.5,0.5);
        cout << endl << "X Shift : " << X.at(0) << "  ;  Y Shift : " << X.at(1) << "  ;  Z Shift : " << X.at(2) << endl;
        cout <<endl << endl << "***************************************************************" << endl;
    }
    for (int evn = EVmin-1; evn<EVmax; evn++)
    {
            cout << endl << "Event #" << evn << "  ->  Shift along X : " << X[evn]  << " cm" << endl;
    }
}

void Edit_Shift_File(string label)
{
    string adress = "/silver/DUNE/peron/simu/sim_10elec_1GeV_less_treshold/" + label + "_Shift_Corrections.txt";
    ofstream Txt(adress.c_str());
    
    if (Txt)
    {
        int EVmin = 1;
        int EVmax = 10;
        double X[EVmax-EVmin+1];
        for (int evn = EVmin-1; evn<EVmax; evn++)
        {
            string evlabel = label + "_" + to_string(evn+1);
            cout << endl << endl << "-----------------  EVENT #" << evn+1 << "  ------------------------"<< endl;
            std::vector<std::vector<double>> Pts = XYZQ(label, evn);
            std::vector<std::vector<double>> SimPts = XYZE(label, evn);
            X[evn] = Opti_xshift(Pts,SimPts,evlabel,1,1,1);
        }
        
        for (int evn = EVmin-1; evn<EVmax; evn++)
        {
            cout << "Event #" << evn+1 << "  ->  Shift along X : " << X[evn]  << " cm" << endl;
            Txt << X[evn] << endl;
        }
    }
    else
    {
        cout << endl << "ERROR : File can't be opened"<<endl;
    }
    
    // pour lire : ifstream Txt(adress.c_str()); 
}


void Shift_Corrections()
{
    string label = "10muon_10GeV";
    int EVmin = 1;
    int EVmax = 10;
    for (int evn = EVmin-1; evn<EVmax; evn++)
    {
        string evlabel = label + "_" + to_string(evn+1);
        cout << endl << endl << "-----------------  EVENT #" << evn+1 << "  ------------------------"<< endl;
        std::vector<std::vector<double>> Pts = XYZQ(label, evn);
        std::vector<std::vector<double>> SimPts = XYZE(label, evn);
        std::vector<double> O = Origin(label,evn);
        plot_SFB(Pts,SimPts,evlabel,-5,1,120,2,2,2);
    }
}
