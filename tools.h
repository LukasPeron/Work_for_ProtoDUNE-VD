// THIS CODE IS A HEADER IN C LANGUAGE //

#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED

#include "import.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "TGraph.h"
#include <string>
#include "TVector3.h"

using namespace art;

//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//Geometry tools

//Return the barycenter coordinates given a list of Points coordinates et associated Weights
std::vector<double> Barycenter(std::vector<std::vector<double>> Pts, std::vector<double> Weight)
{
    std::vector<double> B(3,0);
    int sz = Pts.at(0).size();
    double Wtot(0);
    for (int i=0; i!=sz;++i)
    {
        double x = Pts.at(0).at(i);
        double y = Pts.at(1).at(i);
        double z = Pts.at(2).at(i);
        double w = Weight.at(i);
        Wtot+=w;
        B[0]+= w*x;
        B[1]+= w*y;
        B[2]+= w*z;
    }
    B[0] /= Wtot;
    B[1] /= Wtot;
    B[2] /= Wtot;
    return B;
}

std::vector<double> Barycenter_simu(std::vector<std::vector<double>> Pts, std::vector<double> Weight)
{
    std::vector<double> B(3,0);
    int sz = Pts.at(0).size();
    double Wtot(0);
    double drift_velocity = 0.15*1e6;
    double lifetime = 35*1e-3;
    for (int i=0; i!=sz;++i)
    {
        double x = Pts.at(0).at(i);
        double y = Pts.at(1).at(i);
        double z = Pts.at(2).at(i);
        double t = x/drift_velocity;
        double w = Weight.at(i)*exp(-t/lifetime);
        Wtot+=w;
        B[0]+= w*x;
        B[1]+= w*y;
        B[2]+= w*z;
    }
    B[0] /= Wtot;
    B[1] /= Wtot;
    B[2] /= Wtot;
    return B;
}

//---------------------------------------------------------------------------------------------
// return M1M2 vect with given points M1 and M2 coordinates
vector<double> Vect2Pts(vector<double> M1, vector<double> M2) 
{
    int sz = M1.size();
    vector<double> v(sz);
    for (int i=0; i<sz; i++)
    {
        v.at(i) = M2.at(i) - M1.at(i);
    }
    return v;
}

//---------------------------------------------------------------------------------------------
// return scalar product of given vectors v1 and v2
double scalar(vector<double> v1 , vector<double> v2) 
{
    int sz = v1.size();
    double s(0);
    for (int i=0; i<sz; i++)
    {
        s += (v1.at(i))*(v2.at(i));
    }
    return s;
}

//---------------------------------------------------------------------------------------------
// return distance between points M1 and M2 given their coordinates
double Distance(vector<double> M1, vector<double> M2) 
{
    double d(0);
    for (int i=0, sz=M1.size(); i<sz; i++)
    {
        d+=pow((M2.at(i) - M1.at(i)),1);
    }
    return d;
}

double Distance_algebrique(vector<double> M1, vector<double> M2, vector<double> axis) 
{
    double d = (scalar(Vect2Pts(M1, M2), axis))/sqrt(scalar(axis, axis));
    return d;
}


//---------------------------------------------------------------------------------------------
// return the derivative of a vector Y along a vecteur X
vector<double> Derivative(const vector<double>& X, const vector<double>& Y)
{
    std::vector<double> derivee(Y.size());
    for (int i = 1; i < Y.size() - 1; i++) {
        double h1 = X[i] - X[i - 1];
        double h2 = X[i + 1] - X[i];

        derivee[i] = (Y[i + 1] - Y[i - 1]) / (h1 + h2);
    }

    // Calcul de la dérivée aux extrémités
    double h = X[1] - X[0];
    derivee[0] = (Y[1] - Y[0]) / h;

    h = X[X.size() - 1] - X[X.size() - 2];
    derivee[Y.size() - 1] = (Y[Y.size() - 1] - Y[Y.size() - 2]) / h;

    return derivee;
}


//---------------------------------------------------------------------------------------------
//return new list of coordinates changing the origin :
vector<vector<double>> XYZ_to_OM(vector<vector<double>> Pts , vector<double> O /*New Origin*/)
{
    int sz = Pts.at(0).size();
    vector<vector<double>> OMPts(sz);
    for (int i=0; i<sz; i++)
    {
        vector<double> OM(3);
        OM.at(0) = Pts.at(0).at(i) - O.at(0);
        OM.at(1) = Pts.at(1).at(i) - O.at(1);
        OM.at(2) = Pts.at(2).at(i) - O.at(2);
        OMPts[i] = OM;
    }
    return OMPts;
}

//---------------------------------------------------------------------------------------------
//Change default cartesian coordinates list to axis coordinates list :
vector<vector<double>> XYZ_to_rl(vector<vector<double>> Pts,
                                 vector<double> O, /*Shower Origin*/ 
                                 vector<double> v /*Shower axis*/){
    int sz = Pts.at(0).size();
    vector<vector<double>> OMPts = XYZ_to_OM(Pts,O);
    vector<vector<double>> AxPts(sz);
    for (int i=0; i<sz; i++)
    {
        vector<double> OM = OMPts.at(i);
        double l = scalar(OM,v);
        double OM2 = scalar(OM,OM);
        double r = sqrt(OM2 - (l*l));
        AxPts[0].push_back(l);
        AxPts[1].push_back(r);
    }
    return AxPts;
}

//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//Shift tools

//return the new (X,Y,Z,W=E/Q) after a given shift
std::vector<std::vector<double>> Shift_Pts(std::vector<std::vector<double>> Pts,
                                           double xshift, // default : only a shift along X axis
                                           double yshift=0,
                                           double zshift=0,
                                           double wshift=0) //Energy or Charge Shift
{
    int sz=Pts.at(0).size();
    std::vector<std::vector<double>> Shifted_Pts(4);
    for (int i=0;i<sz;i++)
    {
        Shifted_Pts.at(0).push_back(Pts.at(0).at(i) + xshift);
        Shifted_Pts.at(1).push_back(Pts.at(1).at(i) + yshift);
        Shifted_Pts.at(2).push_back(Pts.at(2).at(i) + zshift);
        Shifted_Pts.at(3).push_back(Pts.at(3).at(i) + wshift);
    }
    return Shifted_Pts;
}

//---------------------------------------------------------------------------------------------

//return the new (X,Y,Z,W=E/Q) after an optimized X-Shift (cf Shift Corrections)
std::vector<std::vector<double>> ShiftXOpti_Pts(std::vector<std::vector<double>> Pts, string label, int EvN=0)
{
    //Read the file with stored calculated optimized X-Shifts
    // string adress = "/silver/DUNE/queau/output/root_" + label + "/" + label + "_Shift_Corrections.txt";
    // ifstream Txt(adress.c_str());
    // int evn(0);
    // double xshift;
    // Txt >> xshift;
    // while (evn!=EvN)
    // {
    //     Txt >> xshift;
    //     evn++;
    // }
    // //Shift with previous code
    return Shift_Pts(Pts,-2.2); //-2.8 for EM shower
}

//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//Plot tools

//Optimized Window size
std::vector<double> WindowSz(vector<vector<double>> Pts /*If RecoPts and SimPts in the window, use SimPts*/)
{
    //Get the min and max coordinates to see everything on screen
    double xmin = Pts.at(0).at(std::min_element(Pts.at(0).begin(),Pts.at(0).end())-Pts.at(0).begin());
    double xmax = Pts.at(0).at(std::max_element(Pts.at(0).begin(),Pts.at(0).end())-Pts.at(0).begin());
    double ymin = Pts.at(1).at(std::min_element(Pts.at(1).begin(),Pts.at(1).end())-Pts.at(1).begin());
    double ymax = Pts.at(1).at(std::max_element(Pts.at(1).begin(),Pts.at(1).end())-Pts.at(1).begin());
    double zmin = Pts.at(2).at(std::min_element(Pts.at(2).begin(),Pts.at(2).end())-Pts.at(2).begin());
    double zmax = Pts.at(2).at(std::max_element(Pts.at(2).begin(),Pts.at(2).end())-Pts.at(2).begin());
       
    std::vector<double> Window(6,0);
    Window[0]= xmin;
    Window[1]= xmax;
    Window[2]= ymin;
    Window[3]= ymax;
    Window[4]= zmin;
    Window[5]= zmax;
    
    return Window;
}

//---------------------------------------------------------------------------------------------
//3D Scatter for given list of Points coordinates
TPolyMarker3D* Make_Poly_3D(std::vector<std::vector<double>> Pts)
{
    TPolyMarker3D* Poly = new TPolyMarker3D(Pts.size());

    for (unsigned j=0 , sz = Pts.at(0).size(); j!=sz;++j)
    {
        double x = Pts.at(0).at(j);
        double y = Pts.at(1).at(j);
        double z = Pts.at(2).at(j);
        Poly->SetPoint(j,x,y,z);
    }
    
    Poly->SetMarkerStyle(2);
    Poly->SetMarkerSize(0.5);
    return Poly;
}

//---------------------------------------------------------------------------------------------
//2D Projections

//Add a scatter plot to 2D projections multigraphs
void plot_proj_pts(vector<vector<double>> Pts, //List of Points Coordinates
               TMultiGraph* mgXY, TMultiGraph* mgYZ,TMultiGraph* mgXZ, //Projections multigraphs
               int color, int marker, //scatter plot settings
               TLegend* legend, string legend_label //legend and entry label
                  ){

    int sz = Pts.at(0).size();
    double X[sz], Y[sz], Z[sz];
    
    for (unsigned j=0; j!=sz;++j)
    {
        X[j] = Pts.at(0).at(j);
        Y[j] = Pts.at(1).at(j);
        Z[j] = Pts.at(2).at(j);
    }
    
    TGraph* gXY = new TGraph(sz,X,Y);
    gXY->SetMarkerColor(color);
    gXY->SetMarkerStyle(marker);
    mgXY->Add(gXY,"P");
    
    TGraph* gYZ = new TGraph(sz,Y,Z);
    gYZ->SetMarkerColor(color);
    gYZ->SetMarkerStyle(marker);
    mgYZ->Add(gYZ,"P");
    
    TGraph* gXZ = new TGraph(sz,X,Z);
    gXZ->SetMarkerColor(color);
    gXZ->SetMarkerStyle(marker);
    mgXZ->Add(gXZ,"P");
    
    TGraph* l = new TGraph();
    l->SetMarkerColor(color);
    
    if (marker == 6)
    {
    	l->SetMarkerStyle(34);
    }
    
    l->SetMarkerSize(3);
    
    legend->AddEntry(l,legend_label.c_str(),"P");
    
}

//---------------------------------------------------------------------------------------------
//Add a point to 2D projections multigraphs
void plot_proj_pt(vector<double> Pt, //Point Coordinates
               TMultiGraph* mgXY, TMultiGraph* mgYZ,TMultiGraph* mgXZ, //Projections multigraphs
               int color, int marker, //point plot settings
               TLegend* legend, string legend_label //legend and entry label
                  ){
    
    double X[1], Y[1], Z[1];

    X[0] = Pt.at(0);
    Y[0] = Pt.at(1);
    Z[0] = Pt.at(2);
    
    TGraph* gXY = new TGraph(1,X,Y);
    gXY->SetMarkerColor(color);
    gXY->SetMarkerStyle(marker);
    mgXY->Add(gXY,"P");
    
    TGraph* gYZ = new TGraph(1,Y,Z);
    gYZ->SetMarkerColor(color);
    gYZ->SetMarkerStyle(marker);
    mgYZ->Add(gYZ,"P");
    
    TGraph* gXZ = new TGraph(1,X,Z);
    gXZ->SetMarkerColor(color);
    gXZ->SetMarkerStyle(marker);
    mgXZ->Add(gXZ,"P");
    
    legend->AddEntry(gXY,legend_label.c_str(),"P");
}


//---------------------------------------------------------------------------------------------

//Add a track to 2D projections multigraphs
void plot_proj_line(vector<vector<double>> Pts, //List of Points Coordinates
               TMultiGraph* mgXY, TMultiGraph* mgYZ,TMultiGraph* mgXZ, //Projections multigraphs
               int color, //line plot settings
               TLegend* legend, string legend_label //legend and entry label
                  ){
    
    int sz = Pts.at(0).size();
    double X[sz], Y[sz], Z[sz];
    
    for (unsigned j=0; j!=sz;++j)
    {
        X[j] = Pts.at(0).at(j);
        Y[j] = Pts.at(1).at(j);
        Z[j] = Pts.at(2).at(j);
    }
    
    TGraph* gXY = new TGraph(sz,X,Y);
    gXY->SetLineColor(color);
    mgXY->Add(gXY,"l");
    
    TGraph* gYZ = new TGraph(sz,Y,Z);
    gYZ->SetLineColor(color);
    mgYZ->Add(gYZ,"l");
    
    TGraph* gXZ = new TGraph(sz,X,Z);
    gXZ->SetLineColor(color);
    mgXZ->Add(gXZ,"l");
    
    legend->AddEntry(gXY,legend_label.c_str(),"l");
}

//---------------------------------------------------------------------------------------------
//Add a line between 2 Points to 2D projections multigraphs
void plot_proj_line_2pts(vector<double> Pt1,vector<double> Pt2, //Points coordinates
               TMultiGraph* mgXY, TMultiGraph* mgYZ,TMultiGraph* mgXZ, //Projections multigraphs
               int color, //line plot settings
               TLegend* legend, string legend_label //legend and entry label
                  ){
    
    double x1(Pt1.at(0)), x2(Pt2.at(0)), y1(Pt1.at(1)),y2(Pt2.at(1)), z1(Pt1.at(2)), z2(Pt2.at(2));
    
    double X[2], Y[2], Z[2];

    X[0] = x1;
    Y[0] = y1;
    Z[0] = z1;
    X[1] = x2;
    Y[1] = y2;
    Z[1] = z2;
    
    TGraph* gXY = new TGraph(2,X,Y);
    gXY->SetLineColor(color);
    mgXY->Add(gXY,"l");
    
    TGraph* gYZ = new TGraph(2,Y,Z);
    gYZ->SetLineColor(color);
    mgYZ->Add(gYZ,"l");
    
    TGraph* gXZ = new TGraph(2,X,Z);
    gXZ->SetLineColor(color);
    mgXZ->Add(gXZ,"l");
    
    legend->AddEntry(gXY,legend_label.c_str(),"l");
}


//---------------------------------------------------------------------------------------------
//2D Hist

//Add 2D projections colormaps to (3) Canvas boxes following box1_nb (inclued)
void plot_proj_Hist(vector<vector<double>> Pts, //List of Points Coordinates
                    TCanvas* c, //Canvas
                    int box1_nb, //Nb of the box for the 1st projection (Y(X))
                    string plotlabel, //Name to distinguish different plots
                    double BinSzX, double BinSzY,double BinSzZ, //Bins dimensions
                    
                    //Default (else) : Detector dimensions
                    //Note : à adapter au réel
                    double xmin=0,double xmax=400, 
		     double ymin = -400 , double ymax = 400, 
		     double zmin = -400, double zmax = 400)
{
    int N = Pts.at(0).size();
    int NbinX = (xmax - xmin) / BinSzX;
    int NbinY = (ymax - ymin) / BinSzY;
    int NbinZ = (zmax - zmin) / BinSzZ;
    
    //Integral for normalisation
    double Wtot(0);
    for (int i = 0; i < N; i++) 
    {
        Wtot+=Pts.at(3).at(i);
    }
    
    string titleXY = plotlabel + " Y(X);X(cm);Y(cm)";
    string hXYlabel = "hXY_" + plotlabel;
    TH2D *hXY = new TH2D(hXYlabel.c_str(),titleXY.c_str(), NbinX, xmin, xmax, NbinY, ymin,ymax);
    
    string titleYZ = plotlabel + " Z(Y);Y(cm);Z(cm)";
    string hYZlabel = "hYZ_" + plotlabel;
    TH2D *hYZ = new TH2D(hYZlabel.c_str(), titleYZ.c_str(), NbinY, ymin, ymax, NbinZ, zmin,zmax);
    
    string titleXZ = plotlabel + " Z(X);X(cm);Z(cm)";
    string hXZlabel = "hXZ_" + plotlabel;
    TH2D *hXZ = new TH2D(hXZlabel.c_str(), titleXZ.c_str(), NbinX, xmin, xmax, NbinZ, zmin,zmax);
    
    for (int i = 0; i < N; i++) 
    {
        double x = Pts.at(0).at(i);
        double y = Pts.at(1).at(i);
        double z = Pts.at(2).at(i);
        double W = Pts.at(3).at(i)/Wtot;
        hXY->Fill(x,y,W);
        hYZ->Fill(y,z,W);
        hXZ->Fill(x,z,W);
    }
    
    c->cd(box1_nb);
    c->cd(box1_nb)->SetLogz();
    hXY->SetStats(kFALSE);
    hXY->Draw("COLZ");
    
    c->cd(box1_nb+1);
    c->cd(box1_nb+1)->SetLogz();
    hXZ->SetStats(kFALSE);
    hXZ->Draw("COLZ");
    
    c->cd(box1_nb+2);
    c->cd(box1_nb+2)->SetLogz();
    hYZ->SetStats(kFALSE);
    hYZ->Draw("COLZ");
}


#endif
