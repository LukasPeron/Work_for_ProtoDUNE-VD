// THIS CODE IS IN C LANGUAGE //

R__ADD_INCLUDE_PATH("lardataobj/RawData/RawDigit.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Wire.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Hit.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Shower.h")
R__ADD_INCLUDE_PATH("gallery/Event.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/MCSFitResult.h")
#include "import.h"
#include "tools.h"
#include <vector>
#include <cmath>
#include <iostream>
#include "TGraph.h"
#include <string>
#include "TVector3.h"

using namespace art;
//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//Shower Analysis

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

double AngleDiffAxe(vector<vector<double>> Pts, vector<vector<double>> SimPts, vector<double> O)
{
    vector<double> v = UnitVectShowerAxe(SimPts,O);
    vector<double> u = UnitVectShowerAxe(Pts,O);
    double cosa = scalar(v,u);
    return acos(cosa);
}

double ShowerMAXLength(vector<vector<double>> Pts, vector<double> O)
{
    vector<double> u = UnitVectShowerAxe(Pts,O);
    vector<vector<double>> AxPts = XYZ_to_rl(Pts,O,u);
    double imax = std::max_element(AxPts.at(0).begin(),AxPts.at(0).end()) - AxPts.at(0).begin();
    double lmax = AxPts.at(0).at(imax);
    return lmax;
}

double ShowerMAXRadius(vector<vector<double>> Pts, vector<double> O)
{
    vector<double> u = UnitVectShowerAxe(Pts,O);
    vector<vector<double>> AxPts = XYZ_to_rl(Pts,O,u);
    double imax = std::max_element(AxPts.at(1).begin(),AxPts.at(1).end()) - AxPts.at(1).begin();
    double rmax = AxPts.at(1).at(imax);
    return rmax;
}

double Energy_Cone(vector<vector<double>> Pts, vector<double> O, double a /*angle en rad*/)
{
    vector<double> v = UnitVectShowerAxe(Pts,O);
    vector<vector<double>> OMPts = XYZ_to_OM(Pts,O);
    double cosa = cos(a);
    double ECone(0);
    int N = Pts.at(0).size();
    for (int i=0; i<N; i++)
    {
        vector<double> OM = OMPts.at(i);
        double r = sqrt(scalar(OM,OM));
        double costheta = scalar(OM,v)/r;
        if (costheta>cosa)
        {
            ECone+=Pts.at(3).at(i);
        }
    }
    return ECone;
}

void plot_AngularDistrib(vector<vector<double>> Pts, vector<double> O, string plotlabel,
                             int NbAngles = 50, double a_min = 0, double a_max = 90)
{
double A[NbAngles], E[NbAngles], Eint[NbAngles];
    for (int i=0; i<NbAngles; i++)
    {
        double a = a_min + (i*(a_max-a_min)/(NbAngles-1));
        A[i] = a;
        Eint[i] = Energy_Cone(Pts,O,a*3.1416/180);
        if (i==0)
        {
            E[i] = Eint[i];
        }
        else
        {
            E[i] = Eint[i]-Eint[i-1];
        }
    }
    string clabel = "ca" + plotlabel;
    TCanvas* c = new TCanvas(clabel.c_str(),(plotlabel+" Distriution angulaire de l'energie").c_str());
    c->Divide(2);
    
    c->cd(1);
    TGraph* g1 = new TGraph(NbAngles,A,Eint);
    g1->SetTitle((plotlabel+" Distriution angulaire integree;angle(°);E").c_str());
    g1->Draw("AC*");
    
    c->cd(2);
    TGraph* g2 = new TGraph(NbAngles,A,E);
    g2->SetTitle((plotlabel+" Distriution angulaire;angle(°);E").c_str());
    g2->Draw("AC*");
}

double ShowerConeAngle(vector<vector<double>> Pts, vector<double> O, double EnergyPart=0.95, double angle_eps = 0.1)
{
    int N=Pts.at(3).size();
    double Wtot(0);
    for (int i = 0; i < N; i++) 
    {
        Wtot+=Pts.at(3).at(i);
    }
    double Ethreshold  = Wtot * EnergyPart;
    double ECone(0);
    double a(0);
    while (ECone<Ethreshold)
    {
        a+=angle_eps;
        ECone = Energy_Cone(Pts,O,a*3.1416/180);
    }
    return a;
}

double Energy_Cyl(vector<vector<double>> Pts, vector<double> O, double r /*rayon (cm)*/, double l /*longueur du cylindre depuis l'origine(cm)*/)
{
    vector<double> v = UnitVectShowerAxe(Pts,O);
    vector<vector<double>> AxPts = XYZ_to_rl(Pts,O,v);
    double ECyl(0);
    int N = Pts.at(0).size();
    for (int i=0; i<N; i++)
    {
        if (AxPts.at(1).at(i)<r && AxPts.at(0).at(i)<l)
        {
            ECyl+=Pts.at(3).at(i);
        }
    }
    return ECyl;
}

vector<double> X0_and_RM(vector<vector<double>> SimPts, vector<double> O, string output_RM="", string output_X0="", int evnt=1, 
               string label="", double r_step=0.1, double l_step=0.5)
{
    int N = SimPts.at(0).size();
    double rmax = ShowerMAXRadius(SimPts,O);
    double lmax = ShowerMAXLength(SimPts,O);
    double Etot(0);
    int Nbin_X0 = lmax/l_step;
    int Nbin_RM = rmax/r_step;
    // // string X0_label = "X0_"+std::to_string(evnt+1);
    // string RM_label = "RM_"+std::to_string(evnt+1);
    // string plot_label_X0 = "Normalized energy versus the shower length event #"+std::to_string(evnt+1);
    // string plot_label_RM = "Normalized energy versus the shower radius event #"+std::to_string(evnt+1);
    // TH1D *plot_X0 = new TH1D(X0_label.c_str(), plot_label_X0.c_str(), Nbin_X0, 0, lmax);
    // TH1D *plot_RM = new TH1D(RM_label.c_str(), plot_label_RM.c_str(), Nbin_RM, 0, rmax);

    // ofstream Flux1(output_RM.c_str(), ios::app);
    // ofstream Flux2(output_X0.c_str(), ios::app);

    for (int i = 0; i < N; i++) 
    {
        Etot+=SimPts.at(3).at(i); //énergie totale
    }
    cout << "Energie totale : " << Etot << endl;
    double l(0);
    double r(0);
    bool test_RM= true;
    bool test_X0 = true;
    double X0(0);
    double RM(0);
    for (int i=0; i<Nbin_X0; i++)
    {
        double E_X0 = Energy_Cyl(SimPts, O, rmax, l);
        // plot_X0->AddBinContent(i, E_X0/Etot);
        if(E_X0/Etot >= (1/std::exp(1.0)) && test_X0)
        {
            X0 = l;
            cout << "Radiation lenght : " << X0 << " cm" << endl;
            test_X0 = false;
            break;
        }
        l+=l_step;
    }

    for (int i=0; i<Nbin_RM; i++)
    {
        double E_RM=Energy_Cyl(SimPts, O, r, lmax);
        // plot_RM->AddBinContent(i, E_RM/Etot);
        if(E_RM/Etot >= 0.9 && test_RM)
        {
            RM = r;
            cout << "Molière Radius : " << RM << " cm" << endl;
            test_RM = false;
            break;
        }
        r+=r_step;
    }

    // Flux1 << RM << ";" << evnt+1 << endl;
    // Flux2 << X0 << ";" << evnt+1 << endl;

    // TCanvas *c1 = new TCanvas(X0_label.c_str());
    // THStack *hs = new THStack("hs", "");
    // TLine *X0_line = new TLine(X0, 0, X0, 1);
    // plot_X0->SetStats(kFALSE);
    // plot_X0->SetLineColor(kRed+2);
    // plot_X0->SetXTitle("Length from the origin [cm]");
    // plot_X0->SetYTitle("E/E_tot");
    // plot_X0->SetLineWidth(2);
    // X0_line->SetLineWidth(2);
    // // plot_X0->Draw();
    // // X0_line->Draw();
    // // string save_name_X0_1 = "/silver/DUNE/peron/simu/radiation_length/electrons/sim/better_calculation/root/"+X0_label+"_"+label+".root";
    // // string save_name_X0_2 = "/silver/DUNE/peron/simu/radiation_length/electrons/sim/better_calculation/svg/"+X0_label+"_"+label+".svg";
    // // string save_name_X0_3 = "/silver/DUNE/peron/simu/radiation_length/electrons/sim/better_calculation/png/"+X0_label+"_"+label+".png";
    // // c1->SaveAs(save_name_X0_1.c_str());
    // // c1->SaveAs(save_name_X0_2.c_str());
    // // c1->SaveAs(save_name_X0_3.c_str());

    // TCanvas *c2 = new TCanvas(RM_label.c_str());
    // THStack *hs2 = new THStack("hs2", "");
    // TLine *RM_line = new TLine(RM, 0, RM, 1);
    // plot_RM->SetStats(kFALSE);
    // plot_RM->SetLineColor(kRed+2);
    // plot_RM->SetXTitle("Radius from the origin [cm]");
    // plot_RM->SetYTitle("E/E_tot");
    // plot_RM->SetLineWidth(2);
    // RM_line->SetLineWidth(2);
    // // plot_RM->Draw();
    // // RM_line->Draw();
    // string save_name_RM_1 = "/silver/DUNE/peron/simu/moliere_radius/electrons/sim/better_calculation/root/"+RM_label+"_"+label+".root";
    // string save_name_RM_2 = "/silver/DUNE/peron/simu/moliere_radius/electrons/sim/better_calculation/svg/"+RM_label+"_"+label+".svg";
    // string save_name_RM_3 = "/silver/DUNE/peron/simu/moliere_radius/electrons/sim/better_calculation/png/"+RM_label+"_"+label+".png";
    // c2->SaveAs(save_name_RM_1.c_str());
    // c2->SaveAs(save_name_RM_2.c_str());
    // c2->SaveAs(save_name_RM_3.c_str());
    
    vector<double> RM_XO_result = {RM, X0};

    return RM_XO_result;
}

vector<vector<vector<double>>> Selec_points_RM_X0(vector<vector<double>> Pts, vector<vector<double>> SimPts, vector<double> O, double RM, double X0){
    
    vector<double> v_Pts = UnitVectShowerAxe(Pts,O);
    vector<double> v_SimPts = UnitVectShowerAxe(SimPts,O);

    vector<vector<double>> Allign_Pts = XYZ_to_rl(Pts, O, v_Pts);
    vector<vector<double>> Allign_SimPts = XYZ_to_rl(SimPts, O, v_SimPts);
    
    vector<vector<double>> Selected_Pts(4);
    vector<vector<double>> Selected_SimPts(4);
    
    int Nhits = Allign_Pts.at(0).size();
    int NSimHits = Allign_SimPts.at(0).size();

    cout << Nhits << endl;
    cout << NSimHits << endl;
    for (int i=0; i < Nhits; i++){
        if (Allign_Pts.at(0).at(i)<=X0 && Allign_Pts.at(1).at(i)<=RM){
            Selected_Pts[0].push_back(Pts.at(0).at(i));
            Selected_Pts[1].push_back(Pts.at(1).at(i));
            Selected_Pts[2].push_back(Pts.at(2).at(i));
            Selected_Pts[3].push_back(Pts.at(3).at(i));
        }
    }
    for (int j=0; j < NSimHits; j++){
        if (Allign_SimPts.at(0).at(j)<=X0 && Allign_SimPts.at(1).at(j)<=RM){
            Selected_SimPts[0].push_back(SimPts.at(0).at(j));
            Selected_SimPts[1].push_back(SimPts.at(1).at(j));
            Selected_SimPts[2].push_back(SimPts.at(2).at(j));
            Selected_SimPts[3].push_back(SimPts.at(3).at(j));
        }
    }

    vector<vector<vector<double>>> res = {Selected_Pts, Selected_SimPts};

    return res;

}

void Viewer3D(
    std::vector<std::vector<double>> Pts, //List of Reco Deposit Coordinates and Charge
    std::vector<std::vector<double>> SimPts, //List of Sim Deposit Coordinates and Charge
    string plotlabel //Name to distinguish different plots
    ){  
    
    string clabel = "c3d_" + plotlabel;
    TCanvas *c3d = new TCanvas(clabel.c_str(),"Reco3D");
    c3d->cd();
    
    //Generate axis
    TH3F *axes = new TH3F("axes", "Reco 3D; X(cm); Y(cm); Z(cm)", 1, 0, 400, 1, -400, 400, 1, -400, 400);
    axes->SetDirectory(0);
    axes->GetXaxis()->SetTitleOffset(2);
    axes->GetYaxis()->SetTitleOffset(2);
    axes->GetZaxis()->SetTitleOffset(1.5);
    std::string Labels = plotlabel + " : simulated and reconstructed SpacePoints";
    axes->SetTitle(Labels.c_str());
    axes->SetStats(kFALSE);
    axes->Draw();
    
    //Reco Deposit plot (cf tools.h)
    TPolyMarker3D* Poly = Make_Poly_3D(Pts);   
    Poly->SetMarkerColor(kBlue);
    Poly->Draw();
    
    //Sim Deposit plot (cf tools.h)
    TPolyMarker3D* SimPoly = Make_Poly_3D(SimPts);  
    SimPoly->SetMarkerColor(kRed);
    SimPoly->Draw();
    
    //Legend
    auto legend = new TLegend();
    legend->SetHeader("Legend","C");
    legend->AddEntry(Poly,"Reco Energy Deposit","P");
    legend->AddEntry(SimPoly,"Sim Energy Deposit","P");
    legend->SetBorderSize(1);
    legend->Draw();
}


double MoliereRadius(vector<vector<double>> Pts, vector<double> O, double EnergyPart=0.90, double r_step = 0.1)
{
    int N=Pts.at(3).size();
    double lmax = ShowerMAXLength(Pts,O);
    double Wtot(0);
    for (int i = 0; i < N; i++) 
    {
        Wtot+=Pts.at(3).at(i);
    }
    double Ethreshold  = Wtot * EnergyPart;
    double ECyl(0);
    double r(0);
    while (ECyl<Ethreshold)
    {
        r+=r_step;
        ECyl = Energy_Cyl(Pts,O,r,lmax);
    }
    return r;
}

double ShowerCylLength(vector<vector<double>> Pts, vector<double> O, double EnergyPart=0.6667, double l_step = 1)
{
    int N=Pts.at(3).size();
    double rmax = ShowerMAXRadius(Pts,O);
    double Wtot(0);
    for (int i = 0; i < N; i++) 
    {
        Wtot+=Pts.at(3).at(i);
    }
    double Ethreshold  = Wtot * EnergyPart;
    double ECyl(0);
    double l(0);
    while (ECyl<Ethreshold)
    {
        l+=l_step;
        ECyl = Energy_Cyl(Pts,O,rmax,l);
    }
    return l;

}

void Slicing(vector<vector<double>> Pts, vector<vector<double>> SimPts, vector<double> O, string plotlabel, bool ShowCylModelLimits = true, int NBin = 100)
{
    //Useful data
    
    vector<double> v = UnitVectShowerAxe(SimPts,O);
    std::vector<double> Charges = Pts.at(3);
    std::vector<double> Energies = SimPts.at(3);
    vector<vector<double>> AxPts = XYZ_to_rl(Pts,O,v);
    vector<vector<double>> AxSimPts = XYZ_to_rl(SimPts,O,v);
    double lmax = ShowerMAXLength(SimPts,O);
    double rmax = ShowerMAXRadius(SimPts,O);

    int N = Pts.at(0).size();
    int Nsim = SimPts.at(0).size();
    
    double Qtot(0);
    for (int i = 0; i < N; i++) 
    {
        Qtot+=Charges.at(i);
    }
    
    double Etot(0);
    for (int i = 0; i < Nsim; i++) 
    {
        Etot+=Energies.at(i);
    }
    
    //Initial objects
    TCanvas* c = new TCanvas(("c_Slicing_"+plotlabel).c_str(),("Slicing "+plotlabel).c_str());
    c->Divide(2,2);
    
    auto legend1 = new TLegend(0,0,1,1);
    legend1->SetHeader("Legend - Energy Distribution along Shower Axis","C");
    auto legend2 = new TLegend(0,0,1,1);
    legend2->SetHeader("Legend - Transverse Energy Distribution","C");
    
    string h1label = "hle_" + plotlabel;
    string h2label = "hlq_" + plotlabel;
    string h3label = "hre_" + plotlabel;
    string h4label = "hrq_" + plotlabel;
    
    TH1D *hle = new TH1D(h1label.c_str(), "Distribution along Shower Axis;l(cm);E/Etot",NBin,0,lmax);
    legend1->AddEntry(hle,"Energy Distribution along Shower Axis","L");
    TH1D *hlq = new TH1D(h2label.c_str(), "Charge Distribution along Shower Axis;l(cm);Q/Qtot",NBin,0,lmax);
    legend1->AddEntry(hlq,"Charge Distribution along Shower Axis","L");
    TH1D *hre = new TH1D(h3label.c_str(), "Global Transverse Energy Distribution;r(cm);E/Etot",NBin,0,rmax);
    legend2->AddEntry(hre,"Global Transverse Energy Distribution","L");
    TH1D *hrq = new TH1D(h4label.c_str(), "Global Transverse Charge Distribution;r(cm);Q/Qtot",NBin,0,rmax);
    legend2->AddEntry(hrq,"Global Transverse Charge Distribution","L");
    
    for (int i=0; i<N; i++)
    {
        double l = AxPts.at(0).at(i);
        double r = AxPts.at(1).at(i);
        double Q = Pts.at(3).at(i) / Qtot;
        hlq->Fill(l,Q);
        hrq->Fill(r,Q);
    }
    
    for (int i=0; i<Nsim; i++)
    {
        double l = AxSimPts.at(0).at(i);
        double r = AxSimPts.at(1).at(i);
        double E = SimPts.at(3).at(i) / Etot;
        hle->Fill(l,E);
        hre->Fill(r,E);
    }
    
    
    c->cd(1);
    hle->SetStats(kFALSE);
    hle->SetLineColor(kRed+2);
    hle->Draw("HIST");
    
    hlq->SetStats(kFALSE);
    hlq->SetLineColor(kCyan+2);
    hlq->Draw("same HIST");
    
    c->cd(2);
    hre->SetStats(kFALSE);
    hre->SetLineColor(kRed+2);
    hre->Draw("HIST");
    
    hrq->SetStats(kFALSE);
    hrq->SetLineColor(kCyan+2);
    hrq->Draw("same HIST");
    
    if (ShowCylModelLimits)
    {
        double rsim = MoliereRadius(SimPts,O);
        double rreco = MoliereRadius(Pts,O);
        double lsim = ShowerCylLength(SimPts,O,rsim);
        double lreco = ShowerCylLength(Pts,O,rreco);
        TLine *lre = new TLine(rsim,0,rsim,hre->GetMaximum());
        TLine *lrq = new TLine(rreco,0,rreco,hre->GetMaximum());
        TLine *lle = new TLine(lsim,0,lsim,hle->GetMaximum());
        TLine *llq = new TLine(lreco,0,lreco,hle->GetMaximum());
        
        cout << endl << "Rayon de Molière Sim : "<< rsim <<"cm";
        cout << endl << "Rayon de Molière Reco : "<< rreco <<"cm";
        cout << endl << "Longueur de Radiation Sim : "<< lsim <<"cm";
        cout << endl << "Longueur de Radiation Reco : "<< lreco <<"cm";
        
        c->cd(1);
        c->cd(2)->SetLogy();
        
        lle->SetLineColor(2);
        legend1->AddEntry(lle,"Radiation Length (sim data)","L");
        lle->Draw("same");
        
        llq->SetLineColor(kBlue);
        legend1->AddEntry(llq,"Radiation Length (reco data)","L");
        llq->Draw("same");
        
        c->cd(2);
        c->cd(2)->SetLogy();
        
        lre->SetLineColor(2);
        legend2->AddEntry(lre,"Moliere Radius (sim data)","L");
        lre->Draw("same");
        
        lrq->SetLineColor(kBlue);
        legend2->AddEntry(lrq,"Moliere Radius (reco data)","L");
        lrq->Draw("same");
        
    }
    
    c->cd(3);
    legend1->SetBorderSize(1);
    legend1->Draw();
    
    c->cd(4);
    legend2->SetBorderSize(1);
    legend2->Draw();
    
    c->Update();
}

//---------------------------------------------------------------------------------------------

void Slicing_allEV(string path, string label, double EVmax, string plotlabel, double DistanceUncertainty=4, int NBin = 100)
{
    //Initialisation
    double rsim(0);
    double rreco(0);
    double lsim(0);
    double lreco(0);
    
    //Limit
    double rmax(0);
    double lmax(0);
    for (int evn = 0; evn<EVmax; evn++)
    {
        std::vector<std::vector<double>> SimPts = XYZE(path, label, "_g4_stage2.root", evn);
        std::vector<double> O = Origin(label,evn);
        double lmaxev = ShowerMAXLength(SimPts,O);
        double rmaxev = ShowerMAXRadius(SimPts,O);
        if (lmaxev>lmax){lmax=lmaxev;}
        if (rmaxev>rmax){rmax=rmaxev;}
    }
    
        //Initial objects
    TCanvas* c = new TCanvas(("c_Slicing_"+plotlabel).c_str(),("Slicing "+plotlabel).c_str());
    c->Divide(2);
    
    auto legend1 = new TLegend();
    legend1->SetHeader("Legend - Energy Distribution along Shower Axis","C");
    auto legend2 = new TLegend();
    legend2->SetHeader("Legend - Transverse Energy Distribution","C");
    
    string h1label = "hle_" + plotlabel;
    string h2label = "hlq_" + plotlabel;
    string h3label = "hre_" + plotlabel;
    string h4label = "hrq_" + plotlabel;
    
    TH1D *hle = new TH1D(h1label.c_str(), "Distribution along Shower Axis;l(cm);E/Etot",NBin,0,lmax);
    legend1->AddEntry(hle,"Energy Distribution along Shower Axis","L");
    TH1D *hlq = new TH1D(h2label.c_str(), "Charge Distribution along Shower Axis;l(cm);Q/Qtot",NBin,0,lmax);
    legend1->AddEntry(hlq,"Charge Distribution along Shower Axis","L");
    TH1D *hre = new TH1D(h3label.c_str(), "Global Transverse Energy Distribution;r(cm);E/Etot",NBin,0,rmax);
    legend2->AddEntry(hre,"Global Transverse Energy Distribution","L");
    TH1D *hrq = new TH1D(h4label.c_str(), "Global Transverse Charge Distribution;r(cm);Q/Qtot",NBin,0,rmax);
    legend2->AddEntry(hrq,"Global Transverse Charge Distribution","L");
    
    for (int evn = 0; evn<EVmax; evn++)
    {
        string ext = "_reco_3sigma.root";
        string evlabel = plotlabel + "_" + to_string(evn+1);
        cout << endl << endl << "-----------------  EVENT #" << evn+1 << "  ------------------------"<< endl;
        std::vector<std::vector<double>> Pts = XYZQ(path, label, ext, evn);
        std::vector<std::vector<double>> SimPts = XYZE(path, label, "_g4_stage2.root", evn);
        std::vector<double> O = Origin(label,evn);
    
        //Useful data
        
        vector<double> v = UnitVectShowerAxe(SimPts,O);
        std::vector<double> Charges = Pts.at(3);
        std::vector<double> Energies = SimPts.at(3);
        vector<vector<double>> AxPts = XYZ_to_rl(Pts,O,v);
        vector<vector<double>> AxSimPts = XYZ_to_rl(SimPts,O,v);
        
        double lmaxev = ShowerMAXLength(SimPts,O);
        double rmaxev = ShowerMAXRadius(SimPts,O);
        if (lmaxev>lmax){lmax=lmaxev;}
        if (rmaxev>rmax){rmax=rmaxev;}
        
        int N = Pts.at(0).size();
        int Nsim = SimPts.at(0).size();
        
        double Qtot(0);
        for (int i = 0; i < N; i++) 
        {
            Qtot+=Charges.at(i);
        }
        
        double Etot(0);
        for (int i = 0; i < Nsim; i++) 
        {
            Etot+=Energies.at(i);
        }
        
        
        for (int i=0; i<N; i++)
        {
            double l = AxPts.at(0).at(i);
            double r = AxPts.at(1).at(i);
            double Q = Pts.at(3).at(i) / Qtot;
            hlq->Fill(l,Q);
            hrq->Fill(r,Q);
        }
        
        for (int i=0; i<Nsim; i++)
        {
            double l = AxSimPts.at(0).at(i);
            double r = AxSimPts.at(1).at(i);
            double E = SimPts.at(3).at(i) / Etot;
            hle->Fill(l,E);
            hre->Fill(r,E);
        }
        
        rsim += MoliereRadius(SimPts,O);
        rreco += MoliereRadius(Pts,O);
        lsim += ShowerCylLength(SimPts,O,rsim);
        lreco += ShowerCylLength(Pts,O,rreco);
        
    }
    
    rsim /= EVmax;
    rreco /= EVmax;
    lsim /= EVmax;
    lreco /= EVmax;
    
    double Emaxl = max(hle->GetMaximum(),hlq->GetMaximum());
    double Emaxr = max(hre->GetMaximum(),hrq->GetMaximum());
    
    c->cd(1);
    hle->SetStats(kFALSE);
    hle->SetLineColor(kRed+2);
    hle->GetYaxis()->SetLimits(0,Emaxl);
    
    hlq->SetStats(kFALSE);
    hlq->SetLineColor(kCyan+2);
    
    
    hlq->Draw("HIST");
    hle->Draw("same HIST");
    
    c->cd(2);
    hre->SetStats(kFALSE);
    hre->SetLineColor(kRed+2);
    hre->GetYaxis()->SetLimits(0,Emaxr);
    
    hrq->SetStats(kFALSE);
    hrq->SetLineColor(kCyan+2);

    hrq->Draw("HIST");
    hre->Draw("same HIST");
    
    TLine *lre = new TLine(rsim,0,rsim,Emaxr);
    TLine *lrq = new TLine(rreco,0,rreco,Emaxr);
    TLine *lle = new TLine(lsim,0,lsim,Emaxl);
    TLine *llq = new TLine(lreco,0,lreco,Emaxl);
    
    double Xul[]={lreco-DistanceUncertainty,lreco-DistanceUncertainty,lreco+DistanceUncertainty,lreco+DistanceUncertainty};
    double Yul[]={0,Emaxl,Emaxl,0};
    TGraph *ul = new TGraph(4,Xul,Yul);
    
    double Xur[]={rreco-DistanceUncertainty,rreco-DistanceUncertainty,rreco+DistanceUncertainty,rreco+DistanceUncertainty};
    double Yur[]={0,Emaxr,Emaxr,0};
    TGraph *ur = new TGraph(4,Xur,Yur);
    
    cout << endl << "Rayon de Molière Sim : "<< rsim <<"cm";
    cout << endl << "Rayon de Molière Reco : "<< rreco <<"cm";
    cout << endl << "Longueur de Radiation Sim : "<< lsim <<"cm";
    cout << endl << "Longueur de Radiation Reco : "<< lreco <<"cm";
        
    c->cd(1);
    c->cd(1)->SetLogy();
    
    lle->SetLineColor(2);
    legend1->AddEntry(lle,"Radiation Length (sim data)","L");
    lle->Draw("same");
    
    llq->SetLineColor(kBlue);
    legend1->AddEntry(llq,"Radiation Length (reco data)","L");
    llq->Draw("same");
    
    ul->SetFillColorAlpha(kBlue,0.1);
    ul->SetFillStyle(3003);
    ul->Draw("same F");
    
    c->cd(2);
    c->cd(2)->SetLogy();
    
    lre->SetLineColor(2);
    legend2->AddEntry(lre,"Moliere Radius (sim data)","L");
    lre->Draw("same");
    
    lrq->SetLineColor(kBlue);
    legend2->AddEntry(lrq,"Moliere Radius (reco data)","L");
    lrq->Draw("same");
    
    ur->SetFillColorAlpha(kBlue,0.3);
    ur->SetFillStyle(3003);
    ur->Draw("same F");
    
    c->cd(1);
    legend1->SetBorderSize(1);
    legend1->Draw();
    
    c->cd(2);
    legend2->SetBorderSize(1);
    legend2->Draw();
    
    c->Update();
}

//----------------------------------------------------------------------------------------------------

void Projections_SimReco(
    std::vector<std::vector<double>> Pts, //List of Reco Deposit Coordinates and Charge
    std::vector<std::vector<double>> SimPts, //List of Sim Deposit Coordinates and Charge
    string plotlabel, //Name to distinguish different plots
    std::vector<double> Origin, //Origin Coordinates to build an axis
    bool OnlyDeposit = false, //if true, no other info than Energy deposit will be plot
    double EminPeak = 0
     ){
    
    //Initialisation
    
    string clabel = "cB_" + plotlabel;
    TCanvas *cB = new TCanvas(clabel.c_str(),"Sim-Reco projections");
    cB->cd();
    cB->Divide(2,2);
    
    TMultiGraph* mgXY = new TMultiGraph();
    TMultiGraph* mgYZ = new TMultiGraph();
    TMultiGraph* mgXZ = new TMultiGraph();
    
    auto legend = new TLegend(0,0,1,1);
    legend->SetHeader("Legend","C");
    legend->SetBorderSize(0);
    
    //Barycenters (cf tools.h)
    std::vector<double> Charges = Pts.at(3);
    std::vector<double> Energies = SimPts.at(3);
    std::vector<double> B = Barycenter(Pts, Charges);
    std::vector<double> Bs = Barycenter(SimPts, Energies);
    
    //Energy Max
    double imax = std::max_element(Energies.begin(),Energies.end()) - Energies.begin();
    double Emax = Energies.at(imax);
    
    //Energy Peaks (over a certain value EminPeak) 
    //-> only if asked by giving a positive value to EminPeak
    std::vector<std::vector<double>> Peaks(3);
    if (EminPeak>0){
        int npk(0);
        for (int i=0 , szs = SimPts.size(); i<szs; i++)
        {
            double E = Energies.at(i);
            if (E>EminPeak)
            {
                Peaks.at(0).push_back(SimPts.at(0).at(i));
                Peaks.at(1).push_back(SimPts.at(1).at(i));
                Peaks.at(2).push_back(SimPts.at(2).at(i));
                npk++;
            }
        }
    }
    
    //Deposit plot
    plot_proj_pts(Pts,mgXY,mgYZ,mgXZ,38,6,legend,"Energy Depot (Reco)");
    plot_proj_pts(SimPts,mgXY,mgYZ,mgXZ,45,6,legend,"Energy Depot (Sim)");
    
    //Specific points plot
    if (!(OnlyDeposit))
    {
        plot_proj_pt(B,mgXY,mgYZ,mgXZ,4,24,legend,"Barycenter (Reco)");
        plot_proj_pt(Bs,mgXY,mgYZ,mgXZ,2,24,legend,"Barycenter (Sim)");
        plot_proj_pt(Origin,mgXY,mgYZ,mgXZ,3,29,legend,"Beam Origin");
        plot_proj_line_2pts(Origin,B,mgXY,mgYZ,mgXZ,4,legend,"Shower axis (Reco)");
        plot_proj_line_2pts(Origin,Bs,mgXY,mgYZ,mgXZ,2,legend,"Shower axis (Sim)");
        if (EminPeak>0){plot_proj_pts(Peaks,mgXY,mgYZ,mgXZ,2,25,legend,"Energy Peak Sim");}
    }
    
    //Drawing the 2D projections
    cB->cd(3);
    string labelXY = plotlabel +" : Y(X);X(cm);Y(cm)";
    mgXY->SetTitle(labelXY.c_str());
    mgXY->Draw("A p l");
    
    cB->cd(2);
    string labelYZ = plotlabel +" : Z(Y);Y(cm);Z(cm)";
    mgYZ->SetTitle(labelYZ.c_str());
    mgYZ->Draw("A p l");

    
    cB->cd(1);
    string labelXZ = plotlabel +" : Z(X);X(cm);Z(cm)";
    mgXZ->SetTitle(labelXZ.c_str());
    mgXZ->Draw("A p l");

    
    cB->cd(4);
    legend->Draw();
    
    cout << endl << "Distance between the Barycenters : " << Distance(B,Bs) << " cm" << endl;
    cout << "Emax = " << Emax << endl;
}


void Geometry_tests(string path, string label, string ext)
{
    string extension_g4_stage2 = "_g4_stage2.root";
    int EVmin = 1;
    int EVmax = 10;
    string output_file_RM_reco = "/silver/DUNE/peron/simu/moliere_radius/electrons/reco/better_calculation/" + label + "3_sigma_MR_G4_v2.txt";
    string output_file_X0_reco = "/silver/DUNE/peron/simu/radiation_length/electrons/reco/better_calculation/" + label + "3_sigma_RL_G4_v2.txt";
    string output_file_RM_sim = "/silver/DUNE/peron/simu/moliere_radius/electrons/sim/better_calculation/" + label + "3_sigma_MR_G4_v2.txt";
    string output_file_X0_sim = "/silver/DUNE/peron/simu/radiation_length/electrons/sim/better_calculation/" + label + "3_sigma_RL_G4_v2.txt";
    // ofstream Flux1(output_file_RM_reco.c_str());
    // ofstream Flux2(output_file_X0_reco.c_str());
    // ofstream Flux3(output_file_RM_sim.c_str());
    // ofstream Flux4(output_file_X0_sim.c_str());
    // Flux1 << "RM ; event" << endl;
    // Flux2 << "X0 ; event" << endl;
    // Flux3 << "RM ; event" << endl;
    // Flux4 << "X0 ; event" << endl;
    for (int evn = EVmin-1; evn<EVmax; evn++)
    {
        string evlabel = label + "_" + to_string(evn+1);
        cout << endl << endl << "-----------------  EVENT #" << evn+1 << "  ------------------------"<< endl;
        std::vector<std::vector<double>> Pts = XYZQ(path, label, ext, evn);
        std::vector<std::vector<double>> SimPts = XYZE(path, label, extension_g4_stage2, evn);
        std::vector<double> O = Origin(label,evn);
        Pts = ShiftXOpti_Pts(Pts,label,evn);
        // cout << endl << "Angle entre l'axe simulé et l'axe reconstruit : " << AngleDiffAxe(Pts, SimPts,O)*180/3.1416 << "° " << endl;
        // cout << endl << "Longueur de la gerbe simulée : " << ShowerMAXLength(SimPts, O) << " cm" << endl << "Longueur de la gerbe reconstruite : " << ShowerMAXLength(Pts, O) << " cm" << endl;
        // cout << endl << "Rayon de la gerbe simulée : " << ShowerMAXRadius(SimPts, O) << " cm" << endl << "Rayon de la gerbe reconstruite : " << ShowerMAXRadius(Pts, O) << " cm" << endl;

        // vector<double> RM_X0_Pts=X0_and_RM(Pts, O, output_file_RM_reco, output_file_X0_reco, evn, label);
        // vector<double> RM_X0_SimPts=X0_and_RM(SimPts, O, output_file_RM_sim, output_file_X0_sim, evn, label);
        // double RM_Pts = RM_X0_Pts.at(0);

        
        string rootfilelabel = path + label + ext;
        vector<string> filename(1, rootfilelabel);
        gallery::Event ev(filename); 

        recob::MCSFitResult FitResult;
        cout << "test" << endl;
        std::vector<float> radlength = FitResult.segmentRadLengths();
        cout << radlength.size() << endl;
        // double RM_SimPts = RM_X0_SimPts.at(0);
        // double X0_Pts = RM_X0_Pts.at(1);
        // double X0_SimPts = RM_X0_SimPts.at(1);

        // cout << "RM Simulation : " << RM_SimPts << " cm" << endl;
        // cout << "X0 Simulation : " << X0_SimPts << " cm" << endl;
        // cout << "RM Reconstruit : " << RM_Pts << " cm" << endl;
        // cout << "X0 Reconstruit : " << X0_Pts << " cm" << endl;

        // vector<vector<double>> Selected_Pts=Selec_points_RM_X0(Pts, SimPts, O, ShowerMAXRadius(Pts, O), ShowerMAXLength(Pts, O)).at(0);
        // vector<vector<double>> Selected_SimPts=Selec_points_RM_X0(Pts, SimPts, O, ShowerMAXRadius(SimPts, O), ShowerMAXLength(SimPts, O)).at(1);

        // cout << "Nombre de points simulés conservés : " << Selected_SimPts.at(0).size() << endl;
        // cout << "Nombre de points reconstruits conservés : " << Selected_Pts.at(0).size() << endl;

        // Viewer3D(Selected_Pts, Selected_SimPts, evlabel);
        // Projections_SimReco(Selected_Pts, Selected_SimPts, evlabel+" with Barycenters", O);

        // plot_AngularDistrib(Pts,O, evlabel,100,0,90);
        // plot_AngularDistrib(SimPts,O, "SIM "+evlabel,100,0,90);
        // cout << endl << "Angle du cone à "<< ep <<" d'énergie Sim : "<< ShowerConeAngle(SimPts,O,ep) <<"°";
        // cout << endl << "Angle du cone à "<< ep <<" d'énergie Reco : "<< ShowerConeAngle(Pts,O,ep) <<"°";

        
        // Slicing(Pts,SimPts,O,evlabel);
        // cout <<endl << endl << "***************************************************************" << endl;
    }
    // Slicing_allEV(label,10,label+"_allEvents");

    return 0;
}
