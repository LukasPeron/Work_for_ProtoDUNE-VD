// THIS CODE IS IN C LANGUAGE //

//Personal libraries
#include "import.h"
#include "tools.h"

//Standard libraries
#include <vector>
#include <cmath>
#include <iostream>
#include "TGraph.h"
#include <string>
#include <typeinfo>

using namespace art;

//Prototype for final "Data Files and Canvas Generator"

//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//3D View

//Plot both Sim and Reco SpacePoints in a 3D Box with dimensions of the detector

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

//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//2D Projections

//Plot both Sim and Reco SpacePoints on 2D projections with Barycenters, Origin of the Beam, Axis (and Energy Peaks over a certain value if asked)

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
    // std::vector<double> Charges = Pts.at(3);
    // std::vector<double> Energies = SimPts.at(3);
    // std::vector<double> B = Barycenter(Pts, Charges);
    // std::vector<double> Bs = Barycenter_simu(SimPts, Energies);
    
    //Energy Max
    // double imax = std::max_element(Energies.begin(),Energies.end()) - Energies.begin();
    // double Emax = Energies.at(imax);
    
    //Energy Peaks (over a certain value EminPeak) 
    //-> only if asked by giving a positive value to EminPeak
    // std::vector<std::vector<double>> Peaks(3);
    // if (EminPeak>0){
    //     int npk(0);
    //     for (int i=0 , szs = SimPts.size(); i<szs; i++)
    //     {
    //         double E = Energies.at(i);
    //         if (E>EminPeak)
    //         {
    //             Peaks.at(0).push_back(SimPts.at(0).at(i));
    //             Peaks.at(1).push_back(SimPts.at(1).at(i));
    //             Peaks.at(2).push_back(SimPts.at(2).at(i));
    //             npk++;
    //         }
    //     }
    // }
    
    //Deposit plot
    plot_proj_pts(Pts,mgXY,mgYZ,mgXZ,38,6,legend,"Energy Depot (Initiaux)"); //reco
    plot_proj_pts(SimPts,mgXY,mgYZ,mgXZ,45,6,legend,"Energy Depot (Selectionnes)"); //sim
    
    //Specific points plot
    if (!(OnlyDeposit))
    {
        // plot_proj_pt(B,mgXY,mgYZ,mgXZ,4,24,legend,"Barycenter (Reco)");
        // plot_proj_pt(Bs,mgXY,mgYZ,mgXZ,2,24,legend,"Barycenter (Sim)");
        // plot_proj_pt(Origin,mgXY,mgYZ,mgXZ,3,29,legend,"Beam Origin");
        // plot_proj_line_2pts(Origin,B,mgXY,mgYZ,mgXZ,4,legend,"Shower axis (Reco)");
        // plot_proj_line_2pts(Origin,Bs,mgXY,mgYZ,mgXZ,2,legend,"Shower axis (Sim)");
        // if (EminPeak>0){plot_proj_pts(Peaks,mgXY,mgYZ,mgXZ,2,25,legend,"Energy Peak Sim");}
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
    
    // cout << endl << "Distance between the Barycenters : " << Distance(B,Bs) << " cm" << endl;
    // cout << "Emax = " << Emax << endl;
}



void Projections_SimReco_compar(
    string save_path,
    string save_end_label,
    std::vector<std::vector<double>> Pts, //List of Reco Deposit Coordinates and Charge
    std::vector<std::vector<double>> Pts_pandora, //List of Reco Deposit Coordinates and Charge via pandora
    std::vector<std::vector<double>> SimPts, //List of Sim Deposit Coordinates and Charge
    std::vector<std::vector<double>> Track, //Pandora track
    string plotlabel, //Name to distinguish different plots
    std::vector<double> Origin, //Origin Coordinates to build an axis
    bool OnlyDeposit = false, //if true, no other info than Energy deposit will be plot
    double EminPeak = 0
     ){
    
    //Initialisation
    
    string clabel = "cB_" + plotlabel;
    TCanvas *cB = new TCanvas(clabel.c_str(),"Sim-Reco projections comparison");
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
    std::vector<double> Charges_pandora = Pts_pandora.at(3);
    std::vector<double> Energies = SimPts.at(3);
    // std::vector<double> B = Barycenter(Pts, Charges);
    // std::vector<double> Bs = Barycenter_simu(SimPts, Energies);
    
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
    plot_proj_pts(Pts,mgXY,mgYZ,mgXZ,4,6,legend,"Energy Depot (Reco3d)");
    plot_proj_pts(Pts_pandora,mgXY,mgYZ,mgXZ,3,6,legend,"Energy Depot (Pandora)");
    plot_proj_pts(SimPts,mgXY,mgYZ,mgXZ,2,6,legend,"Energy Depot (Sim)");
    plot_proj_pts(Track,mgXY,mgYZ,mgXZ,1,6,legend,"Track");
    // plot_proj_pts(Track_start,mgXY,mgYZ,mgXZ,6,6,legend,"Track start");
    // plot_proj_pts(Track_end,mgXY,mgYZ,mgXZ,6,6,legend,"Track end");
    
    //Specific points plot
    // if (!(OnlyDeposit))
    // {
    //     plot_proj_pt(B,mgXY,mgYZ,mgXZ,4,24,legend,"Barycenter (Reco)");
    //     plot_proj_pt(Bs,mgXY,mgYZ,mgXZ,2,24,legend,"Barycenter (Sim)");
    //     plot_proj_pt(Origin,mgXY,mgYZ,mgXZ,3,29,legend,"Beam Origin");
    //     plot_proj_line_2pts(Origin,B,mgXY,mgYZ,mgXZ,4,legend,"Shower axis (Reco)");
    //     plot_proj_line_2pts(Origin,Bs,mgXY,mgYZ,mgXZ,2,legend,"Shower axis (Sim)");
    //     if (EminPeak>0){plot_proj_pts(Peaks,mgXY,mgYZ,mgXZ,2,25,legend,"Energy Peak Sim");}
    // }
    
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

    // for (unsigned j = 0; j<Track.at(0).size(); j+=2){
    //     cB->cd(3);
    //     TLine* Track_yx = new TLine(Track.at(0).at(j), Track.at(1).at(j), Track.at(0).at(j+1), Track.at(1).at(j+1));
    //     Track_yx->SetLineColor(6);
    //     Track_yx->SetLineWidth(2);
    //     Track_yx->Draw();
    //     cB->cd(2);
    //     TLine* Track_zy = new TLine(Track.at(1).at(j), Track.at(2).at(j), Track.at(1).at(j+1), Track.at(2).at(j+1));
    //     Track_zy->SetLineColor(6);
    //     Track_zy->SetLineWidth(2);
    //     Track_zy->Draw();
    //     cB->cd(1);
    //     TLine* Track_zx = new TLine(Track.at(0).at(j), Track.at(2).at(j), Track.at(0).at(j+1), Track.at(2).at(j+1));
    //     Track_zx->SetLineColor(6);
    //     Track_zx->SetLineWidth(2);
    //     Track_zx->Draw();
    // }

    cB->cd(4);
    legend->Draw();
    
    string pdf_ext = ".pdf";
    string png_ext = ".png";
    string svg_ext = ".svg";

    string save_name = save_path+plotlabel+save_end_label+pdf_ext;
    string save_name2 = save_path+plotlabel+save_end_label+png_ext;
    string save_name3 = save_path+plotlabel+save_end_label+svg_ext;

    cB->SaveAs(save_name.c_str());
    cB->SaveAs(save_name2.c_str());
    cB->SaveAs(save_name3.c_str());
    
    // cout << endl << "Distance between the Barycenters : " << Distance(B,Bs) << " cm" << endl;
    // cout << "Emax = " << Emax << endl;
}






//----------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//Energy Distribution

//plot 2D projections of a colormap of the spatial Energy distribution in Bins of given size
void EnergyMap (
    std::vector<std::vector<double>> Pts, //List of Reco Deposit Coordinates and Charge
    std::vector<std::vector<double>> SimPts, //List of Sim Deposit Coordinates and Charge
    string plotlabel, //Name to distinguish different plots
    double BinSzX, double BinSzY, double BinSzZ, //Bins dimensions
    bool OptWindow = true, //Option to optimize the dimensions of the window for a best view of the shower
    
    //Default (else) : Detector dimensions
    //Note : à adapter au réel (y compris dans le tool associé)
    double xmin=0,double xmax=400, 
    double ymin = -400 , double ymax = 400, 
    double zmin = -400, double zmax = 400)
{
    string clabel = "cc_" + plotlabel;
    cout << "test 1" << clabel << endl;
    TCanvas *c1 = new TCanvas(clabel.c_str(),"Simulated - Reconstructed diff Energy");
    c1->Divide(3,2);
    
    //Optimized Window (cf tools.h)
    if (OptWindow)
    {
        std::vector<double> W = WindowSz(SimPts);
        xmin = W[0];
        xmax = W[1];
        ymin = W[2];
        ymax = W[3];
        zmin = W[4];
        zmax = W[5];
    }
    
    //plot Energy colormap 2D Histograms (cf tools.h)
    plot_proj_Hist(SimPts,c1,1,plotlabel+" - Normalized Energy (Simulation)",
                   BinSzX,BinSzY,BinSzZ,xmin,xmax,ymin,ymax,zmin,zmax);
    plot_proj_Hist(Pts,c1,4,plotlabel+" - Normalized Charge (Reconstruction)",
                   BinSzX,BinSzY,BinSzZ,xmin,xmax,ymin,ymax,zmin,zmax);
}

//---------------------------------------------------------------------------------------------

//Give lots of data in lists and plots about the energy map

// Note : Passer ce programme en création de fichier récapitulatif de chaque événement + générateur de Canvas spectre en énergie

void EnergyRecoData (
    std::vector<std::vector<double>> Pts, //List of Reco Deposit Coordinates and Charge
    std::vector<std::vector<double>> SimPts, //List of Sim Deposit Coordinates and Charge
    string plotlabel, //Name to distinguish different plots
    double BinSzX, double BinSzY, double BinSzZ, 
    double xmin=0,double xmax=400, 
    double ymin = -400 , double ymax = 400, 
    double zmin = -400, double zmax = 400,
    string label = "", int n_sigma = 0)
{
    std::vector<double> Charges = Pts.at(3);
    std::vector<double> Energies = SimPts.at(3);
    //Useful const
    int const N = Pts.at(0).size();
    int const Nsim = SimPts.at(0).size();
    int const NbinX = (xmax - xmin) / BinSzX;
    int const NbinY = (ymax - ymin) / BinSzY;
    int const NbinZ = (zmax - zmin) / BinSzZ;
    
    //
    cout << endl << "***" << endl;
    cout << "Nb of Simulated Depositions : " << Nsim << endl;
    cout << "Nb of Reconstructed Depositions : " << N << endl;
    
    cout << endl << "***" << endl << "Bins Info" << endl;
    cout << "X Direction : " << NbinX << "bins of size " << BinSzX << " cm (" << xmin << "cm < x < " << xmax << "cm)" << endl;
    cout << "Y Direction : " << NbinY << "bins of size " << BinSzY << " cm (" << ymin << "cm < y < " << ymax << "cm)" << endl;
    cout << "Z Direction : " << NbinZ << "bins of size " << BinSzZ << " cm (" << zmin << "cm < z < " << zmax << "cm)" << endl;
    
    cout << endl << "***" << endl << "Normalization process : " << endl;
    
    
    //Charge and Energy Integration for the Normalisation
    double Qtot(0);
    for (int i = 0; i < N; i++) 
    {
        Qtot+=Charges.at(i);
    }
    
    ofstream FluxCharge("/silver/DUNE/peron/simu/Total_charge.txt", ios::app);
    ofstream FluxHits("/silver/DUNE/peron/simu/Total_hits.txt", ios::app);

    if(FluxCharge && FluxHits)    
    {
        FluxCharge << label +" " << n_sigma << " sigma total charge : " << Qtot << " (méthode SummedADC)" << endl;
        FluxHits << label +" " << n_sigma << " sigma total hits : " << N << " (méthode SummedADC)" << endl;
    }
    else
    {
        cout << "ERROR: Unable to open the files." << endl;
    }

    double Etot(0);
        for (int i = 0; i < Nsim; i++) 
    {
        Etot+=Energies.at(i);
    }
    
    //
    cout << "Integrated Charge (Reco) : " << Qtot << " ADC" << endl;
    cout << "Integrated Energy (Sim) : " << Etot << " MeV" << endl;

    cout << endl << "***" << endl << "Filling Bins..." << endl;
    
    
    //Objects Initialisation
    string clabel = "c_" + plotlabel;
    TCanvas *c1 = new TCanvas(clabel.c_str(),"Simulated - Reconstructed diff Energy");
    
    string hdifflabel = "hdiff_" + plotlabel;
    TH3D *hdiff = new TH3D(hdifflabel.c_str(), "Energy diff;X(cm);Y(cm);Z(cm)", NbinX, xmin, xmax, NbinY, ymin,ymax , NbinZ ,zmin , zmax);
    string helabel = "he_" + plotlabel;
    TH3D *he = new TH3D(helabel.c_str(), "Energy;X(cm);Y(cm);Z(cm)", NbinX, xmin, xmax, NbinY, ymin,ymax , NbinZ ,zmin , zmax);
    string hqlabel = "hq_" + plotlabel;
    TH3D *hq = new TH3D(hqlabel.c_str(), "Charge;X(cm);Y(cm);Z(cm)", NbinX, xmin, xmax, NbinY, ymin,ymax , NbinZ ,zmin , zmax);
    
    //Initialisation of vectors of bin IDs (for comparison of filled bins)
    vector<int> QFilledBins;
    vector<int> EFilledBins;
    vector<int> SharedFilledBins;
    vector<int> AllFilledBins;

    //Filling with Reco
    for (int i = 0; i < N; i++) {
        double x = Pts.at(0).at(i);
        double y = Pts.at(1).at(i);
        double z = Pts.at(2).at(i);
        double Q = Charges.at(i)/Qtot; //Normalized charge
        hdiff->Fill(x,y,z,-Q); // Charge added in negative in the shared hist "hdiff"
        hq->Fill(x,y,z,Q);
        
        //Adding the bin to the list...
        int bin = hdiff->FindBin(x,y,z);
        bool New = true;
        for (int k = 0, sz = QFilledBins.size(); k < sz; k++){
            int i = QFilledBins.at(k);
            if (i==bin){New = false;}
        }
        //... only if the bin is new
        if (New){QFilledBins.push_back(bin);}
    }
    
    //Filling with Sim
    for (int i = 0; i < Nsim; i++) {
        double x = SimPts.at(0).at(i);
        double y = SimPts.at(1).at(i);
        double z = SimPts.at(2).at(i);
        double E = Energies.at(i) / Etot; //Normalized Energy
        hdiff->Fill(x,y,z,E);
        he->Fill(x,y,z,E);
        
        //Adding the bin to the list if new
        int bin = hdiff->FindBin(x,y,z);
        bool New = true;
        for (int k = 0, sz = EFilledBins.size(); k < sz; k++){
            int i = EFilledBins.at(k);
            if (i==bin){New = false;}
        }
        if (New){EFilledBins.push_back(bin);}
    }
    cout << "Filling process completed" << endl;
    
//     //Possibility to draw 3D Bins but hard to handle:  
//     string c3dbinlabel = "c3dbin_" + plotlabel;
//     TCanvas *c3dbin = new TCanvas(c3dbinlabel.c_str(),"Reconstructed 3D Bins");
//     c3dbin->cd();
//     hq->Draw("BOX2Z");
    
    cout << "Number of bins filled with Simulated Energy Deposition : " << EFilledBins.size()<< endl;
    cout << "Number of bins filled with Reconstructed Energy Deposition : " << QFilledBins.size() << endl;
    
    //Counting bins filled with both Reco and Sim Deposit (SharedFilledBins)
    for (int k = 0, sz = QFilledBins.size(); k < sz; k++){
        int i = QFilledBins.at(k);
        bool NoShared = true;
        for (int l = 0, sze = EFilledBins.size(); l < sze; l++)
        {
            int j = EFilledBins.at(l);
            if (i==j)
            {
                //mutual bin -> added to SharedFilledBins
                SharedFilledBins.push_back(i);
                NoShared = false;
                break;
            }
        }
        if (NoShared)
        {
        //The bins filled with Reco but not Sim are added to AllFilledBins so it will miss simply all bins filled with sim
            AllFilledBins.push_back(i);
        }
    }
    cout << "Number of bins filled with both Simulated and Reconstructed Energy Deposition : " << SharedFilledBins.size() << endl;

    //Adding the missing bins to have all bins filled either with Reco or Sim Deposit
    for (int k = 0, sz = EFilledBins.size(); k < sz; k++){
        int i = EFilledBins.at(k);
        AllFilledBins.push_back(i);
    }

    cout << "Total of bins filled with Simulated/Reconstructed Energy Deposition : " << AllFilledBins.size() << endl;
    
    //Generation Energy Spectrum and linked datas (explained in each title)
    string hnelabel = "hne_" + plotlabel;
    string hnqlabel = "hnq_" + plotlabel;
    string hlabel = "h_" + plotlabel;
    string hrlabel = "hr_" + plotlabel;
    
    TH1D *hne = new TH1D(hnelabel.c_str(), "Energy Spectrum;E/Etot;Nbins",100,0,0.05);
    TH1D *hnq = new TH1D(hnqlabel.c_str(), "Charge Spectrum;Q/Qtot;Nbins",100,0,0.05);
    TH1D *hee = new TH1D(hnelabel.c_str(), "Energy Spectrum (cells weighted with energy);E/Etot;E/Etot",100,0,0.05);
    TH1D *hqq = new TH1D(hnqlabel.c_str(), "Charge Spectrum (cells weighted with charge);Q/Qtot;Q/Qtot",100,0,0.05);
    TH1D *hr = new TH1D(hlabel.c_str(), "Distribution of the energy  by logarithmic energy relative difference in bins;log(|E/Etot - Q/Qtot|/E/Etot);E/Etot",100,-2,1);
    TH1D *h = new TH1D(hrlabel.c_str(), "Distribution of the energy  by logarithmic energy difference in bins;log(|E/Etot - Q/Qtot|);Nbins",100,-6,-1);
    
    double Qnorm(0), Enorm(0),Qshnorm(0), Eshnorm(0);
    
    for (int k = 0, sz = AllFilledBins.size(); k < sz; k++){
        int i = AllFilledBins.at(k);
        double Ebin = he->GetBinContent(i);
        double Qbin = hq->GetBinContent(i);
        double Ediffbin = hdiff->GetBinContent(i);
        double Ereldiffbin = Ediffbin / Ebin;
        hne->Fill((Ebin),1);
        hnq->Fill((Qbin),1);
        hee->Fill((Ebin),Ebin);
        hqq->Fill((Qbin),Qbin);
        
        Qnorm+=Qbin;
        Enorm+=Ebin;
    }
    for (int k = 0, sz = SharedFilledBins.size(); k < sz; k++)
    {
        int i = SharedFilledBins.at(k);
        double Ebin = he->GetBinContent(i);
        double Qbin = hq->GetBinContent(i);
        double Ediffbin = hdiff->GetBinContent(i);
        double Ereldiffbin = Ediffbin / Ebin;
        h->Fill(log10(abs(Ediffbin)),Ebin);
        hr->Fill(log10(abs(Ereldiffbin)),Ebin);
        
        Qshnorm+=Qbin;
        Eshnorm+=Ebin;
    }
            
    
    cout << endl << endl << "***" << endl << "Check Normalization : E.integral = " << Enorm << " ; Q.integral = " << Qnorm << endl << endl;
    cout << endl << endl << "***" << endl << "Shared part : E.integral = " << Eshnorm << " ; Q.integral = " << Qshnorm << endl << endl;
    
    c1->cd();
    c1->Divide(2,2);
    
    c1->cd(1);
    c1->cd(1)->SetLogy();
    hne->SetStats(kFALSE);
    hne->SetLineColor(kBlue);
    hne->Draw("HIST");
    
    hnq->SetStats(kFALSE);
    hnq->SetLineColor(2);
    hnq->SetLineStyle(2);
    hnq->Draw("same HIST");
    
    c1->cd(2);
    auto legend = new TLegend(0,0,1,1);
    legend->SetHeader((plotlabel+" - Energy Spectrum").c_str(),"C");
    legend->AddEntry(hne,"Energy Spectrum (sim data)","L");
    legend->AddEntry(hnq,"Charge Spectrum (reco data)","L");
    legend->Draw();
    
    c1->cd(3);
//     c1->cd(3)->SetLogy();
    hee->SetStats(kFALSE);
    hee->SetLineColor(kBlue);
    hee->Draw("HIST");
    
    hqq->SetStats(kFALSE);
    hqq->SetLineColor(2);
    hqq->SetLineStyle(2);
    hqq->Draw("same HIST");
    
//     c1->cd(3);
//     
//     h->SetStats(kFALSE);
//     h->Draw("HIST");
    
    c1->cd(4);
    hr->SetStats(kFALSE);
    hr->Draw("HIST");
    
    c1->Update();
}




//---------------------------------------------------------------------------------------------





void EnergySpectrum (string path,
                     string label,
                     string ext_reco_file,
                     string plotlabel,
                     double EVmax,
                     double BinSzX, double BinSzY, double BinSzZ, 
                     double xmin=0,double xmax=400, 
                     double ymin = -400 , double ymax = 400, 
                     double zmin = -400, double zmax = 400,
                     int NBinE = 100
                    )
{
    string ext_g4_stage2 = "_g4_stage2.root";
    //Bining
    int const NbinX = (xmax - xmin) / BinSzX;
    int const NbinY = (ymax - ymin) / BinSzY;
    int const NbinZ = (zmax - zmin) / BinSzZ;
    
    //Objects Initialisation
    string clabel = "c_" + plotlabel;
    TCanvas *c1 = new TCanvas(clabel.c_str(),"Simulated - Reconstructed diff Energy");
    
    string hdifflabel = "hdiff_" + plotlabel;
    TH3D *hdiff = new TH3D(hdifflabel.c_str(), "Energy diff;X(cm);Y(cm);Z(cm)", NbinX, xmin, xmax, NbinY, ymin,ymax , NbinZ ,zmin , zmax);
    string helabel = "he_" + plotlabel;
    TH3D *he = new TH3D(helabel.c_str(), "Energy;X(cm);Y(cm);Z(cm)", NbinX, xmin, xmax, NbinY, ymin,ymax , NbinZ ,zmin , zmax);
    string hqlabel = "hq_" + plotlabel;
    TH3D *hq = new TH3D(hqlabel.c_str(), "Charge;X(cm);Y(cm);Z(cm)", NbinX, xmin, xmax, NbinY, ymin,ymax , NbinZ ,zmin , zmax);
    
    //Initialisation of vectors of bin IDs (for comparison of filled bins)
    vector<int> QFilledBins;
    vector<int> EFilledBins;
    vector<int> SharedFilledBins;
    vector<int> AllFilledBins;
    
    for (int evn = 0; evn<EVmax; evn++)
    {
        string evlabel = plotlabel + "_" + to_string(evn+1);
        cout << endl << endl << "-----------------  EVENT #" << evn+1 << "  ------------------------"<< endl;
        std::vector<std::vector<double>> Pts = XYZQ(path, label, ext_reco_file, evn);
        std::vector<std::vector<double>> SimPts = XYZE(path, label, ext_g4_stage2, evn);
        std::vector<double> O = Origin(evlabel,evn);
        
        
        int const N = Pts.at(0).size();
        int const Nsim = SimPts.at(0).size();
        
        //Charge and Energy Integration for the Normalisation
        double Qtot(0);
        for (int i = 0; i < N; i++) 
        {
            Qtot+=Pts.at(3).at(i);
        }
        
        double Etot(0);
            for (int i = 0; i < Nsim; i++) 
        {
            Etot+=SimPts.at(3).at(i);
        }

        //Filling with Reco
        for (int i = 0; i < N; i++) {
            double x = Pts.at(0).at(i);
            double y = Pts.at(1).at(i);
            double z = Pts.at(2).at(i);
            double Q = Pts.at(3).at(i)/Qtot; //Normalized charge
            hdiff->Fill(x,y,z,-Q); // Charge added in negative in the shared hist "hdiff"
            hq->Fill(x,y,z,Q);
            
            //Adding the bin to the list...
            int bin = hdiff->FindBin(x,y,z);
            bool New = true;
            for (int k = 0, sz = QFilledBins.size(); k < sz; k++){
                int i = QFilledBins.at(k);
                if (i==bin){New = false;}
            }
            //... only if the bin is new
            if (New){QFilledBins.push_back(bin);}
        }
        
        //Filling with Sim
        for (int i = 0; i < Nsim; i++) {
            double x = SimPts.at(0).at(i);
            double y = SimPts.at(1).at(i);
            double z = SimPts.at(2).at(i);
            double E = SimPts.at(3).at(i) / Etot; //Normalized Energy
            hdiff->Fill(x,y,z,E);
            he->Fill(x,y,z,E);
            
            //Adding the bin to the list if new
            int bin = hdiff->FindBin(x,y,z);
            bool New = true;
            for (int k = 0, sz = EFilledBins.size(); k < sz; k++){
                int i = EFilledBins.at(k);
                if (i==bin){New = false;}
            }
            if (New){EFilledBins.push_back(bin);}
        }
    }
        
    //Counting bins filled with both Reco and Sim Deposit (SharedFilledBins)
    for (int k = 0, sz = QFilledBins.size(); k < sz; k++){
        int i = QFilledBins.at(k);
        bool NoShared = true;
        for (int l = 0, sze = EFilledBins.size(); l < sze; l++)
        {
            int j = EFilledBins.at(l);
            if (i==j)
            {
                //mutual bin -> added to SharedFilledBins
                SharedFilledBins.push_back(i);
                NoShared = false;
                break;
            }
        }
        if (NoShared)
        {
        //The bins filled with Reco but not Sim are added to AllFilledBins so it will miss simply all bins filled with sim
            AllFilledBins.push_back(i);
        }
    }


    //Adding the missing bins to have all bins filled either with Reco or Sim Deposit
    for (int k = 0, sz = EFilledBins.size(); k < sz; k++){
        int i = EFilledBins.at(k);
        AllFilledBins.push_back(i);
    }

    //Generation Energy Spectrum and linked datas (explained in each title)
    string hnelabel = "hne_" + plotlabel;
    string hnqlabel = "hnq_" + plotlabel;
    string hlabel = "h_" + plotlabel;
    string hrlabel = "hr_" + plotlabel;
    
    TH1D *hne = new TH1D(hnelabel.c_str(), "Energy Spectrum;Ebin/Etot;Nbins",NBinE,0,0.15);
    TH1D *hnq = new TH1D(hnqlabel.c_str(), "Charge Spectrum;Qbin/Qtot;Nbins",NBinE,0,0.15);
    TH1D *hee = new TH1D(hnelabel.c_str(), "Energy Spectrum (cells weighted with energy);E_bin/Etot;E/Etot",NBinE,0,0.15);
    TH1D *hqq = new TH1D(hnqlabel.c_str(), "Charge Spectrum (cells weighted with charge);Q_bin/Qtot;Q/Qtot",NBinE,0,0.15);
    TH1D *hr = new TH1D(hlabel.c_str(), "Distribution of the energy  by logarithmic energy relative difference in bins;log(|E/Etot - Q/Qtot|/E/Etot);E/Etot",NBinE,-2,1);
    TH1D *h = new TH1D(hrlabel.c_str(), "Distribution of the energy  by logarithmic energy difference in bins;log(|E/Etot - Q/Qtot|);Nbins",NBinE,-6,-1);
    
    double Qnorm(0), Enorm(0),Qshnorm(0), Eshnorm(0);
    
    for (int k = 0, sz = AllFilledBins.size(); k < sz; k++){
        int i = AllFilledBins.at(k);
        double Ebin = he->GetBinContent(i);
        double Qbin = hq->GetBinContent(i);
        double Ediffbin = hdiff->GetBinContent(i);
        double Ereldiffbin = Ediffbin / Ebin;
        hne->Fill((Ebin),1);
        hnq->Fill((Qbin),1);
        hee->Fill((Ebin),Ebin);
        hqq->Fill((Qbin),Qbin);
    }
    for (int k = 0, sz = SharedFilledBins.size(); k < sz; k++)
    {
        int i = SharedFilledBins.at(k);
        double Ebin = he->GetBinContent(i);
        double Qbin = hq->GetBinContent(i);
        double Ediffbin = hdiff->GetBinContent(i);
        double Ereldiffbin = Ediffbin / Ebin;
        h->Fill(log10(abs(Ediffbin)),Ebin);
        hr->Fill(log10(abs(Ereldiffbin)),Ebin);
        
        Qshnorm+=Qbin;
        Eshnorm+=Ebin;
    }
            
    c1->cd();
    c1->Divide(2,2);
    
    c1->cd(2);
    c1->cd(2)->SetLogy();
    hne->SetStats(kFALSE);
    hne->SetLineColor(2);
    hne->Draw("HIST");
    
    hnq->SetStats(kFALSE);
    hnq->SetLineColor(kBlue);
    hnq->SetLineStyle(2);
    hnq->Draw("same HIST");
    
    
    
    c1->cd(1);
//     c1->cd(3)->SetLogy();
    hee->SetStats(kFALSE);
    hee->SetLineColor(2);
    hee->Draw("HIST");
    
    hqq->SetStats(kFALSE);
    hqq->SetLineColor(kBlue);
    hqq->SetLineStyle(2);
    hqq->Draw("same HIST");
    
    c1->cd(4);
    hr->SetStats(kFALSE);
    hr->Draw("HIST");
    
    c1->cd(1);
    auto legend = new TLegend();
//     legend->SetHeader((plotlabel+" - Energy Spectrum").c_str(),"C");
    legend->SetHeader("Legend","C");
    legend->AddEntry(hne,"Energy Spectrum (sim data)","L");
    legend->AddEntry(hnq,"Charge Spectrum (reco data)","L");
    legend->Draw();
    
    c1->Update();
}

//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------


void GradEVSr(string path, string label, string ext_reco_file, string plotlabel,
              std::vector<std::vector<double>> Pts,
              std::vector<std::vector<double>> SimPts,
              std::vector<double> O,
              double BinSzX, double BinSzY, double BinSzZ, 
              double xmin=0,double xmax=400, 
              double ymin = -400 , double ymax = 400, 
              double zmin = -400, double zmax = 400,
              int NBinE = 100
              )
{

    int const NbinX = (xmax - xmin) / BinSzX;
    int const NbinY = (ymax - ymin) / BinSzY;
    int const NbinZ = (zmax - zmin) / BinSzZ;

    int const N = Pts.at(0).size();
    int const Nsim = SimPts.at(0).size();

    string clabel="c_"+plotlabel;
    TCanvas *c1 = new TCanvas(clabel.c_str(), "Test");

    string Qgradientlabel = "Qgradient_"+plotlabel;
    string Egradientlabel = "Egradient_"+plotlabel;
    TH1D *Qgrad_hist = new TH1D(Qgradientlabel.c_str(), "test", 100, 0, 400);
    TH1D *Egrad_hist = new TH1D(Egradientlabel.c_str(), "test", 100, 0, 400);

    vector<int>EFilledBins;
    vector<int> QFilledBins;
    vector<int> SharedFilledBins;
    vector<int> AllFilledBins;
    
    string ext_g4_stage2 = "_g4_stage2.root";
    for (int evn = 0; evn<1; evn++)
    {
        string evlabel = plotlabel + "_" + to_string(evn+1);
        cout << endl << endl << "-----------------  EVENT #" << evn+1 << "  ------------------------"<< endl;
        std::vector<std::vector<double>> Pts = XYZQ(path, label, ext_reco_file, evn);
        std::vector<std::vector<double>> SimPts = XYZE(path, label, ext_g4_stage2, evn);

            //Calculate derivative for Q and E
            vector<double> X = Pts.at(0);
            vector<double> Y = Pts.at(1);
            vector<double> Z = Pts.at(2);
            vector<double> Q = Pts.at(3);
            vector<double> Xsim = SimPts.at(0);
            vector<double> Ysim = SimPts.at(1);
            vector<double> Zsim = SimPts.at(2);
            vector<double> E = SimPts.at(3);

            // vector<double> Grad_E = Derivative(Xsim,E) + Derivative(Ysim, E) + Derivative(Zsim, E);
            // vector<double> Grad_Q = Derivative(X, Q) + Derivative(Y, Q) + Derivative(Z, Q);

            //Fill with reco
            for (int i=0; i<N; i++){
                double x = Pts.at(0).at(i);
                double y = Pts.at(1).at(i);
                double z = Pts.at(2).at(i);
                double r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
                double grad_q = Derivative(X, Q).at(i) + Derivative(Y, Q).at(i) + Derivative(Z, Q).at(i);
                Qgrad_hist->Fill(r,grad_q);

                //Adding the bin to the list...
                int bin = Qgrad_hist->FindBin(x);
                bool New = true;
                for (int k = 0, sz = QFilledBins.size(); k < sz; k++){
                    int i = QFilledBins.at(k);
                    if (i==bin){New = false;}
                }
                //... only if the bin is new
                if (New){QFilledBins.push_back(bin);}
            }

            //Fill with sim
            for (int i=0; i<N; i++){
                double xsim = SimPts.at(0).at(i);
                double ysim = SimPts.at(1).at(i);
                double zsim = SimPts.at(2).at(i);
                double rsim = sqrt(pow(xsim,2)+pow(ysim,2)+pow(zsim,2));
                double grad_e = Derivative(Xsim,E).at(i) + Derivative(Ysim, E).at(i) + Derivative(Zsim, E).at(i);
                Egrad_hist->Fill(rsim,grad_e);

                //Adding the bin to the list...
                int bin = Egrad_hist->FindBin(xsim);
                bool New = true;
                for (int k = 0, sz = EFilledBins.size(); k < sz; k++){
                    int i = EFilledBins.at(k);
                    if (i==bin){New = false;}
                }
                //... only if the bin is new
                if (New){EFilledBins.push_back(bin);}
            }
            //Counting bins filled with both Reco and Sim Deposit (SharedFilledBins)
            for (int k = 0, sz = QFilledBins.size(); k < sz; k++){
                int i = QFilledBins.at(k);
                bool NoShared = true;
                for (int l = 0, sze = EFilledBins.size(); l < sze; l++)
                {
                    int j = EFilledBins.at(l);
                    if (i==j)
                    {
                        //mutual bin -> added to SharedFilledBins
                        SharedFilledBins.push_back(i);
                        NoShared = false;
                        break;
                    }
                }
                if (NoShared)
                {
                //The bins filled with Reco but not Sim are added to AllFilledBins so it will miss simply all bins filled with sim
                    AllFilledBins.push_back(i);
                }
            }


            //Adding the missing bins to have all bins filled either with Reco or Sim Deposit
            for (int k = 0, sz = EFilledBins.size(); k < sz; k++){
                int i = EFilledBins.at(k);
                AllFilledBins.push_back(i);
            }

    c1->cd();
    c1->Divide(2,2);

    c1->cd(1);
    Qgrad_hist->SetStats(kFALSE);
    Qgrad_hist->SetLineColor(2);
    Qgrad_hist->Draw("HIST");
    
    c1->cd(2);
    Egrad_hist->SetStats(kFALSE);
    Egrad_hist->SetLineColor(kBlue);
    Egrad_hist->SetLineStyle(2);
    Egrad_hist->Draw("same HIST");
    }

}

vector<double> UnitVectShowerAxe(std::vector<std::vector<double>> Pts,std::vector<double> O)
{
    std::vector<double> B = Barycenter(Pts,Pts.at(3));
    vector<double> v = Vect2Pts(O,B);
    double vNorm = sqrt(scalar(v,v));
    v[0]/=vNorm;
    v[1]/=vNorm;
    v[2]/=vNorm;    
    return v;
}

void Efficiency_criteria(
    string output_file_label,
    std::vector<std::vector<double>> SimPts, //Simulated Points
    std::vector<std::vector<double>> Track, //Pandora Track
    std::vector<double> O, //Origin
    int evn, //number of current event
    std::vector<double> theta_phi //theta and phi angle for simulated track
)
{
    ofstream Flux(output_file_label.c_str(), ios::app);
    double dist = 0;
    double RMS = 0;
    int number_track_point = 0;
    double mean = 0;
    std::vector<double> simu_axis = UnitVectShowerAxe(SimPts, O); // DEJA UNITAIRE
    std::vector<double> dist_data;
    double dist_each_point (0);
    for(int i(0);i<Track.at(0).size();i++)
    {
        std::vector<double> Point_on_Track = {Track.at(0).at(i), Track.at(1).at(i), Track.at(2).at(i)};
        std::vector<double> Vect_on_Track = Vect2Pts(O, Point_on_Track);
        dist_each_point = sqrt(scalar(Vect_on_Track, Vect_on_Track) - pow(scalar(Vect_on_Track, simu_axis),2));
        mean += dist_each_point;
        // sigma += dist_each_point*dist_each_point;
        cout << dist_each_point << endl;
        dist_data.push_back(dist_each_point);

    }
    mean /= Track.at(0).size();

    double sigma = 0;
    for(int i(0);i<dist_data.size();i++)
    {
        sigma += pow(dist_data.at(i)-mean,2);
    }
    sigma/=dist_data.size();

    cout << Track.at(0).size() << endl;
    Flux << mean << ";" << sqrt(sigma) << ";" << theta_phi.at(0) << ";" << theta_phi.at(1) << ";" << evn+1 << endl;
    return;
}

std::vector<double> ThetaPhi(std::vector<std::vector<double>> SimPts)
{

    double Dx = SimPts.at(0).at(SimPts.at(0).size()-1) - SimPts.at(0).at(0);
    double Dy = SimPts.at(1).at(SimPts.at(0).size()-1) - SimPts.at(1).at(0);
    double Dz = SimPts.at(2).at(SimPts.at(0).size()-1) - SimPts.at(2).at(0);

    double norm = sqrt(pow(Dy, 2)+pow(Dz,2));

    if(norm==0) return {0, 0};
    
    double theta = 90 - atan(Dx/norm)*180/3.1415;
    double phi = 0;
    if(Dy>=0)
    {
        phi = atan(Dz/Dy)*180/3.1415;
    }
    else if (Dy<0 && Dz<=0)
    {
        phi = atan(Dz/Dy)*180/3.1415 - 180;
    }
    else if(Dy<0 && Dz>0)
    {
        phi = atan(Dz/Dy)*180/3.1415 + 180;
    }

    std::vector<double> res = {theta, phi};
    return res;
}

void Event_SimVSReco(string path, string label, string ext_reco_file, string save_path, string save_end_label, string ang)
{
    string extension_g4_stage2 = "_g4_stage2.root";
    int n_sigma = 3;
//     EnergySpectrum (label,label,10,2,2,2);
//     EnergySpectrum (label,label+"_1cm_BinShift",10,2,2,2,0+1,400+1,-400+1,400+1,-400+1,400+1);
    int EVmin = 1;
    int EVmax = 1;
    // string output_file_q_vs_ne = "/silver/DUNE/peron/simu/charge_vs_ne/corrected/muon/charge_vs_ne_" + label + "3_sigma_summed_charge_rectified2.txt";
    // ofstream Flux1(output_file_q_vs_ne.c_str());
    // Flux1 << "q ; ne; event" << endl;
    // ofstream Flux1((path+label+"_track_size_"+ang+".txt").c_str());
    // ofstream Flux2((path+label+"_track_lenght_"+ang+".txt").c_str());
    // Flux1 << "Event;Track_size" << endl;
    // Flux2 << "Event;Track_length;Chi2" << endl;
    // string efficiency_label = label+"_efficiency_"+ang;
    // ofstream Flux((path+efficiency_label+".txt").c_str());
    // Flux << "<d>;sigma;theta;phi;event" << endl;
    for (int evn = EVmin-1; evn<EVmax; evn++)
    {
        string evlabel = label + "_" + to_string(evn+1);
        cout << endl << endl << "-----------------  EVENT #" << evn+1 << "  ------------------------"<< endl;
        std::vector<std::vector<double>> Pts = XYZQ(path, label, ext_reco_file, evn);
        std::vector<std::vector<double>> Pts_pandora = XYZQ_pandora(path, label, ext_reco_file, evn);
        std::vector<std::vector<double>> SimPts = XYZE(path, label, extension_g4_stage2, evn);
        // std::vector<std::vector<double>> Track = GetTrack(path, label, ext_reco_file, evn, ang);
        // std::vector<std::vector<double>> OrigSimDepo = GetOriginSimDepo(path, label, ext_reco_file, evn);
        
        cout << Pts.size()<< endl;
        cout << SimPts.size() << endl;
        cout << Pts.at(0).size()<< endl;
        cout << SimPts.at(0).size() << endl;
        // cout << Track.size() << endl;
        // cout << Track.at(0).size() << endl;
        // cout << OrigSimDepo.size() << endl;
        // cout << OrigSimDepo.at(0).size() << endl;

        std::vector<double> O = Origin(label,evn);
        // double ne(0);
        // double drift_velocity = 0.15*1e6; // cm/s
        // double elec_lifetime = 35*1e-3; // s
        // for (int i(0); i<SimPts.at(3).size(); i++)
        // {
        //     double drift_time = SimPts.at(0).at(i)/drift_velocity;
        //     ne+=SimPts.at(3).at(i)*exp(-drift_time/elec_lifetime);
        //     // cout << "Decay ratio : " << exp(-drift_time/elec_lifetime) << endl;
        //     ne+=SimPts.at(3).at(i);
        // }
        // cout << "event_sim_reco ne =" << ne << endl;
        // double qtot(0);
        // for (int i(0); i<Pts.at(3).size(); i++)
        // {
        //     qtot+=Pts.at(3).at(i);
        // }
        // cout << "event_sim_reco qtot = " << qtot << endl;
        // Flux1 << qtot << ";" << ne << ";" << evn+1 << endl;
        
        // Projections_SimReco(Pts,SimPts,evlabel+" without X Correction",O,true);
        // Pts = ShiftXOpti_Pts(Pts,label,evn);
        // std::vector<double> theta_phi = ThetaPhi(SimPts);
        // if(evn<10){Projections_SimReco_compar(save_path, save_end_label, Pts, Pts_pandora, SimPts, Track, evlabel,O);}
        // Efficiency_criteria(path+efficiency_label+".txt", SimPts, Track, O, evn, theta_phi);
        // Projections_SimReco(XYZQ_pandora,OrigSimDepo,evlabel,O,true);
        // Viewer3D(Pts,SimPts,evlabel);
        // Viewer3D(SimPts, OrigSimDepo, evlabel);
//         cout << endl << endl;
        EnergyMap(Pts_pandora,SimPts,evlabel,2,2,2);
        // EnergyRecoData(Pts,SimPts,evlabel,2,2,2,-400,400,-400,400,0,400, label, n_sigma);
// //         cout << endl << endl << "*************" << endl << endl;
// //         EnergyRecoData(Pts,SimPts,evlabel+"_1cmShiftedBins",1,1,1,-400+0.5,400+0.5,-400,400,0,400);
// //         cout <<endl << endl << "***************************************************************" << endl;
    
    // GradEVSr(label, "test", Pts, SimPts, O, 1,1,1);

    }
}
