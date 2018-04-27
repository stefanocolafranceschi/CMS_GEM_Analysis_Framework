//
//  VisualizeUniformity.cpp
//  
//
//  Created by Brian L Dorney on 18/02/16.
//
//

#include "DetectorMPGD.h"
#include "VisualizeUniformity.h"
#include <TAxis.h>
#include <iomanip>      // setprecision
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream, std::stringbuf
#include <TMath.h>
#include "TH3.h"

using std::cout;
using std::endl;
using std::make_shared;
using std::map;
using std::shared_ptr;
using std::string;
using std::tuple;
using std::vector;

using QualityControl::Timing::getString;
using QualityControl::Timing::printROOTFileStatus;

using namespace QualityControl::Uniformity;

//Default Constructor
VisualizeUniformity::VisualizeUniformity(){
    m_bSaveCanvases = false;
    //strCanvIdent = strCanvIdentNoSpec = "Ana";
} //End Default Constructor

//Constructor with Setup & Detector inputs
VisualizeUniformity::VisualizeUniformity(Uniformity::AnalysisSetupUniformity inputSetup, Uniformity::DetectorMPGD inputDet){
    aSetup          = inputSetup;
    m_bSaveCanvases   = false;
    detMPGD         = inputDet;
    //strCanvIdent = strCanvIdentNoSpec = "Ana";
} //End Constructor with Setup & Detector inputs

//Makes a 2D plot of a given observable in the detector's active area
//Takes a std::string which stores the physical filename as input
void VisualizeUniformity::storeCanvasHisto2D(std::string & strOutputROOTFileName, std::string strOption, std::string strObsName, std::string strDrawOption){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    TFile * ptr_fileOutput = new TFile(strOutputROOTFileName.c_str(), strOption.c_str(),"",1);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !ptr_fileOutput->IsOpen() || ptr_fileOutput->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","makeAndStoreCanvasHisto2D","Error: File I/O");
        printROOTFileStatus(ptr_fileOutput);
        printClassMethodMsg("VisualizeUniformity","makeAndStoreCanvasHisto2D", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","makeAndStoreCanvasHisto2D", "\tExiting; Nothing has been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Call the method below
    storeCanvasHisto2D(ptr_fileOutput, strObsName, strDrawOption);
    
    //Close the File
    //------------------------------------------------------
    ptr_fileOutput->Close();
    
    return;
} //End VisualizeUniformity::storeCanvasHisto2D()

//Makes a 2D plot of a given observable in the detector's active area
//Takes a TFile *, which the canvas is writtent to, as input
void VisualizeUniformity::storeCanvasHisto2D(TFile * file_InputRootFile, std::string strObsName, std::string strDrawOption){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    int iNumEta = detMPGD.getNumEtaSectors();
    
    shared_ptr<TH1F> hObs; //Observable to be drawn
    
    ReadoutSectorEta etaSector;
    
    vector<tuple<float,float,float> > vec_tup3DPt;
    
    TGraph2D *g2DObs = new TGraph2D();   //Two dimmensional graph
    
    //Make the Canvas
    //------------------------------------------------------
    TCanvas canv_DetSum( ("canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "2D_AllEta" ).c_str(), ( strObsName + " for All Eta" ).c_str(), 800, 800);
    //canv_DetSum.SetRightMargin(0.5);
    canv_DetSum.SetLeftMargin(0.25);
    canv_DetSum.SetBottomMargin(0.15);
    
    //Set the name of the g2DObs
    //------------------------------------------------------
    g2DObs->SetName( ("g2D_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "_AllEta").c_str() );
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","makeAndStoreCanvasHisto2D","Error: File I/O");
        printROOTFileStatus(file_InputRootFile);
        printClassMethodMsg("VisualizeUniformity","makeAndStoreCanvasHisto2D", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","makeAndStoreCanvasHisto2D", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Get/Make the Summary Directory
    //------------------------------------------------------
    //Check to see if dir_Summary exists already, if not create it
    TDirectory *dir_Summary = file_InputRootFile->GetDirectory("Summary", false, "GetDirectory" );
    if (dir_Summary == nullptr) { //Case: Directory did not exist in file, CREATE
        dir_Summary = file_InputRootFile->mkdir("Summary");
    } //End Case: Directory did not exist in file, CREATE
    
    //Loop Over the detector's Eta Sectors to fill vec_tup3DPt
    //------------------------------------------------------
    for (int iEta=1; iEta <= iNumEta; ++iEta) {
        //Get the histogram & draw it
        etaSector = detMPGD.getEtaSector(iEta);
        hObs = getObsHisto(strObsName, etaSector);
        
        //Fill the points of the 2D graph
        //cout<<"iEta\ti\tL.B.\tU.B.\tPx\tPy\tPz\n";
        float fPx=0, fPy=etaSector.fPos_Y, fObs=0;
        for (int i = ( 1 + (iEta-1) * hObs->GetNbinsX() ); i < ( 1 + iEta * hObs->GetNbinsX() ); ++i) { //Loop Over Points of hObs
            fPx = hObs->GetBinCenter( (iEta * hObs->GetNbinsX() ) - i );
            fObs = hObs->GetBinContent( (iEta * hObs->GetNbinsX() ) - i );
            
            //Need to prevent (0,Y,0)'s from being drawn (causes coincidence points in drawing, bad)
            //First we store all non (0,Y,0) points then in another loop we fill g2DObs
            if( !(fPx == fObs ) ) vec_tup3DPt.push_back( std::make_tuple(fPx, fPy, fObs) );
        } //End Loop Over Points of gObs
    } //End Loop Over Detector's Eta Sector
    
    //Loop Over vec_tup3DPt and set members to g2DObs
    //------------------------------------------------------
    for( int i=0; i < vec_tup3DPt.size(); ++i){ //Loop over vec_tup3DPt
        g2DObs->SetPoint(i,std::get<0>(vec_tup3DPt[i]), std::get<1>(vec_tup3DPt[i]), std::get<2>(vec_tup3DPt[i]));
    } //End Loop over vec_tup3DPt
    
    //Make the color Palette
    //------------------------------------------------------
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    
    //Draw the g2DObs
    //------------------------------------------------------
    canv_DetSum.cd();
    g2DObs->SetTitle("");
    //g2DObs->GetZaxis()->SetTitle( hObs->GetYaxis()->GetTitle() );
    //g2DObs->UseCurrentStyle();
    //g2DObs->SetNpx(200);
    //g2DObs->SetNpy(200);
    g2DObs->Draw( strDrawOption.c_str() );
    canv_DetSum.SetTheta(90);
    canv_DetSum.SetPhi(0.001);
    
    //Setup the TLatex for "CMS Preliminary"
    //------------------------------------------------------
    TLatex latex_CMSPrelim;
    latex_CMSPrelim.SetTextSize(0.05);
    latex_CMSPrelim.DrawLatexNDC(0.1, 0.905, "CMS Preliminary" );
    
    //Write the Canvas to the File
    //------------------------------------------------------
    dir_Summary->cd();
    canv_DetSum.Write();
    if (m_bSaveCanvases) { save2png(canv_DetSum); }
    g2DObs->Write();
    
    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End VisualizeUniformity::makeAndStoreCanvasHisto2D()

void VisualizeUniformity::storeCanvasData(std::string & strOutputROOTFileName, std::string strOption, std::string strObsName, std::string strDrawOption, bool bShiftMean){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    TFile * ptr_fileOutput = new TFile(strOutputROOTFileName.c_str(), strOption.c_str(),"",1);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !ptr_fileOutput->IsOpen() || ptr_fileOutput->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasData","Error: File I/O");
        printROOTFileStatus(ptr_fileOutput);
        printClassMethodMsg("VisualizeUniformity","storeCanvasData", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasData", "\tExiting; Nothing has been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Call the method below
    storeCanvasData(ptr_fileOutput, strObsName, strDrawOption, bShiftMean);
    
    //Close the File
    //------------------------------------------------------
    ptr_fileOutput->Close();
    
    return;
} //End VisualizeUniformity::storeCanvasData()

//Draws the distribution of an observable onto a single pad of canvas
//Takes a TFile * which the histograms are written to as input
void VisualizeUniformity::storeCanvasData(TFile * file_InputRootFile, std::string strObsName, std::string strDrawOption, bool bShiftMean){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    int iNumEta = detMPGD.getNumEtaSectors();
    
    float fMu, fErr_Mu;		//Mean
    float fNormChi2;		//Normalized Chi2
    float fSigma, fErr_Sigma;	//Sigma
    float fPerErr, fErr_PerErr;	//Percent Error
    
    SummaryStatistics statObs;
    
    //Make the Canvas
    //------------------------------------------------------
    TCanvas canv_DetSum( ("canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "Dataset_AllEta" ).c_str(), ( strObsName + " for All Eta" ).c_str(), 800, 800);
    //canv_DetSum.SetRightMargin(0.5);
    canv_DetSum.SetLeftMargin(0.25);
    canv_DetSum.SetBottomMargin(0.15);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph","Error: File I/O");
        printROOTFileStatus(file_InputRootFile);
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Get/Make the Summary Directory
    //------------------------------------------------------
    //Check to see if dir_Summary exists already, if not create it
    TDirectory *dir_Summary = file_InputRootFile->GetDirectory("Summary", false, "GetDirectory" );
    if (dir_Summary == nullptr) { //Case: Directory did not exist in file, CREATE
        dir_Summary = file_InputRootFile->mkdir("Summary");
    } //End Case: Directory did not exist in file, CREATE
    
    //Get the requested dataset
    //------------------------------------------------------
    statObs = getObsData(strObsName);
    
    //Make Shifted Distribution
    //------------------------------------------------------
    int iNbinsX = statObs.hDist->GetNbinsX();
    shared_ptr<TGraphErrors> gDistShifted = std::make_shared<TGraphErrors>( TGraphErrors( iNbinsX+2) );
    for (int i=1; i<=iNbinsX; ++i) {
        //Set values
        gDistShifted->SetPoint( i-1, statObs.hDist->GetBinCenter(i) - statObs.hDist->GetMean()+statObs.hDist->GetMean(), statObs.hDist->GetBinContent(i) );
        gDistShifted->SetPointError( i-1, 0., statObs.hDist->GetBinError(i) );
    } //End Loop Over bins of hSummaryDist
    
    //Add two fake points so we can change the x-axis at will
    gDistShifted->SetPoint(iNbinsX, 0, 0);
    gDistShifted->SetPoint(iNbinsX+1, 3000, 0);
    
    //Make Shifted Fit
    //------------------------------------------------------
    string strShifted = "Shifted";
    TF1 * fitDistShifted = (TF1*) statObs.fitDist->Clone( (statObs.fitDist->GetName() + strShifted).c_str() );
    float fXmax= statObs.fitDist->GetXmax();
    float fXmin= statObs.fitDist->GetXmin();
    fitDistShifted->SetRange(fXmin - statObs.hDist->GetMean(), fXmax - statObs.hDist->GetMean()+statObs.hDist->GetMean() );
    fitDistShifted->SetParameter(1, statObs.fitDist->GetParameter(1) - statObs.hDist->GetMean()+statObs.hDist->GetMean() );
    float shiftedMu;
    //Draw Observable
    TString plottitle="";
    //------------------------------------------------------
    canv_DetSum.cd();
    if (bShiftMean) {
        gDistShifted->Draw( (strDrawOption + "AP").c_str() );
        fitDistShifted->Draw("same");
        
        if (strcmp(strObsName.c_str(), "ResponseFitPkPos")==0 ) {
            plottitle = "Peak pos. of cluster charge fits [ADC counts]";
            gDistShifted->GetXaxis()->SetLabelSize(0.03);
            shiftedMu= statObs.fitDist->GetParameter(1); //Shift removed-statObs.hDist->GetMean();
            //if (bNormalize) {
            //   plottitle = "Normalized Peak pos. of cluster charge fits";
            //   gDistShifted->GetXaxis()->SetLabelSize(0.03);
            //}
        }
        if (strcmp(strObsName.c_str(), "ResponseFitChi2")==0 ) {
            plottitle = "#chi^{2} of cluster charge fits";
            //if (bNormalize) plottitle = "#chi^{2} of cluster charge fits";
        }
        if (strcmp(strObsName.c_str(), "ResponseFitPkRes")==0 ) {
            plottitle = "Rel. E resolution from cluster charge peak fits";
            //if (bNormalize) plottitle = "Normalized Rel. E resolution from cluster charge peak fits";
            gDistShifted->GetXaxis()->SetRangeUser(-2,1);
            canv_DetSum.Update();
        }
        gDistShifted->GetXaxis()->SetTitle(plottitle);
        gDistShifted->GetYaxis()->SetTitle("Number of fits");
        gDistShifted->GetYaxis()->SetTitleOffset(2.5);
        gDistShifted->SetTitle("");
    }
    else{
        statObs.hDist->Draw( strDrawOption.c_str() );
        statObs.fitDist->Draw("same");
        
        if (strcmp(strObsName.c_str(), "ResponseFitPkPos")==0 ) {
            plottitle = "Peak pos. of cluster charge fits [ADC counts]";
            statObs.hDist->GetXaxis()->SetLabelSize(0.03);
            
            //if (bNormalize ) {
            //   plottitle = "Normalized Peak pos. of cluster charge fits";
            //   statObs.hDist->GetXaxis()->SetLabelSize(0.03);
            //}
        }
        if (strcmp(strObsName.c_str(), "ResponseFitChi2")==0 ) {
            plottitle = "#chi^{2} of cluster charge fits";
            //if (bNormalize) plottitle = "#chi^{2} of cluster charge fits";
        }
        if (strcmp(strObsName.c_str(), "ResponseFitPkRes")==0 ) {
            plottitle = "Rel. E resolution from cluster charge peak fits";
            //if (bNormalize) plottitle = "Normalized Rel. E resolution from cluster charge peak fits";
            statObs.hDist->GetXaxis()->SetRangeUser(-2,1);
            canv_DetSum.Update();
        }
        statObs.hDist->GetXaxis()->SetTitle(plottitle);
        statObs.hDist->GetYaxis()->SetTitle("Number of fits");
        statObs.hDist->GetYaxis()->SetTitleOffset(2.5);
        statObs.hDist->SetTitle("");
    }
    
    //Setup the TLatex for "CMS Preliminary"
    //------------------------------------------------------
    fMu		= statObs.fitDist->GetParameter(1);
    if (strcmp(strObsName.c_str(), "ResponseFitPkPos")==0 ) fMu=shiftedMu;
    fNormChi2	= statObs.fitDist->GetChisquare() / statObs.fitDist->GetNDF();
    fSigma	= statObs.fitDist->GetParameter(2);
    fPerErr	= fSigma / fMu;
    int           precision = 1, newprecision = 1;
    fErr_Mu	= statObs.fitDist->GetParError(1);
    fErr_Sigma	= statObs.fitDist->GetParError(2);
    fErr_PerErr	= sqrt( pow(fErr_Sigma / fMu, 2) + pow( ( fSigma * fErr_Mu ) / pow(fMu, 2), 2) );
    
    std::stringstream stream0,stream1,stream2,stream3,stream4,stream5,stream6;
    
    // Clearing streams
    stream0.clear();
    stream1.clear();
    stream2.clear();
    stream3.clear();
    stream4.clear();
    stream5.clear();
    stream6.clear();
    
    TLatex latex_CMSPrelim, latex_Obs_Mean, latex_Obs_NormChi2, latex_Obs_StdDev, latex_Obs_PerErr;
    latex_CMSPrelim.SetTextSize(0.05);
    latex_CMSPrelim.DrawLatexNDC(0.2, 0.905, "CMS Preliminary" );
    
    latex_Obs_NormChi2.SetTextSize(0.03);
    
    // Algorithm to understand number of meaningfull digits for Chi2
    for( int a = 10; a >= 1; a = a - 1 ) {
        if ( abs(fNormChi2 * pow (10, (a+1)) )>1 ) precision=a+1;
    }
    if (precision<5) {
        stream0 << std::fixed << std::setprecision(precision-1) << fNormChi2;
    }
    else {
        newprecision=2;
        stream0 << std::fixed << std::scientific << std::setprecision(newprecision) << fNormChi2;
    }
    latex_Obs_NormChi2.DrawLatexNDC(0.28, 0.85, ( "#chi^{2} / NDF = " + stream0.str() ).c_str() );
    latex_Obs_Mean.SetTextSize(0.03);
    for( int a = 10; a >= 1; a = a - 1 ) {
        if ( abs(fErr_Mu * pow (10, (a+1)) ) >=1 ) precision=a+1;
    }
    
    
    // Algorithm to understand number of meaningfull digits for Mu
    if (precision<5) {
        stream1 << std::fixed << std::setprecision(precision) << fErr_Mu;
    }
    else {
        newprecision=2;
        stream1 << std::fixed << std::scientific << std::setprecision(newprecision) << fErr_Mu;
    }
    
    if (precision<5) {
        stream4 << std::fixed << std::setprecision(precision) << fMu;
    }
    else {
        newprecision=2;
        stream4 << std::fixed << std::scientific << std::setprecision(newprecision) << fMu;
    }
    
    latex_Obs_Mean.DrawLatexNDC(0.28, 0.8, ("#mu = " + stream4.str() + " #pm " + stream1.str() ).c_str() );
    latex_Obs_StdDev.SetTextSize(0.03);
    
    
    // Algorithm to understand number of meaningfull digits for Sigma
    for( int a = 10; a >= 1; a = a - 1 ) {
        if ( abs(fErr_Sigma * pow (10, (a+1)) ) >1 ) precision=a+1;
    }
    if (precision<5) {
        stream5 << std::fixed << std::setprecision(precision) << fSigma;
    }
    else {
        newprecision=2;
        stream5 << std::fixed << std::scientific << std::setprecision(newprecision) << fSigma;
    }
    
    if (precision<5) {
        stream2 << std::fixed << std::setprecision(precision) << fErr_Sigma;
    }
    else {
        newprecision=2;
        stream2 << std::fixed << std::scientific << std::setprecision(newprecision) << fErr_Sigma;
    }
    
    
    latex_Obs_StdDev.DrawLatexNDC(0.28, 0.75, ("#sigma = " + stream5.str() + " #pm " + stream2.str() ).c_str() );
    latex_Obs_PerErr.SetTextSize(0.03);
    
    
    // Algorithm to understand number of meaningfull digits for Per
    for( int a = 10; a >= 1; a = a - 1 ) {
        if ( abs(fErr_PerErr * pow (10, (a+1)) ) >1 ) precision=a+1;
    }
    if (precision<5) {
        stream3 << std::fixed << std::setprecision(precision) << fErr_PerErr;
    }
    else {
        newprecision=2;
        stream3 << std::fixed << std::scientific << std::setprecision(newprecision) << fErr_PerErr;
    }
    
    if (precision<5) {
        stream6 << std::fixed << std::setprecision(precision) << fPerErr;
    }
    else {
        newprecision=2;
        stream6 << std::fixed << std::scientific << std::setprecision(newprecision) << fPerErr;
    }
    
    latex_Obs_PerErr.DrawLatexNDC(0.28, 0.7, ("#frac{#sigma}{#mu} = " + stream6.str() + " #pm " + stream3.str() ).c_str() );
    
    //Write the Canvas to the File
    //------------------------------------------------------
    dir_Summary->cd();
    canv_DetSum.Write();
    if (m_bSaveCanvases) { save2png(canv_DetSum); }
    if (bShiftMean) {
        gDistShifted->Write();
        fitDistShifted->Write();
    }
    else{
        statObs.hDist->Write();
        statObs.fitDist->Write();
    }
    
    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End VisualizeUniformity::storeCanvasData()

void VisualizeUniformity::storeCanvasFits(std::string & strOutputROOTFileName, std::string strOption, std::string strDrawOption){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    TFile * ptr_fileOutput = new TFile(strOutputROOTFileName.c_str(), strOption.c_str(),"",1);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !ptr_fileOutput->IsOpen() || ptr_fileOutput->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasFits","Error: File I/O");
        printROOTFileStatus(ptr_fileOutput);
        printClassMethodMsg("VisualizeUniformity","storeCanvasFits", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasFits", "\tExiting; Nothing has been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Call the method below
    storeCanvasFits(ptr_fileOutput, strDrawOption);
    
    //Close the File
    //------------------------------------------------------
    ptr_fileOutput->Close();
    
    return;
} //End VisualizeUniformity::storeCanvasFits()

void VisualizeUniformity::storeCanvasFits(TFile * file_InputRootFile, std::string strDrawOption){
    //Variable Declaration
    int iNumEta = detMPGD.getNumEtaSectors();
    
    ReadoutSectorEta etaSector = detMPGD.getEtaSector(1);
    
    TH2F *hFitSucess2D = new TH2F("h_Summary_FitSuccess","",3,0.5,etaSector.map_sectorsPhi.size()+0.5, iNumEta, 0.5, iNumEta+0.5 );
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph","Error: File I/O");
        printROOTFileStatus(file_InputRootFile);
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    for (int iEta = 1; iEta <= iNumEta; ++iEta) {
        etaSector = detMPGD.getEtaSector(iEta);
        
        //Get Directory
        //-------------------------------------
        //Check to see if dir_SectorEta directory exists already, if not create it
        TDirectory *dir_SectorEta = file_InputRootFile->GetDirectory( ( "SectorEta" + getString( iEta ) ).c_str(), false, "GetDirectory" );
        if (dir_SectorEta == nullptr) { //Case: Directory did not exist in file, CREATE
            dir_SectorEta = file_InputRootFile->mkdir( ( "SectorEta" + getString( iEta ) ).c_str() );
        } //End Case: Directory did not exist in file, CREATE
        
        for (auto iterPhi = etaSector.map_sectorsPhi.begin(); iterPhi != etaSector.map_sectorsPhi.end(); ++iterPhi) {
            
            //Get Directory
            //-------------------------------------
            //Check to see if dir_SectorPhi directory exists already, if not create it
            TDirectory *dir_SectorPhi = dir_SectorEta->GetDirectory( ( "SectorPhi" + getString( (*iterPhi).first ) ).c_str(), false, "GetDirectory"  );
            if (dir_SectorPhi == nullptr) { //Case: Directory did not exist in file, CREATE
                dir_SectorPhi = dir_SectorEta->mkdir( ( "SectorPhi" + getString( (*iterPhi).first ) ).c_str() );
            } //End Case: Directory did not exist in file, CREATE
            
            hFitSucess2D->SetBinContent( (*iterPhi).first, iEta, (*iterPhi).second.fNFitSuccess / aSetup.iUniformityGranularity );
            
            for (auto iterSlice = (*iterPhi).second.map_slices.begin(); iterSlice != (*iterPhi).second.map_slices.end(); ++iterSlice) {
                //Check to see if dir_Slice exists already, if not create it
                TDirectory *dir_Slice = dir_SectorPhi->GetDirectory( ( "Slice" + getString( (*iterSlice).first ) ).c_str(), false, "GetDirectory"  );
                if (dir_Slice == nullptr) { //Case: Directory did not exist in file, CREATE
                    dir_Slice = dir_SectorPhi->mkdir( ( "Slice" + getString( (*iterSlice).first ) ).c_str() );
                } //End Case: Directory did not exist in file, CREATE
                
                if ( (*iterSlice).second.fitSlice_ClustADC == nullptr ) {
                    continue;
                }
                
                TCanvas *canv_SliceSum = (TCanvas*) (getCanvasSliceFit( (*iterSlice).second, iEta, (*iterPhi).first, (*iterSlice).first, false ) )->Clone();
                TCanvas *canv_SliceDataOverFit = (TCanvas*) (getCanvasSliceFit( (*iterSlice).second, iEta, (*iterPhi).first, (*iterSlice).first, true ) )->Clone();
                
                dir_Slice->cd();
                canv_SliceSum->Write();
                canv_SliceDataOverFit->Write();
            } //End Loop Over iSlice
        } //End Loop Over iPhi
    } //End Loop Over iEta
    
    //Make the Canvas
    //------------------------------------------------------
    TCanvas canv_DetSum( ("canv_" + detMPGD.getNameNoSpecial() + "_FitSuccess" ).c_str(), "Fit Successes by Sector", 800, 800);
    canv_DetSum.SetRightMargin(0.25);
    canv_DetSum.SetLeftMargin(0.25);
    canv_DetSum.SetBottomMargin(0.15);
    
    //Get/Make the Summary Directory
    //------------------------------------------------------
    //Check to see if dir_Summary exists already, if not create it
    TDirectory *dir_Summary = file_InputRootFile->GetDirectory("Summary", false, "GetDirectory" );
    if (dir_Summary == nullptr) { //Case: Directory did not exist in file, CREATE
        dir_Summary = file_InputRootFile->mkdir("Summary");
    } //End Case: Directory did not exist in file, CREATE
    
    //Draw mgraph_Obs
    //------------------------------------------------------
    canv_DetSum.cd();
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("5.2f");
    hFitSucess2D->GetZaxis()->SetRangeUser(0.,1.0);
    hFitSucess2D->GetXaxis()->SetTitle("phi sector");
    hFitSucess2D->GetYaxis()->SetTitle("eta sector");
    hFitSucess2D->SetMarkerSize(2.0);
    hFitSucess2D->Draw( strDrawOption.c_str() );
    TAxis* a = hFitSucess2D->GetXaxis();
    TAxis* b = hFitSucess2D->GetYaxis();
    a->SetNdivisions(3);
    b->SetNdivisions(8);
    //Setup the TLatex for "CMS Preliminary"
    //------------------------------------------------------
    TLatex latex_CMSPrelim;
    latex_CMSPrelim.SetTextSize(0.05);
    latex_CMSPrelim.DrawLatexNDC(0.1, 0.905, "CMS Preliminary" );
    
    //Write the Canvas to the File
    //------------------------------------------------------
    dir_Summary->cd();
    canv_DetSum.Write();
    if (m_bSaveCanvases) { save2png(canv_DetSum); }
    hFitSucess2D->Write();
    
    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End VisualizeUniformity::storeCanvasFits()

//Draws a given observable onto a single pad of canvas
//Takes a std::string which stores the physical filename as input
void VisualizeUniformity::storeCanvasGraph(std::string & strOutputROOTFileName, std::string strOption, std::string strObsName, std::string strDrawOption, bool bShowPhiSegmentation){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    TFile * ptr_fileOutput = new TFile(strOutputROOTFileName.c_str(), strOption.c_str(),"",1);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !ptr_fileOutput->IsOpen() || ptr_fileOutput->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph","Error: File I/O");
        printROOTFileStatus(ptr_fileOutput);
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Call the method below
    storeCanvasGraph(ptr_fileOutput, strObsName, strDrawOption, bShowPhiSegmentation);
    
    //Close the File
    //------------------------------------------------------
    ptr_fileOutput->Close();
    
    return;
} //End VisualizeUniformity::storeCanvasGraph()

//Draws a given observable onto a single pad of canvas
//Takes a TFile * which the histograms are written to as input
void VisualizeUniformity::storeCanvasGraph(TFile * file_InputRootFile, std::string strObsName, std::string strDrawOption, bool bShowPhiSegmentation){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    int iNumEta = detMPGD.getNumEtaSectors();
    
    shared_ptr<TGraphErrors> gObs; //Observable to be drawn
    
    ReadoutSectorEta etaSector;
    
    TLegend *legObs = new TLegend(0.8,0.8,0.95,0.95);//  0.2,0.2,0.6,0.4);
    
    TMultiGraph *mgraph_Obs = new TMultiGraph( ( "mgraph_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "_AllEta" ).c_str(), "");
    
    vector<shared_ptr<TGraphErrors> > vec_gObs;
    
    //Make the Canvas
    //------------------------------------------------------
    TCanvas canv_DetSum( ("canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "_AllEta" ).c_str(), ( strObsName + " for All Eta" ).c_str(), 800, 800);
    //canv_DetSum.SetRightMargin(0.5);
    canv_DetSum.SetLeftMargin(0.25);
    canv_DetSum.SetBottomMargin(0.15);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph","Error: File I/O");
        printROOTFileStatus(file_InputRootFile);
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Get/Make the Summary Directory
    //------------------------------------------------------
    //Check to see if dir_Summary exists already, if not create it
    TDirectory *dir_Summary = file_InputRootFile->GetDirectory("Summary", false, "GetDirectory" );
    if (dir_Summary == nullptr) { //Case: Directory did not exist in file, CREATE
        dir_Summary = file_InputRootFile->mkdir("Summary");
    } //End Case: Directory did not exist in file, CREATE
    
    //Setup the Legend
    //------------------------------------------------------
    legObs->SetNColumns(2);
    legObs->SetFillColor(kWhite);
    legObs->SetLineColor(kBlack);
    
    //Loop Over the detector's Eta Sectors
    //------------------------------------------------------
    for (int iEta=1; iEta <= iNumEta; ++iEta) {
        //Get the histogram & draw it
        etaSector = detMPGD.getEtaSector(iEta);
        gObs = getObsGraph(strObsName, etaSector);
        
        gObs->SetLineColor( Timing::getCyclicColor(iEta) );
        gObs->SetMarkerColor( Timing::getCyclicColor(iEta) );
        
        legObs->AddEntry(gObs.get(), ( "i#eta = " + getString(iEta) ).c_str(), "LPE");
        
        vec_gObs.push_back(gObs);			//Need to keep this pointer alive outside of Loop?
        
        mgraph_Obs->Add( gObs.get() );
    } //End Loop Over Detector's Eta Sector
    
    //Draw mgraph_Obs
    //------------------------------------------------------
    canv_DetSum.cd();
    mgraph_Obs->Draw( strDrawOption.c_str() );
    mgraph_Obs->GetHistogram()->GetXaxis()->SetTitleOffset(1.6);
    mgraph_Obs->GetHistogram()->GetYaxis()->SetTitleOffset(1.6);
    mgraph_Obs->GetHistogram()->GetZaxis()->SetTitleOffset(1.6);
    mgraph_Obs->GetXaxis()->SetTitle("Cluster position X [mm]");
    
    TString plottitle="";
    if (strcmp(strObsName.c_str(), "ResponseFitChi2")==0 ) {
        plottitle = "#chi^{2} of cluster charge fits";
    }
    if (strcmp(strObsName.c_str(), "ResponseFitNormChi2")==0 ) {
        plottitle = "#chi^{2} of cluster charge fits";
    }
    if (strcmp(strObsName.c_str(), "ResponseFitPkPos")==0 ) {
        plottitle = "Peak pos. of cluster charge fits [ADC counts]";
        //if (bNormalize) plottitle = "Normalized Peak pos. of cluster charge fits [ADC counts]";
    }
    if (strcmp(strObsName.c_str(), "ResponseFitPkRes")==0 ) {
        plottitle = "Rel. E resolution from cluster charge peak fits";
        //if (bNormalize) plottitle = "Normalized Rel. E resolution from cluster charge peak fits";
    }
    mgraph_Obs->GetYaxis()->SetTitle(plottitle);
    // int binmax = h->GetMaximumBin();
    //double x = h->GetXaxis()->GetBinCenter(binmax);
    
    float fYmax=mgraph_Obs->GetYaxis()->GetXmax();
    float fYmin= mgraph_Obs->GetYaxis()->GetXmin();
    //cout << " a " << fYmin << " " << fYmax << endl;
    
    if (strcmp(strObsName.c_str(), "ResponseFitPkRes")==0 ) {
        fYmax=1;
        fYmin=0;
    }
    
    mgraph_Obs->GetYaxis()->SetRangeUser(fYmin,fYmax*1.1);
    canv_DetSum.Modified();
    mgraph_Obs->GetYaxis()->SetTitle(plottitle);
    gPad->Modified();
    
    //Setup the TLatex for "CMS Preliminary"
    //------------------------------------------------------
    TLatex latex_CMSPrelim;
    latex_CMSPrelim.SetTextSize(0.05);
    latex_CMSPrelim.DrawLatexNDC(0.1, 0.905, "CMS Preliminary" );
    
    //Setup the iPhi designation
    //------------------------------------------------------
    //etaSector should be set here based on the last iteration of the above loop
    if(bShowPhiSegmentation){ //Case: Show iPhi Segmentation
        for(auto iterPhi = etaSector.map_sectorsPhi.begin(); iterPhi != etaSector.map_sectorsPhi.end(); ++iterPhi){
            //Ensure the canvas is the active canvas (it should be already but who knows...)
            canv_DetSum.cd();
            
            //Declare the TLatex
            TLatex latex_PhiSector;
            
            //Determine the iPhi index
            int iPhiPos = std::distance( etaSector.map_sectorsPhi.begin(), iterPhi);
            
            //Draw the TLatex
            latex_PhiSector.SetTextSize(0.05);
            //latex_PhiSector.DrawLatexNDC(0.125 + 0.875 * ( (iPhiPos) / (float)etaSector.map_sectorsPhi.size() ), 0.8, ( "i#phi = " + getString(iPhiPos+1) ).c_str() );
            
            //Segment the Plot with lines
            if (iPhiPos < (etaSector.map_sectorsPhi.size() - 1) ) { //Case: Not the Last Phi Segment Yet
                TLine line_PhiSeg;
                
                line_PhiSeg.SetLineStyle(2);
                line_PhiSeg.SetLineWidth(2);
                
                //line_PhiSeg.DrawLineNDC( ( (iPhiPos+1) / (float)etaSector.map_sectorsPhi.size() ), 0., ( (iPhiPos+1) / (float)etaSector.map_sectorsPhi.size() ), 1. );
            } //End Case: Not the Last Phi Segment Yet
        } //End Loop Over Sector Phi
    } //End Case: Show iPhi Segmentation
    
    //Draw the Legend
    //------------------------------------------------------
    legObs->Draw("same");
    
    //Write the Canvas to the File
    //------------------------------------------------------
    dir_Summary->cd();
    canv_DetSum.Write();
    if (m_bSaveCanvases) { save2png(canv_DetSum); }
    mgraph_Obs->Write();
    
    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End VisualizeUniformity::storeCanvasGraph()

//Makes a 2D plot of a given observable in the detector's active area
//Takes a std::string which stores the physical filename as input
void VisualizeUniformity::storeCanvasGraph2D(std::string & strOutputROOTFileName, std::string strOption, std::string strObsName, std::string strDrawOption, bool bNormalize){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    TFile * ptr_fileOutput = new TFile(strOutputROOTFileName.c_str(), strOption.c_str(),"",1);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !ptr_fileOutput->IsOpen() || ptr_fileOutput->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph2D","Error: File I/O");
        printROOTFileStatus(ptr_fileOutput);
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph2D", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph2D", "\tExiting; Nothing has been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Call the method below
    storeCanvasGraph2D(ptr_fileOutput, strObsName, strDrawOption, bNormalize);
    
    //Close the File
    //------------------------------------------------------
    ptr_fileOutput->Close();
    
    return;
} //End VisualizeUniformity::storeCanvasGraph2D()

//Makes a 2D plot of a given observable in the detector's active area
//Takes a TFile *, which the canvas is writtent to, as input
void VisualizeUniformity::storeCanvasGraph2D(TFile * file_InputRootFile, std::string strObsName, std::string strDrawOption, bool bNormalize){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    double dAvg = 0.;
    int iNumEta = detMPGD.getNumEtaSectors();
    shared_ptr<TGraphErrors> gObs; //Observable to be drawn
    ReadoutSectorEta etaSector;
    bool logar;
    
    //std::vector<shared_ptr<TGraphErrors> > vec_gObs;
    vector<tuple<double,double,double> > vec_tup3DPt;
    
    // Building the Y scale with variable binninb
    const const Int_t YBINS = 8;
    Double_t yEdges[YBINS + 1] = {1200, 1300, 1400, 1520, 1640, 1780, 1920, 2085, 2245};
    
    TGraph2D *g2DObs = new TGraph2D();   //Two dimmensional TGraph2D
    TH2F *maphisto = new TH2F("maphisto","",aSetup.iUniformityGranularity*3-1,-5.075,5.075,YBINS, yEdges); //8,1200,2200);
    
    float xcorrected, ycorrected, correction, newr, newphi, zvalue=0, oldminimum, oldmaximum;
    
    //Make the Canvas
    //------------------------------------------------------
    TCanvas canv_DetSum( ("canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "2D_AllEta_1" ).c_str(), ( strObsName + " for All Eta" ).c_str(), 800, 800);
    canv_DetSum.SetRightMargin(0.25);
    canv_DetSum.SetLeftMargin(0.25);
    canv_DetSum.SetBottomMargin(0.15);
    
    TCanvas canv_DetSum2( ("canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "2D_AllEta_2lin" ).c_str(), ( strObsName + " for All Eta" ).c_str(), 800, 800);
    canv_DetSum2.SetRightMargin(0.25);
    canv_DetSum2.SetLeftMargin(0.25);
    canv_DetSum2.SetBottomMargin(0.15);
    
    TCanvas canv_DetSum3( ("canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "2D_AllEta_2log" ).c_str(), ( strObsName + " for All Eta" ).c_str(), 800, 800);
    canv_DetSum3.SetRightMargin(0.25);
    canv_DetSum3.SetLeftMargin(0.25);
    canv_DetSum3.SetBottomMargin(0.15);
    
    if (bNormalize) {
        canv_DetSum.SetName( ("canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "Normalized2D_AllEta_1" ).c_str() );
        canv_DetSum2.SetName( ("canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "Normalized2D_AllEta_2lin" ).c_str() );
        canv_DetSum3.SetName( ("canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "Normalized2D_AllEta_2log" ).c_str() );
    }
    
    //2Set the name of the g2DObs
    //------------------------------------------------------
    if (bNormalize) {
        g2DObs->SetName( ("g2D_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "Normalized_AllEta").c_str() );
    }
    else{
        g2DObs->SetName( ("g2D_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "_AllEta").c_str() );
    }
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph2D","Error: File I/O");
        printROOTFileStatus(file_InputRootFile);
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph2D", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasGraph2D", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Get/Make the Summary Directory
    //------------------------------------------------------
    //Check to see if dir_Summary exists already, if not create it
    TDirectory *dir_Summary = file_InputRootFile->GetDirectory("Summary", false, "GetDirectory" );
    if (dir_Summary == nullptr) { //Case: Directory did not exist in file, CREATE
        dir_Summary = file_InputRootFile->mkdir("Summary");
    } //End Case: Directory did not exist in file, CREATE
    
    //Loop Over the detector's Eta Sectors to fill vec_tup3DPt
    //------------------------------------------------------
    for (int iEta=1; iEta <= iNumEta; ++iEta) {
        //Get the histogram & draw it
        etaSector = detMPGD.getEtaSector(iEta);
        gObs = getObsGraph(strObsName, etaSector);
        
        //Fill the points of the 2D graph
        double dPx=0, dPy=etaSector.fPos_Y, dObs=0;
        
        for (int i = ( (iEta-1) * gObs->GetN() ); i < (iEta * gObs->GetN() ); ++i) { //Loop Over Points of gObs
            gObs->GetPoint( (iEta * gObs->GetN() ) - i,dPx, dObs);
            
            //Need to prevent (0,Y,0)'s from being drawn (causes coincidence points in drawing, bad)
            //First we store all non (0,Y,0) points then in another loop we fill g2DObs
            
            if( !(dPx == dObs ) ) vec_tup3DPt.push_back( std::make_tuple(dPx, iEta, dObs) );
        } //End Loop Over Points of gObs
    } //End Loop Over Detector's Eta Sector
    
    //Normalize elements in vec_tup3DPt to the average of the elements of vec_tup3DPt
    //------------------------------------------------------
    //    if ( bNormalize ) { //Case: User wants normalization
    //        dAvg = getObsData(strObsName).fMean;
    //
    //        //Normalize elements of vec_tup3DPt to average
    //        std::transform(vec_tup3DPt.begin(), vec_tup3DPt.end(), vec_tup3DPt.begin(), Uniformity::divides(dAvg) );
    //    } //End Case: User wants normalization
    
    zvalue=0;
    //Loop Over vec_tup3DPt and set members to g2DObs
    //------------------------------------------------------
    for( int i=0; i < vec_tup3DPt.size(); ++i){ //Loop over vec_tup3DPt
        
        xcorrected = std::get<0>(vec_tup3DPt[i]);
        ycorrected =       (detMPGD.getEtaSector(std::get<1>(vec_tup3DPt[i]))).fPos_Y             + 153*TMath::Tan(  TMath::Pi() *84.925/180);
        //cout << " sector values y " << ycorrected << endl;
        zvalue = zvalue + std::get<2>(vec_tup3DPt[i]);
        
        newr = pow( ( pow(xcorrected,2) + pow(ycorrected,2) ), 0.5);
        //newr = std::get<1>(vec_tup3DPt[i]);
        newphi = TMath::ATan( xcorrected / ycorrected    );
        
        g2DObs->SetPoint(i, newphi*180/TMath::Pi(), newr, std::get<2>(vec_tup3DPt[i]));
        
        //for( int kk=0; kk < std::get<2>(vec_tup3DPt[i]); ++kk){
        //   maphisto->Fill(newphi*180/TMath::Pi(), newr);
        //}
        maphisto->Fill(newphi*180/TMath::Pi(), newr, std::get<2>(vec_tup3DPt[i])      );
    } //End Loop over vec_tup3DPt
    
    zvalue = 1/(zvalue / vec_tup3DPt.size());
    cout <<	"Normalization factor is " << zvalue << endl;
    
    if ( bNormalize ) {
        //float factor = 1 / maphisto->Integral(); //1 instead of getentries
        maphisto->Scale(zvalue);
    }
    
    
    //Make the color Palette
    //------------------------------------------------------
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    
    float axis_offset=3;
    
    //Draw the g2DObs
    TString plottitle="";
    //------------------------------------------------------
    if (strcmp(strObsName.c_str(), "ResponseFitPkPos")==0 ) {
        plottitle = ";Cluster position Phi [deg]; Cluster position R [mm];Peak pos. of cluster charge fits [ADC counts]";
        if (bNormalize) plottitle = ";Cluster position Phi [deg]; Cluster position R [mm];Normalized Peak pos. of cluster charge fits";
    }
    if (strcmp(strObsName.c_str(), "ResponseFitPkRes")==0 ) {
        plottitle = ";Cluster position Phi [deg]; Cluster position R [mm];Rel. E resolution from cluster charge peak fits";
        if (bNormalize) plottitle = ";Cluster position Phi [deg]; Cluster position R [mm];Normalized Rel. E resolution from cluster charge peak fits";
        //g2DObs->GetHistogram()->GetZaxis()->SetRangeUser(0,0.8);
        //maphisto->GetZaxis()->SetRangeUser(0,0.8);
    }
    //Setup the TLatex for "CMS Preliminary"
    //------------------------------------------------------
    TLatex latex_CMSPrelim;
    latex_CMSPrelim.SetTextSize(0.05);
    
    dir_Summary->cd();
    
    g2DObs->SetTitle(plottitle);
    g2DObs->GetHistogram()->GetXaxis()->SetTitleOffset(2);
    g2DObs->GetHistogram()->GetYaxis()->SetTitleOffset(2.5);
    g2DObs->GetHistogram()->GetZaxis()->SetTitleOffset(2);
    //g2DObs->UseCurrentStyle();
    //g2DObs->SetNpx(200);
    //g2DObs->SetNpy(200);
    canv_DetSum.cd();
    latex_CMSPrelim.DrawLatexNDC(0.1, 0.915, "CMS Preliminary" );
    g2DObs->Draw( strDrawOption.c_str() );
    canv_DetSum.SetTheta(90);
    canv_DetSum.SetPhi(0.001);
    
    // Preventing triangularized plot from being saved
    if (!(bNormalize)) {
        canv_DetSum.Write();
        if (m_bSaveCanvases) save2png(canv_DetSum);
    }
    
    canv_DetSum2.cd();
    latex_CMSPrelim.DrawLatexNDC(0.1, 0.915, "CMS Preliminary" );
    
    //g2DObs->Draw( "HCOLZ" );
    maphisto->SetTitle(plottitle);
    maphisto->GetXaxis()->SetTitleOffset(2);
    maphisto->GetYaxis()->SetTitleOffset(2.5);
    maphisto->GetZaxis()->SetTitleOffset(2);
    
    
    if (bNormalize) {
        oldminimum = maphisto->GetMinimum();
        oldmaximum = maphisto->GetMaximum();
        maphisto->SetMinimum(0);
        maphisto->SetMaximum(2);
        maphisto->Draw("HCOLZ"); //COL POL");1
    }
    else {
        maphisto->Draw("HCOLZ"); //COL POL");
    }
    
    canv_DetSum2.SetTheta(90);
    canv_DetSum2.SetPhi(0.001);
    canv_DetSum2.Write();
    //maphisto->SetMinimum(oldminimum);
    //maphisto->SetMaximum(oldmaximum);
    
    if (!(bNormalize)) {
        if (strcmp(strObsName.c_str(), "ResponseFitPkRes")==0 ) {
            maphisto->SetMinimum(0);
            maphisto->SetMaximum(0.5);
        }
    }
    if (m_bSaveCanvases) save2png(canv_DetSum2);
    
    canv_DetSum3.cd();
    latex_CMSPrelim.DrawLatexNDC(0.1, 0.915, "CMS Preliminary" );
    gPad->SetLogz();
    
    if (bNormalize) {
        oldminimum = maphisto->GetMinimum();
        oldmaximum = maphisto->GetMaximum();
        maphisto->SetMinimum(0.1);
        maphisto->SetMaximum(10);
        maphisto->Draw("HCOLZ");
    }
    else {
        maphisto->Draw("HCOLZ");
    }
    canv_DetSum3.Write();
    //maphisto->SetMinimum(oldminimum);
    //maphisto->SetMaximum(oldmaximum);
    
    if (m_bSaveCanvases) save2png(canv_DetSum3);
    
    
    //Write the Canvas to the File
    //------------------------------------------------------
    g2DObs->Write();
    
    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End VisualizeUniformity::storeCanvasGraph2D

//Draws a given observable onto a single pad of canvas
//Takes a std::string which stores the physical filename as input
void VisualizeUniformity::storeCanvasHisto(std::string & strOutputROOTFileName, std::string strOption, std::string strObsName, std::string strDrawOption, bool bShowPhiSegmentation){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    TFile * ptr_fileOutput = new TFile(strOutputROOTFileName.c_str(), strOption.c_str(),"",1);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !ptr_fileOutput->IsOpen() || ptr_fileOutput->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasHisto","Error: File I/O");
        printROOTFileStatus(ptr_fileOutput);
        printClassMethodMsg("VisualizeUniformity","storeCanvasHisto", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasHisto", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Call the method below
    //------------------------------------------------------
    storeCanvasHisto(ptr_fileOutput, strObsName, strDrawOption, bShowPhiSegmentation);
    
    //Close the File
    //------------------------------------------------------
    ptr_fileOutput->Close();
    
    return;
} //End VisualizeUniformity::storeCanvasHisto()

//Draws a given observable onto a single pad of canvas
//Takes a TFile * which the histograms are written to as input
void VisualizeUniformity::storeCanvasHisto(TFile * file_InputRootFile, std::string strObsName, std::string strDrawOption, bool bShowPhiSegmentation){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    int iNumEta = detMPGD.getNumEtaSectors();
    
    shared_ptr<TH1F> hObs; //Observable to be drawn
    
    ReadoutSectorEta etaSector;
    
    std::vector<shared_ptr<TH1F> > vec_hObs;
    
    TLegend *legObs = new TLegend(0.8,0.8,0.95,0.95); //0.2 0.2 .6 .4
    
    //Make the Canvas
    //------------------------------------------------------
    TCanvas canv_DetSum( ("canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "_AllEta" ).c_str(), ( strObsName + " for All Eta" ).c_str(), 800, 800);
    //canv_DetSum.SetRightMargin(0.5);
    canv_DetSum.SetLeftMargin(0.25);
    canv_DetSum.SetBottomMargin(0.15);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasHisto","Error: File I/O");
        printROOTFileStatus(file_InputRootFile);
        printClassMethodMsg("VisualizeUniformity","storeCanvasHisto", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasHisto", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Get/Make the Summary Directory
    //------------------------------------------------------
    //Check to see if dir_Summary exists already, if not create it
    TDirectory *dir_Summary = file_InputRootFile->GetDirectory("Summary", false, "GetDirectory" );
    if (dir_Summary == nullptr) { //Case: Directory did not exist in file, CREATE
        dir_Summary = file_InputRootFile->mkdir("Summary");
    } //End Case: Directory did not exist in file, CREATE
    
    //Setup the Legend
    //------------------------------------------------------
    legObs->SetNColumns(2);
    legObs->SetFillColor(kWhite);
    legObs->SetLineColor(kBlack);
    
    //Loop Over the detector's Eta Sectors
    //------------------------------------------------------
    for (int iEta=1; iEta <= iNumEta; ++iEta) {
        //Get the histogram & draw it
        etaSector = detMPGD.getEtaSector(iEta);
        hObs = getObsHisto(strObsName, etaSector);
        
        hObs->SetLineColor( Timing::getCyclicColor(iEta) );
        hObs->SetMarkerColor( Timing::getCyclicColor(iEta) );
        legObs->AddEntry(hObs.get(), ( "i#eta = " + getString(iEta) ).c_str(), "LPE");
        
        vec_hObs.push_back(hObs);			//Need to keep this pointer alive outside of Loop?
        
        canv_DetSum.cd();
        if( 1 == iEta ){
            vec_hObs[iEta-1]->Draw( strDrawOption.c_str() );
        }
        else{
            vec_hObs[iEta-1]->Draw( (strDrawOption + "same").c_str() );
        }
    } //End Loop Over Detector's Eta Sector
    
    //Setup the TLatex for "CMS Preliminary"
    //------------------------------------------------------
    TLatex latex_CMSPrelim;
    latex_CMSPrelim.SetTextSize(0.05);
    latex_CMSPrelim.DrawLatexNDC(0.1, 0.905, "CMS Preliminary" );
    
    //Setup the iPhi designation
    //------------------------------------------------------
    //etaSector should be set here based on the last iteration of the above loop
    if(bShowPhiSegmentation){ //Case: Show iPhi Segmentation
        for(auto iterPhi = etaSector.map_sectorsPhi.begin(); iterPhi != etaSector.map_sectorsPhi.end(); ++iterPhi){
            //Ensure the canvas is the active canvas (it should be already but who knows...)
            canv_DetSum.cd();
            
            //Declare the TLatex
            TLatex latex_PhiSector;
            
            //Determine the iPhi index
            int iPhiPos = std::distance( etaSector.map_sectorsPhi.begin(), iterPhi);
            
            //Draw the TLatex
            latex_PhiSector.SetTextSize(0.05);
            //latex_PhiSector.DrawLatexNDC(0.125 + 0.875 * ( (iPhiPos) / (float)etaSector.map_sectorsPhi.size() ), 0.8, ( "i#phi = " + getString(iPhiPos+1) ).c_str() );
            
            //Segment the Plot with lines
            if (iPhiPos < (etaSector.map_sectorsPhi.size() - 1) ) { //Case: Not the Last Phi Segment Yet
                TLine line_PhiSeg;
                
                line_PhiSeg.SetLineStyle(2);
                line_PhiSeg.SetLineWidth(2);
                
                //line_PhiSeg.DrawLineNDC( ( (iPhiPos+1) / (float)etaSector.map_sectorsPhi.size() ), 0., ( (iPhiPos+1) / (float)etaSector.map_sectorsPhi.size() ), 1. );
            } //End Case: Not the Last Phi Segment Yet
        } //End Loop Over Sector Phi
    } //End Case: Show iPhi Segmentation
    
    //Draw the Legend
    //------------------------------------------------------
    legObs->Draw("same");
    
    //Write the Canvas to the File
    //------------------------------------------------------
    dir_Summary->cd();
    canv_DetSum.Write();
    if (m_bSaveCanvases) { save2png(canv_DetSum); }
    
    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End VisualizeUniformity::storeCanvasHisto()

//This method is longer than I'd like it to be
//But it seems that TCanvas doesn't perpetuate its drawn TObject's
//So passing it to another method by reference keeps the TCanvas alive, but ends up being blank with nothing drawn on it =/
//Draws the distribution pointed to by inputObjPtr on a pad of inputCanvas
//inputCanvas is split into two columns;
//The Pad is created when this method is called; iEta and iNumEta define the pad position automatically
//Odd (even) values of iEta are on the left (right)
//The ReadoutSectorEta is used to determine the location of the ReadoutSectorPhi's
//Takes a std::string which stores the physical filename as input
void VisualizeUniformity::storeCanvasHistoSegmented(std::string & strOutputROOTFileName, std::string strOption, std::string strObsName, std::string strDrawOption, bool bShowPhiSegmentation){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    TFile * ptr_fileOutput = new TFile(strOutputROOTFileName.c_str(), strOption.c_str(),"",1);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !ptr_fileOutput->IsOpen() || ptr_fileOutput->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasHistoSegmented","Error: File I/O");
        printROOTFileStatus(ptr_fileOutput);
        printClassMethodMsg("VisualizeUniformity","storeCanvasHistoSegmented", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasHistoSegmented", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Call the method below
    //------------------------------------------------------
    storeCanvasHistoSegmented(ptr_fileOutput, strObsName, strDrawOption, bShowPhiSegmentation);
    
    //Close the File
    //------------------------------------------------------
    ptr_fileOutput->Close();
    
    return;
} //End VisualizeUniformity::storeCanvasHistoSegmented()

//This method is longer than I'd like it to be
//But it seems that TCanvas doesn't perpetuate its drawn TObject's
//So passing it to another method by reference keeps the TCanvas alive, but ends up being blank with nothing drawn on it =/
//Draws the distribution pointed to by inputObjPtr on a pad of inputCanvas
//inputCanvas is split into two columns;
//The Pad is created when this method is called; iEta and iNumEta define the pad position automatically
//Odd (even) values of iEta are on the left (right)
//The ReadoutSectorEta is used to determine the location of the ReadoutSectorPhi's
//Takes a TFile * which the histograms are written to as input
void VisualizeUniformity::storeCanvasHistoSegmented(TFile * file_InputRootFile, std::string strObsName, std::string strDrawOption, bool bShowPhiSegmentation){
    //Variable Declaration
    int iNumEta = detMPGD.getNumEtaSectors();
    Int_t dxbins, dstartbin, dxendbin;
    TString plottitle;
    float fMaxBinVal = -1;
    
    shared_ptr<TH1F> hObs; //Observable to be drawn
    
    ReadoutSectorEta etaSector;
    
    std::vector<shared_ptr<TH1F> > vec_hObs;
    std::vector<TPad *> vec_padSectorObs;
    
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    
    //Make the Canvas
    //------------------------------------------------------
    TCanvas canv_DetSum( ("canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "_AllEta_Segmented" ).c_str(), ( strObsName + " for All Eta" ).c_str(), 2400, 1000);
    TCanvas canv_DetSum2( ("canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "_AllEta2D" ).c_str(), ( strObsName + " for All Eta" ).c_str(), 2400, 1000);
    //canv_DetSum.SetRightMargin(0.5);
    canv_DetSum.SetLeftMargin(0.25);
    canv_DetSum.SetBottomMargin(0.25);
    canv_DetSum2.SetLeftMargin(0.1);
    canv_DetSum2.SetRightMargin(0.2);
    canv_DetSum2.SetBottomMargin(0.15);
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasHistoSegmented","Error: File I/O");
        printROOTFileStatus(file_InputRootFile);
        printClassMethodMsg("VisualizeUniformity","storeCanvasHistoSegmented", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasHistoSegmented", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Check to see if dir_Summary exists, if not create it
    //------------------------------------------------------
    TDirectory *dir_Summary = file_InputRootFile->GetDirectory("Summary", false, "GetDirectory" );
    if (dir_Summary == nullptr) { //Case: Directory did not exist in file, CREATE
        dir_Summary = file_InputRootFile->mkdir("Summary");
    } //End Case: Directory did not exist in file, CREATE
    
    //Loop Over the detector's Eta sectors to determine the Y-Axis range to be used for all plots
    //------------------------------------------------------
    for (int iEta=1; iEta <= iNumEta; ++iEta) {
        //Get the histogram & draw it
        etaSector = detMPGD.getEtaSector(iEta);
        hObs = getObsHisto(strObsName, etaSector);
        vec_hObs.push_back(hObs);			//Need to keep this pointer alive outside of Loop?
        
        if ( hObs->GetBinContent( hObs->GetMaximumBin() ) > fMaxBinVal ) { //Case: Check for the Maximum Bin Value
            fMaxBinVal = hObs->GetBinContent( hObs->GetMaximumBin() );
        } //End Case: Check for the Maximum Bin Value
    } //End Loop Over Detector's Eta Sectors
    
    //Round fMaxBinVal to the nearest power of ten
    //------------------------------------------------------
    //fMaxBinVal = Uniformity::ceilPowerTen(fMaxBinVal, 0);
    //fMaxBinVal = Uniformity::ceilPowerTen(fMaxBinVal, 1, 0);
    fMaxBinVal = std::ceil(1.1 * (fMaxBinVal / 1e3) ) * 1e3;
    
    
    int histobin;
    
    float bindistance, maxdistance, testbin, oldtestbin;
    for (int iEta=1; iEta <= iNumEta; ++iEta) {
        
        maxdistance = vec_hObs[iEta-1]->GetSize();
        //for (Int_t i=0;i<histobin;i++) {
        //if (i>0) {
        //   oldtestbin = vec_hObs[iEta-1]->GetBinContent(i-1);
        //   testbin = vec_hObs[iEta-1]->GetBinContent(i);
        //if (i==1) {
        //   if (iEta==1) maxdistance = abs(testbin-oldtestbin);
        //}
        //bindistance = abs(testbin-oldtestbin);
        //if (bindistance>maxdistance) maxdistance = bindistance;
        //cout << "i = " << i << " testbin= " << testbin << " oldtestbin = " << oldtestbin << " distance << " << bindistance << " max << " << maxdistance << endl;
        //}
        //}
    }
    
    
    if (strcmp(strObsName.c_str(), "ClustADC")==0 ) {
        dxbins=fMaxBinVal/4000;
        dstartbin=0;
        dxendbin=6000;
        plottitle="Cluster charge [ADC counts]";
    }
    if (strcmp(strObsName.c_str(), "ClustMulti")==0 ) {
        dxbins=7;
        dstartbin=1;
        dxendbin=8;
        plottitle="Cluster Multiplicity";
    }
    if (strcmp(strObsName.c_str(), "ClustPos")==0 ) {
        dstartbin=-200;
        dxendbin=200;
        dxbins=maxdistance;
        plottitle="Cluster Position #left(mm#right)";
    }
    if (strcmp(strObsName.c_str(), "ClustSize")==0 ) {
        dxbins=7;
        dstartbin=1;
        dxendbin=8;
        plottitle="Strip multiplicity of clusters #left(N_{strips}#right)";
    }
    if (strcmp(strObsName.c_str(), "ClustTime")==0 ) {
        dxbins=30;
        dstartbin=1;
        dxendbin=30;
        plottitle="Time Bin [25ns]";
    }
    if (strcmp(strObsName.c_str(), "HitADC")==0 ) {
        dxbins=fMaxBinVal/1000;
        dstartbin=0;
        dxendbin=2000;
        plottitle="Hit charge [ADC counts]";
    }
    if (strcmp(strObsName.c_str(), "HitMulti")==0 ) {
        dxbins=7;
        dstartbin=1;
        dxendbin=8;
        plottitle="Number of strips with hits";
    }
    if (strcmp(strObsName.c_str(), "HitPos")==0 ) {
        dxbins=384;
        dstartbin=1;
        dxendbin=384;
        plottitle="Hit Position #left(<ch>#right)";
    }
    if (strcmp(strObsName.c_str(), "HitSize")==0 ) {
        dxbins=7;
        dstartbin=1;
        dxendbin=8;
        plottitle="Hit multiplicity of clusters #left(N_{strips}#right)";
    }
    if (strcmp(strObsName.c_str(), "HitTime")==0 ) {
        dxbins=30;
        dstartbin=1;
        dxendbin=30;
        plottitle="Time Bin [25ns]";
    }
    
    TH2F *cluster2dmap = new TH2F("cluster2dmap","",dxbins,dstartbin,dxendbin,8,0.5,8.5);
    TAxis* byy = cluster2dmap->GetYaxis();
    byy->SetNdivisions(8);
    
    cluster2dmap->GetXaxis()->SetTitleOffset(1.6);
    gStyle->SetOptStat(0);
    
    //Loop Over the detector's Eta Sectors to make the TCanvas
    //----------------------------------------------------
    
    for (int iEta=1; iEta <= iNumEta; ++iEta) {
        //Get the ReadoutSectorEta
        etaSector = detMPGD.getEtaSector(iEta);
        
        TPad *pad_SectorObs = (TPad *) getPadEta(iEta, iNumEta)->Clone( getNameByIndex(iEta, -1, -1, "pad", "Obs" + getString(iEta) ).c_str() );
        vec_padSectorObs.push_back(pad_SectorObs);	//Need to keep this pointer alive outside of Loop?
        
        canv_DetSum.cd();
        vec_padSectorObs[iEta-1]->Draw();
        vec_padSectorObs[iEta-1]->cd();
        
        //Draw the histogram
        vec_hObs[iEta-1]->GetYaxis()->SetRangeUser(1e-1, fMaxBinVal);
        vec_hObs[iEta-1]->Draw( strDrawOption.c_str() );
        vec_hObs[iEta-1]->GetXaxis()->SetTitleOffset(1.6);
        vec_hObs[iEta-1]->GetYaxis()->SetTitleOffset(2);
        
        histobin = vec_hObs[iEta-1]->GetSize();
        for (Int_t i=0;i<histobin;i++) {
            //cout << " histo contenuto bin " <<  i << " value " << vec_hObs[iEta-1]->GetBinContent(i) << " x=" << vec_hObs[iEta-1]->GetBinCenter(i) <<endl;
            for (Int_t j=0;j<vec_hObs[iEta-1]->GetBinContent(i);j++) {
                cluster2dmap->Fill( vec_hObs[iEta-1]->GetBinCenter(i) ,iEta);
            }
        }
        
        //Setup the TLatex for "CMS Preliminary"
        TLatex latex_CMSPrelim;
        latex_CMSPrelim.SetTextSize(0.05);
        if( 1 == iEta){
            latex_CMSPrelim.DrawLatexNDC(0.1, 0.905, "CMS Preliminary" );
        }
        
        //Setup the TLatex for this iEta sector
        TLatex latex_EtaSector;
        latex_EtaSector.SetTextSize(0.05);
        latex_EtaSector.DrawLatexNDC(0.125, 0.85, ( "i#eta = " + getString(iEta) ).c_str() );
        
        //Setup the iPhi designation
        if(bShowPhiSegmentation){ //Case: Show iPhi Segmentation
            for(auto iterPhi = etaSector.map_sectorsPhi.begin(); iterPhi != etaSector.map_sectorsPhi.end(); ++iterPhi){
                //Ensure the pad is the active pad (it should be already but who knows...)
                vec_padSectorObs[iEta-1]->cd();
                
                //Declare the TLatex
                TLatex latex_PhiSector;
                
                //Determine the iPhi index
                int iPhiPos = std::distance( etaSector.map_sectorsPhi.begin(), iterPhi);
                
                //Draw the TLatex
                latex_PhiSector.SetTextSize(0.05);
                latex_PhiSector.DrawLatexNDC(0.125 + 0.875 * ( (iPhiPos) / (float)etaSector.map_sectorsPhi.size() ), 0.7, ( "i#phi = " + getString(iPhiPos+1) ).c_str() );
                
                //Segment the Plot with lines
                if (iPhiPos < (etaSector.map_sectorsPhi.size() - 1) ) { //Case: Not the Last Phi Segment Yet
                    TLine line_PhiSeg;
                    
                    line_PhiSeg.SetLineStyle(2);
                    line_PhiSeg.SetLineWidth(2);
                    
                    line_PhiSeg.DrawLineNDC( 0.37, 0.1 , 0.37, 0.90 );
                    line_PhiSeg.DrawLineNDC( 0.63, 0.1 , 0.63, 0.90 );
                } //End Case: Not the Last Phi Segment Yet
            } //End Loop Over Sector Phi
        } //End Case: Show iPhi Segmentation
    } //End Loop Over Detector's Eta Sector
    
    //Write the Canvas to the File
    //------------------------------------------------------
    dir_Summary->cd();
    canv_DetSum.Write();
    if (m_bSaveCanvases) { save2png(canv_DetSum); }
    
    canv_DetSum2.cd();
    cluster2dmap->GetXaxis()->SetTitle(plottitle);
    cluster2dmap->GetYaxis()->SetTitle("Eta Sector");
    cluster2dmap->GetZaxis()->SetTitle("N");
    
    cluster2dmap->Draw("hcolz");
    canv_DetSum2.Write();
    if (m_bSaveCanvases) { save2png(canv_DetSum2); }
    
    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End VisualizeUniformity::storeCanvasHistoSegmented()

//Draws the run history of a given observable onto a segment canvas
//Takes a std::string which stores the physical filename as input
void VisualizeUniformity::storeCanvasHisto2DHistorySegmented(std::string & strOutputROOTFileName, std::string strOption, std::string strObsName, std::string strDrawOption, bool bIsEta){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    TFile * ptr_fileOutput = new TFile(strOutputROOTFileName.c_str(), strOption.c_str(),"",1);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !ptr_fileOutput->IsOpen() || ptr_fileOutput->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasHisto2DHistorySegmented","Error: File I/O");
        printROOTFileStatus(ptr_fileOutput);
        printClassMethodMsg("VisualizeUniformity","storeCanvasHisto2DHistorySegmented", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasHisto2DHistorySegmented", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Call the method below
    //------------------------------------------------------
    storeCanvasHisto2DHistorySegmented(ptr_fileOutput, strObsName, strDrawOption, bIsEta);
    
    //Close the File
    //------------------------------------------------------
    ptr_fileOutput->Close();
    
    return;
} //End VisualizeUniformity::storeCanvasHisto2DHistorySegmented()

//Draws the run history of a given observable onto a segment canvas
//Takes a std::string which stores the physical filename as input
void VisualizeUniformity::storeCanvasHisto2DHistorySegmented(TFile * file_InputRootFile, std::string strObsName, std::string strDrawOption, bool bIsEta){
    //Variable Declaration
    int iNumEta = detMPGD.getNumEtaSectors();
    
    float fMaxBinVal = -1;
    
    ReadoutSectorEta etaSector;
    
    map<int, shared_ptr<TH2F> > map_hObs2DRunHistory;
    shared_ptr<TH2F> hObs2D;    //Maybe this should be a map?
    
    map<int, map<int, shared_ptr<TH2F> > > map_hObs2DSummaries; //Outer key -> iEta; inner key-> iPhi
    std::vector<TPad *> vec_padSectorObs;
    
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Make the canvas
    //------------------------------------------------------
    string strCanvName, strCanvName2;
    if (bIsEta) {
        strCanvName = "canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "_AllEta_Segmented";
    }
    else {
        strCanvName = "canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "_AllPhi_Segmented1";
        strCanvName2 = "canv_" + detMPGD.getNameNoSpecial() + "_" + strObsName + "_AllPhi_Segmented2";
    }
    
    TCanvas canv_DetSum( strCanvName.c_str(), strObsName.c_str(), 2400, 1000);
    TCanvas canv_DetSum2( strCanvName2.c_str(), strObsName.c_str(), 2400, 1000);
    //canv_DetSum.SetRightMargin(0.5);
    canv_DetSum.SetLeftMargin(0.15);
    canv_DetSum.SetBottomMargin(0.15);
    canv_DetSum2.SetLeftMargin(0.15);
    canv_DetSum2.SetBottomMargin(0.15);
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeCanvasHisto2DHistorySegmented","Error: File I/O");
        printROOTFileStatus(file_InputRootFile);
        printClassMethodMsg("VisualizeUniformity","storeCanvasHisto2DHistorySegmented", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeCanvasHisto2DHistorySegmented", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Check to see if dir_Summary exists, if not create it
    //------------------------------------------------------
    TDirectory *dir_Summary = file_InputRootFile->GetDirectory("Summary", false, "GetDirectory" );
    if (dir_Summary == nullptr) { //Case: Directory did not exist in file, CREATE
        dir_Summary = file_InputRootFile->mkdir("Summary");
    } //End Case: Directory did not exist in file, CREATE
    
    //Loop Over the detector's readout sectors to get the TH2F to be plotted
    //Also determine maximum bin val
    //------------------------------------------------------
    for (int iEta=1; iEta <= iNumEta; ++iEta) {
        //Get the histogram & draw it
        etaSector = detMPGD.getEtaSector(iEta);
        
        //create temporary container
        map<int, shared_ptr<TH2F> > map_hObs2DTemp;
        
        if (bIsEta) { //Case: run history @ eta level
            map_hObs2DRunHistory = getMapObsHisto2D(strObsName, etaSector);
            
            //TLatex latex_EtaSector;
            //        latex_EtaSector.SetTextSize(0.05);
            //        latex_EtaSector.DrawLatexNDC(0.125, 0.85, ( "i#eta = " + getString(iEta) ).c_str() );
            
            hObs2D = getSummarizedRunHistoryHisto2D(map_hObs2DRunHistory, iEta, -1 );
            
            //Setup the TLatex for this iEta sector
            //        TLatex latex_EtaSector;
            //      latex_EtaSector.SetTextSize(0.05);
            //    latex_EtaSector.DrawLatexNDC(0.125, 0.85, ( "i#eta = " + getString(iEta) ).c_str() );
            
            //Check if there is a new fMaxBinVal
            if (fMaxBinVal < hObs2D->GetBinContent(hObs2D->GetMaximumBin() ) ) {
                fMaxBinVal = hObs2D->GetBinContent(hObs2D->GetMaximumBin() );
            }
            
            map_hObs2DTemp[-1]=hObs2D;
            map_hObs2DSummaries[iEta]=map_hObs2DTemp;
        } //End Case: run history @ eta level
        else { //Case: run history @ phi level
            for (auto iterPhi = etaSector.map_sectorsPhi.begin(); iterPhi != etaSector.map_sectorsPhi.end(); ++iterPhi) { //Loop over phi sectors
                map_hObs2DRunHistory = getMapObsHisto2D(strObsName, (*iterPhi).second );
                
                hObs2D = getSummarizedRunHistoryHisto2D(map_hObs2DRunHistory, iEta, (*iterPhi).first );
                
                //Check if there is a new fMaxBinVal
                if (fMaxBinVal < hObs2D->GetBinContent(hObs2D->GetMaximumBin() ) ) {
                    fMaxBinVal = hObs2D->GetBinContent(hObs2D->GetMaximumBin() );
                }
                
                map_hObs2DTemp[(*iterPhi).first]=hObs2D;
            } //End Loop over phi sectors
            
            map_hObs2DSummaries[iEta]=map_hObs2DTemp;
        } //End Case: run history @ phi level
    } //End Loop Over Detector's Eta Sectors
    
    //Round fMaxBinVal to the nearest power of ten
    //------------------------------------------------------
    //fMaxBinVal = Uniformity::ceilPowerTen(fMaxBinVal, 0);
    //fMaxBinVal = Uniformity::ceilPowerTen(fMaxBinVal, 1, 0);
    fMaxBinVal = std::ceil(1.1 * (fMaxBinVal / 1e3) ) * 1e3;
    int a=-1;
    //Loop Over map_hObs2DSummaries to make the TCanvas
    for (auto iterEta = map_hObs2DSummaries.begin(); iterEta != map_hObs2DSummaries.end(); ++iterEta) { //Loop over eta sectors
        
        for (auto iterPhi = (*iterEta).second.begin(); iterPhi != (*iterEta).second.end(); ++iterPhi) { //Loop over histos stored in phi sectors
            //Declare the TPad
            TPad *pad_SectorObs;
            a++;
            if (bIsEta) { //Case: Arrange Pads in eta grid
                pad_SectorObs = (TPad *) getPadEta( (*iterEta).first, iNumEta)->Clone( getNameByIndex((*iterEta).first, -1, -1, "pad", "Obs" + getString( (*iterEta).first ) ).c_str() );
            } //End Case: Arrange Pads in eta grid
            else { //Case: Arrange Pads in phi grid
                pad_SectorObs = (TPad *) getPadPhi( (*iterEta).first, iNumEta, (*iterPhi).first, (*iterEta).second.size() )->Clone( getNameByIndex((*iterEta).first, (*iterPhi).first, -1, "pad", "Obs" + getString( (*iterEta).first ) + getString( (*iterPhi).first ) ).c_str() );
            } //End Case: Arrange Pads in phi grid
            
            vec_padSectorObs.push_back(pad_SectorObs);	//Need to keep this pointer alive outside of Loop?
            
            canv_DetSum.cd();
            if (a>=12) canv_DetSum2.cd();
            vec_padSectorObs[vec_padSectorObs.size()-1]->Draw();
            vec_padSectorObs[vec_padSectorObs.size()-1]->cd();
            (*iterPhi).second->GetZaxis()->SetRangeUser(1e-1,fMaxBinVal);
            (*iterPhi).second->Draw( strDrawOption.c_str() );
            if (bIsEta) (*iterPhi).second->GetYaxis()->SetTitleOffset(2.5);
            //Setup the TLatex for "CMS Preliminary"
            std::stringstream streameta, streamphi;
            int phi;
            TLatex latex_CMSPrelim, latex_CMSeta;
            latex_CMSPrelim.SetTextSize(0.05);
            //if( 1 == (*iterEta).first ){
            latex_CMSPrelim.DrawLatexNDC(0.1, 0.905, "CMS Preliminary" );
            if (!bIsEta) streameta << int(a/3)+1;
            if (bIsEta) streameta << (a+1);
            if ( (a==1) || (a==4) || (a==7) || (a==10) ) phi=2;
            if ( (a==2) || (a==5) || (a==8) || (a==11) ) phi=3;
            if ( (a==3) || (a==6) || (a==9) || (a==12) ) phi=1;
            
            if ( (a==13) || (a==16) || (a==19) || (a==22) ) phi=2;
            if ( (a==14) || (a==17) || (a==20) || (a==23) ) phi=3;
            if ( (a==15) || (a==18) || (a==21) || (a==24) ) phi=1;
            
            streamphi << phi;
            if (!bIsEta) latex_CMSeta.DrawLatexNDC(0.2, 0.805, ( "Eta = " + streameta.str() + " , Phi = " + streamphi.str() ).c_str() );
            if (bIsEta) latex_CMSeta.DrawLatexNDC(0.2, 0.805, ( "Eta = " + streameta.str()).c_str() );
            //}
            
        } //End Loop over histos stored in phi sectors
    } //End Loop over eta sectors
    
    //Write the Canvas to the File
    //------------------------------------------------------
    dir_Summary->cd();
    canv_DetSum.Write();
    canv_DetSum2.Write();
    if (m_bSaveCanvases) { save2png(canv_DetSum); }
    if (m_bSaveCanvases) { save2png(canv_DetSum2); }
    
    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End VisualizeUniformity::storeCanvasHisto2DHistorySegmented()

//For each member of the input map storeListOfCanvasesGraph is called
//Takes a std::string which stores the physical filename as input
//map_strObsNameAndDrawOpt
//  first   -> Obs Name
//  second  -> Draw option
void VisualizeUniformity::storeListOfCanvasesGraph(std::string & strOutputROOTFileName, std::string strOption, std::map<std::string, std::string> & map_strObsNameAndDrawOpt, bool bShowPhiSegmentation){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    TFile * ptr_fileOutput = new TFile(strOutputROOTFileName.c_str(), strOption.c_str(),"",1);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !ptr_fileOutput->IsOpen() || ptr_fileOutput->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph","Error: File I/O");
        printROOTFileStatus(ptr_fileOutput);
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    storeListOfCanvasesGraph(ptr_fileOutput, map_strObsNameAndDrawOpt, bShowPhiSegmentation);
    
    //Close the File
    //------------------------------------------------------
    ptr_fileOutput->Close();
    
    return;
} //End VisualizeUniformity::storeListOfCanvasesGraph

//For each member of the input map storeListOfCanvasesGraph is called
//Takes a TFile * which the histograms are written to as input
//map_strObsNameAndDrawOpt
//  first   -> Obs Name
//  second  -> Draw option
void VisualizeUniformity::storeListOfCanvasesGraph(TFile * file_InputRootFile, std::map<std::string, std::string> & map_strObsNameAndDrawOpt, bool bShowPhiSegmentation){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph","Error: File I/O");
        printROOTFileStatus(file_InputRootFile);
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    for (auto iterMap = map_strObsNameAndDrawOpt.begin(); iterMap != map_strObsNameAndDrawOpt.end(); ++iterMap) { //Loop over input observables
        
        storeCanvasGraph(file_InputRootFile, (*iterMap).first, (*iterMap).second, bShowPhiSegmentation);
    } //End Loop over input observables
    
    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End VisualizeUniformity::storeListOfCanvasesGraph

//For each member of the input map storeListOfCanvasesHisto is called
//Takes a std::string which stores the physical filename as input
//map_strObsNameAndDrawOpt
//  first   -> Obs Name
//  second  -> Draw option
void VisualizeUniformity::storeListOfCanvasesHisto(std::string & strOutputROOTFileName, std::string strOption, std::map<std::string, std::string> & map_strObsNameAndDrawOpt, bool bShowPhiSegmentation){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    TFile * ptr_fileOutput = new TFile(strOutputROOTFileName.c_str(), strOption.c_str(),"",1);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !ptr_fileOutput->IsOpen() || ptr_fileOutput->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph","Error: File I/O");
        printROOTFileStatus(ptr_fileOutput);
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    storeListOfCanvasesHisto(ptr_fileOutput, map_strObsNameAndDrawOpt, bShowPhiSegmentation);
    
    //Close the File
    //------------------------------------------------------
    ptr_fileOutput->Close();
    
    return;
} //End VisualizeUniformity::storeListOfCanvasesHisto

//For each member of the input map storeListOfCanvasesHisto is called
//Takes a TFile * which the histograms are written to as input
//map_strObsNameAndDrawOpt
//  first   -> Obs Name
//  second  -> Draw option
void VisualizeUniformity::storeListOfCanvasesHisto(TFile * file_InputRootFile, std::map<std::string, std::string> & map_strObsNameAndDrawOpt, bool bShowPhiSegmentation){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph","Error: File I/O");
        printROOTFileStatus(file_InputRootFile);
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    for (auto iterMap = map_strObsNameAndDrawOpt.begin(); iterMap != map_strObsNameAndDrawOpt.end(); ++iterMap) { //Loop over input observables
        
        storeCanvasHisto(file_InputRootFile, (*iterMap).first, (*iterMap).second, bShowPhiSegmentation);
    } //End Loop over input observables
    
    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End VisualizeUniformity::storeListOfCanvasesHisto

//For each member of the input map storeListOfCanvasesHistoSegmented is called
//Takes a std::string which stores the physical filename as input
//map_strObsNameAndDrawOpt
//  first   -> Obs Name
//  second  -> Draw option
void VisualizeUniformity::storeListOfCanvasesHistoSegmented(std::string & strOutputROOTFileName, std::string strOption, std::map<std::string, std::string> & map_strObsNameAndDrawOpt, bool bShowPhiSegmentation){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    TFile * ptr_fileOutput = new TFile(strOutputROOTFileName.c_str(), strOption.c_str(),"",1);
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !ptr_fileOutput->IsOpen() || ptr_fileOutput->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph","Error: File I/O");
        printROOTFileStatus(ptr_fileOutput);
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    storeListOfCanvasesHistoSegmented(ptr_fileOutput, map_strObsNameAndDrawOpt, bShowPhiSegmentation);
    
    //Close the File
    //------------------------------------------------------
    ptr_fileOutput->Close();
    
    return;
} //End VisualizeUniformity::storeListOfCanvasesHistoSegmented

//For each member of the input map storeListOfCanvasesHistoSegmented is called
//Takes a TFile * which the histograms are written to as input
//map_strObsNameAndDrawOpt
//  first   -> Obs Name
//  second  -> Draw option
void VisualizeUniformity::storeListOfCanvasesHistoSegmented(TFile * file_InputRootFile, std::map<std::string, std::string> & map_strObsNameAndDrawOpt, bool bShowPhiSegmentation){
    //TFile does not manage objects
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    
    //Check if File Failed to Open Correctly
    //------------------------------------------------------
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie()  ) {
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph","Error: File I/O");
        printROOTFileStatus(file_InputRootFile);
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("VisualizeUniformity","storeListOfCanvasesGraph", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    for (auto iterMap = map_strObsNameAndDrawOpt.begin(); iterMap != map_strObsNameAndDrawOpt.end(); ++iterMap) { //Loop over input observables
        
        storeCanvasHistoSegmented(file_InputRootFile, (*iterMap).first, (*iterMap).second, bShowPhiSegmentation);
    } //End Loop over input observables
    
    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End VisualizeUniformity::storeListOfCanvasesHistoSegmented

TCanvas * VisualizeUniformity::getCanvasSliceFit(SectorSlice & inputSlice, int iEta, int iPhi, int iSlice, bool bDataOverFit){
    //Variable Declaration
    string strName = ""; (bDataOverFit) ? strName = "clustADCDataOverFit" : strName = "clustADCfit";
    string strFitStatus; (inputSlice.bFitAccepted) ? strFitStatus = "Passed" : strFitStatus = "Failed";
    
    TCanvas * ret_Canvas = new TCanvas( getNameByIndex(iEta, iPhi, iSlice, "canv", strName).c_str(), "Cluster ADC Fit", 800, 800 );
    
    TF1 *func_fitSig, *func_fitBkg;
    TH1F *hDataOverFit;
    
    TLegend *legSlice;
    
    //Draw the histogram and slice
    //------------------------------------------------------
    ret_Canvas->cd();
    if (!bDataOverFit) { //Case: Data & Fit
        //Setup Legend
        legSlice = new TLegend(0.4,0.8,0.95,0.95);
        legSlice->AddEntry(inputSlice.hSlice_ClustADC.get(), "Data", "LPE" );
        legSlice->AddEntry(inputSlice.fitSlice_ClustADC.get(), "Fit", "L" );
        
        //Draw Data & Fit
        inputSlice.hSlice_ClustADC->Draw();
        inputSlice.fitSlice_ClustADC->Draw("same");
        
        //Draw the Signal Only?
        if (aSetup.fitSetup_clustADC.m_strFit_Formula_Sig.length() > 0) {
            func_fitSig = new TF1(
                                  "func_fitSig",
                                  aSetup.fitSetup_clustADC.m_strFit_Formula_Sig.c_str(),
                                  inputSlice.fitSlice_ClustADC->GetXmin(),
                                  inputSlice.fitSlice_ClustADC->GetXmax()
                                  );
            
            //Set the fit parameters from the slice fit
            for (int i=aSetup.fitSetup_clustADC.m_pair_iParamRange_Sig.first;
                 i <= aSetup.fitSetup_clustADC.m_pair_iParamRange_Sig.second;
                 ++i) {
                func_fitSig->SetParameter(i-aSetup.fitSetup_clustADC.m_pair_iParamRange_Sig.first, inputSlice.fitSlice_ClustADC->GetParameter(i) );
            }
            
            func_fitSig->SetLineStyle(2);
            func_fitSig->SetLineWidth(2);
            func_fitSig->SetLineColor(kGreen);
            
            func_fitSig->Draw("same");
            legSlice->AddEntry(func_fitSig, "Fit - Signal Only", "L");
        }
        
        //Draw the Background Only?
        if (aSetup.fitSetup_clustADC.m_strFit_Formula_Bkg.length() > 0) {
            func_fitBkg = new TF1(
                                  "func_fitBkg",
                                  aSetup.fitSetup_clustADC.m_strFit_Formula_Bkg.c_str(),
                                  inputSlice.fitSlice_ClustADC->GetXmin(),
                                  inputSlice.fitSlice_ClustADC->GetXmax()
                                  );
            
            //Set the fit parameters from the slice fit
            for (int i=aSetup.fitSetup_clustADC.m_pair_iParamRange_Bkg.first;
                 i <= aSetup.fitSetup_clustADC.m_pair_iParamRange_Bkg.second;
                 ++i) {
                func_fitBkg->SetParameter(i-aSetup.fitSetup_clustADC.m_pair_iParamRange_Bkg.first, inputSlice.fitSlice_ClustADC->GetParameter(i) );
            }
            
            func_fitBkg->SetLineStyle(2);
            func_fitBkg->SetLineWidth(2);
            func_fitBkg->SetLineColor(kBlue);
            
            func_fitBkg->Draw("same");
            legSlice->AddEntry(func_fitBkg, "Fit - Background Only", "L");
        }
        
        //Draw the legend
        legSlice->Draw("same");
    } //End Case: Data & Fit
    else { //Case: Data / Fit
        //Draw the histogram
        hDataOverFit = (TH1F*) inputSlice.hSlice_ClustADC->Clone( getNameByIndex(iEta, iPhi, iSlice, "canv", strName).c_str() );
        
        hDataOverFit->GetXaxis()->SetRangeUser(inputSlice.fitSlice_ClustADC->GetXmin(), inputSlice.fitSlice_ClustADC->GetXmax() );
        hDataOverFit->GetYaxis()->SetTitle("Data Over Fit #left(A.U.#right)");
        hDataOverFit->GetYaxis()->SetRangeUser(0.5,1.5);
        hDataOverFit->Divide( inputSlice.fitSlice_ClustADC.get() );
        hDataOverFit->Draw();
        
        //Make the line
        TLine line_DataOverFit;
        
        line_DataOverFit.SetLineStyle(2);
        line_DataOverFit.SetLineWidth(2);
        
        line_DataOverFit.DrawLine( inputSlice.fitSlice_ClustADC->GetXmin(), 1., inputSlice.fitSlice_ClustADC->GetXmax(), 1. );
    } //End Case: Data / Fit
    
    //Setup the TLatex for "CMS Preliminary"
    //------------------------------------------------------
    TLatex latex_CMSPrelim;
    latex_CMSPrelim.SetTextSize(0.05);
    latex_CMSPrelim.DrawLatexNDC(0.1, 0.905, "CMS Preliminary" );
    
    //Setup the TLatex for Slice Position
    //------------------------------------------------------
    TLatex latex_SlicePos;
    latex_SlicePos.SetTextSize(0.03);
    latex_SlicePos.DrawLatexNDC(0.55, 0.7, ("Slice Pos #left(mm#right) = " + getString(inputSlice.fPos_Center) ).c_str() );
    
    //Setup the TLatex's for Fit Status
    //------------------------------------------------------
    TLatex latex_FitStatus;
    latex_FitStatus.SetTextSize(0.03);
    latex_FitStatus.DrawLatexNDC(0.55, 0.65, ("Slice Fit = " + strFitStatus ).c_str() );
    
    TLatex latex_MinuitCode;
    latex_MinuitCode.SetTextSize(0.03);
    latex_MinuitCode.DrawLatexNDC(0.55, 0.6, ("Minuit Status Code = " + getString(inputSlice.iMinuitStatus) ).c_str() );
    
    //Setup the TLatex's for Our Internal Checks
    //------------------------------------------------------
    TLatex latex_NormChi2;
    latex_NormChi2.SetTextSize(0.03);
    latex_NormChi2.DrawLatexNDC(0.55, 0.55, ("#chi^{2} / NDF = " + getString(inputSlice.fitSlice_ClustADC->GetChisquare() / inputSlice.fitSlice_ClustADC->GetNDF() ) ).c_str() );
    
    //Get the Peak Position
    auto iterParamPEAK = std::find(aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.begin(), aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.end(), "PEAK");
    int iIdxPk = std::distance(aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.begin(), iterParamPEAK);
    float fPkPos = inputSlice.fitSlice_ClustADC->GetParameter(iIdxPk);
    float fPkPosErr = inputSlice.fitSlice_ClustADC->GetParError(iIdxPk);
    
    TLatex latex_PkPos;
    latex_PkPos.SetTextSize(0.03);
    latex_PkPos.DrawLatexNDC(0.55, 0.50, ("Fitted PkPos = " + getString(fPkPos) + " #pm " + getString(fPkPosErr) ).c_str() );
    
    double dPkLimit_Min, dPkLimit_Max;
    inputSlice.fitSlice_ClustADC->GetParLimits(iIdxPk, dPkLimit_Min, dPkLimit_Max);
    TLatex latex_PkLimits;
    latex_PkLimits.SetTextSize(0.03);
    latex_PkLimits.DrawLatexNDC(0.55,0.45, ("PkPos Limits #epsilon #left(" + getString( dPkLimit_Min ) + "," + getString( dPkLimit_Max ) + "#right)" ).c_str() );
    
    return ret_Canvas;
} //End VisualizeUniformity::getCanvasSliceFit()

SummaryStatistics VisualizeUniformity::getObsData(std::string strObsName){
    //Variable Declaration
    SummaryStatistics ret_stat;
    
    std::transform(strObsName.begin(),strObsName.end(),strObsName.begin(),toupper);
    
    //=======================Fit Result Parameters=======================
    if (0 == strObsName.compare("RESPONSEFITCHI2") ) { //Case: Fit Pk Pos
        ret_stat = detMPGD.getStatNormChi2();
    } //End Case: Fit Pk Pos
    else if (0 == strObsName.compare("RESPONSEFITPKPOS") ) { //Case: Fit Pk Pos
        ret_stat = detMPGD.getStatPkPos();
    } //End Case: Fit Pk Pos
    else if (0 == strObsName.compare("RESPONSEFITPKRES") ) { //Case: Fit Pk Resolution
        ret_stat = detMPGD.getStatPkRes();
    } //End Case: Fit Pk Resolution
    //=======================Unrecognized Parameters=======================
    else{ //Case: Unrecognized Parameter
        cout<<"QualityControl::Uniformity::VisualizeUniformity::getObsHisto() - Parameter " << strObsName.c_str() << " not recognized!!!\n";
    } //End Case: Unrecognized Parameter
    
    return ret_stat;
} //End VisualizeUniformity::getObsGraph()

std::shared_ptr<TGraphErrors> VisualizeUniformity::getObsGraph(std::string strObsName, Uniformity::ReadoutSectorEta &inputEta){
    //Variable Declaration
    std::shared_ptr<TGraphErrors> ret_graph;
    
    std::transform(strObsName.begin(),strObsName.end(),strObsName.begin(),toupper);
    
    //=======================Fit Result Parameters=======================
    if (0 == strObsName.compare("RESPONSEFITCHI2") ) { //Case: Fit Norm Chi2
        ret_graph = inputEta.gEta_ClustADC_Fit_NormChi2;
    } //End Case: Fit Norm Chi2
    else if (0 == strObsName.compare("RESPONSEFITPKPOS") ) { //Case: Fit Pk Pos
        ret_graph = inputEta.gEta_ClustADC_Fit_PkPos;
    } //End Case: Fit Pk Pos
    else if (0 == strObsName.compare("RESPONSEFITPKRES") ) { //Case: Fit Pk Resolution
        ret_graph = inputEta.gEta_ClustADC_Fit_PkRes;
    } //End Case: Fit Pk Resolution
    //=======================Unrecognized Parameters=======================
    else{ //Case: Unrecognized Parameter
        cout<<"QualityControl::Uniformity::VisualizeUniformity::getObsHisto() - Parameter " << strObsName.c_str() << " not recognized!!!\n";
    } //End Case: Unrecognized Parameter
    
    return ret_graph;
} //End VisualizeUniformity::getObsGraph()

std::shared_ptr<TH1F> VisualizeUniformity::getObsHisto(std::string strObsName, Uniformity::ReadoutSector &inputSector){
    //Variable Declaration
    std::shared_ptr<TH1F> ret_histo;
    
    std::transform(strObsName.begin(),strObsName.end(),strObsName.begin(),toupper);
    
    //=======================Cluster Parameters=======================
    if (0 == strObsName.compare("CLUSTADC") ) { //Case: Cluster ADC's
        ret_histo = inputSector.clustHistos.hADC;
    } //End Case: Cluster ADC's
    else if (0 == strObsName.compare("CLUSTMULTI") ) { //Case: Cluster Multi
        ret_histo = inputSector.clustHistos.hMulti;
    } //End Case: Cluster Multi
    else if (0 == strObsName.compare("CLUSTPOS") ) { //Case: Cluster Position
        ret_histo = inputSector.clustHistos.hPos;
    } //End Case: Cluster Position
    else if (0 == strObsName.compare("CLUSTSIZE") ) { //Case: Cluster Size
        ret_histo = inputSector.clustHistos.hSize;
    } //End Case: Cluster Size
    else if (0 == strObsName.compare("CLUSTTIME") ) { //Case: Cluster Time
        ret_histo = inputSector.clustHistos.hTime;
    } //End Case: Cluster Time
    //=======================Fit Parameters=======================
    //=======================Hit Parameters=======================
    else if (0 == strObsName.compare("HITADC") ) { //Case: Hit ADC
        ret_histo = inputSector.hitHistos.hADC;
    } //End Case: Hit ADC
    else if (0 == strObsName.compare("HITMULTI") ) { //Case: Cluster Multi
        ret_histo = inputSector.hitHistos.hMulti;
    } //End Case: Cluster Multi
    else if (0 == strObsName.compare("HITPOS") ) { //Case: Hit Position
        ret_histo = inputSector.hitHistos.hPos;
    } //End Case: Hit Position
    else if (0 == strObsName.compare("HITTIME") ) { //Case: Hit Time
        ret_histo = inputSector.hitHistos.hTime;
    } //End Case: Hit Time
    //=======================Unrecognized Parameters=======================
    else{ //Case: Unrecognized Parameter
        cout<<"QualityControl::Uniformity::VisualizeUniformity::getObsHisto() - Parameter " << strObsName.c_str() << " not recognized!!!\n";
    } //End Case: Unrecognized Parameter
    
    return ret_histo;
} //End VisualizeUniformity::getObsHisto()

std::map<int, std::shared_ptr<TH2F> > VisualizeUniformity::getMapObsHisto2D(std::string strObsName, Uniformity::ReadoutSector &inputSector){
    //Variable Declaration
    map<int, std::shared_ptr<TH2F> > ret_map;
    
    std::transform(strObsName.begin(),strObsName.end(),strObsName.begin(),toupper);
    
    //=======================Cluster Parameters=======================
    if (0 == strObsName.compare("HISTORYCLUSTADC") ) { //Case: Cluster ADC's
        ret_map = inputSector.clustHistos.map_hADC_v_EvtNum_by_Run;
    } //End Case: Cluster ADC's
    else if (0 == strObsName.compare("HISTORYCLUSTTIME") ) { //Case: Cluster Time
        ret_map = inputSector.clustHistos.map_hTime_v_EvtNum_by_Run;
    } //End Case: Cluster Time
    //=======================Fit Parameters=======================
    //=======================Hit Parameters=======================
    else if (0 == strObsName.compare("HISTORYHITADC") ) { //Case: Hit ADC
        ret_map = inputSector.hitHistos.map_hADC_v_EvtNum_by_Run;
    } //End Case: Hit ADC
    else if (0 == strObsName.compare("HISTORYHITTIME") ) { //Case: Hit Time
        ret_map = inputSector.hitHistos.map_hTime_v_EvtNum_by_Run;
    } //End Case: Hit Time
    //=======================Unrecognized Parameters=======================
    else{ //Case: Unrecognized Parameter
        cout<<"QualityControl::Uniformity::VisualizeUniformity::getMapObsHisto2D() - Parameter " << strObsName.c_str() << " not recognized!!!\n";
    } //End Case: Unrecognized Parameter
    
    return ret_map;
} //End VisualizeUniformity::getMapObsHisto2D()

//Summarizes a map of run history histograms
//There are no checks to make sure the map does not return a null pointer
//At this stage if it returns a null pointer there was a problem elsewhere!
shared_ptr<TH2F> VisualizeUniformity::getSummarizedRunHistoryHisto2D(map<int, shared_ptr<TH2F> > inputMapHisto2D, int iEta, int iPhi ){
    //Variable Declaration
    Timing::HistoSetup setupHisto_RunHistory_X, setupHisto_RunHistory_Y;
    
    map<int, shared_ptr<TH2F> >::iterator iterTempHisto2D = inputMapHisto2D.begin();
    
    string strObsName;
    
    shared_ptr<TH2F> ret_histo2D;
    
    vector<string> vec_strParsedName = Timing::getCharSeparatedList( (*iterTempHisto2D).second->GetName(), '_');
    
    //Initialize setupHisto_RunHistory_X
    //------------------------------------------------------
    setupHisto_RunHistory_X.strHisto_Name	= "RunNo";
    setupHisto_RunHistory_X.fHisto_xLower 	= inputMapHisto2D.begin()->first;
    setupHisto_RunHistory_X.fHisto_xUpper 	= inputMapHisto2D.rbegin()->first+1;
    setupHisto_RunHistory_X.iHisto_nBins 	= setupHisto_RunHistory_X.fHisto_xUpper - setupHisto_RunHistory_X.fHisto_xLower;
    setupHisto_RunHistory_X.strHisto_Title_X	= "Run Number";
    gStyle->SetOptStat(0);
    //Initialize setupHisto_RunHistory_Y
    //------------------------------------------------------
    //Check for observable name; will have "clustX" or "hitX" in title
    for (int i=0; i < vec_strParsedName.size(); ++i) {
        if ( vec_strParsedName[i].find("clust") != std::string::npos || vec_strParsedName[i].find("hit") != std::string::npos ) {
            setupHisto_RunHistory_Y.strHisto_Name = vec_strParsedName[i];
            setupHisto_RunHistory_Y.strHisto_Title_X = vec_strParsedName[i];
            
            if ( vec_strParsedName[i]=="clustTime" ) {
                setupHisto_RunHistory_Y.strHisto_Title_X = "Cluster time bin [25 ns]";
            }
            if ( vec_strParsedName[i]=="clustADC" ) {
                setupHisto_RunHistory_Y.strHisto_Title_X = "ADC Cluster time bin [25 ns]";
            }
            
            break;
        }
    } //End Loop over parsed histogram name
    
    setupHisto_RunHistory_Y.iHisto_nBins = (*iterTempHisto2D).second->GetNbinsY();
    setupHisto_RunHistory_Y.fHisto_xLower = (*iterTempHisto2D).second->GetYaxis()->GetBinLowEdge(1);
    setupHisto_RunHistory_Y.fHisto_xUpper = (*iterTempHisto2D).second->GetYaxis()->GetBinLowEdge(setupHisto_RunHistory_Y.iHisto_nBins+1);
    
    //Initialize ret_histo2D
    //------------------------------------------------------
    ret_histo2D = make_shared<TH2F>( getHistogram2D(iEta, iPhi, setupHisto_RunHistory_X, setupHisto_RunHistory_Y ) );
    
    //Loop over map to summarize the stored TH2F's
    //------------------------------------------------------
    for (auto iterHisto2D = inputMapHisto2D.begin(); iterHisto2D != inputMapHisto2D.end(); ++iterHisto2D) {
        //Project the i^th element of the map
        std::shared_ptr<TH1F> hTempHisto = make_shared<TH1F>( *( (TH1F*) (*iterHisto2D).second->ProjectionY( "hTempHisto",1, setupHisto_RunHistory_Y.iHisto_nBins, "" ) ) );
        
        //Fill the slice of ret_histo2D that this corresponds too
        int iBinX = (*iterHisto2D).first - setupHisto_RunHistory_X.fHisto_xLower + 1;
        for (int j=1; j <= hTempHisto->GetNbinsX(); ++j) {
            ret_histo2D->SetBinContent(iBinX, j, hTempHisto->GetBinContent(j) );
            ret_histo2D->SetBinError(iBinX, j, hTempHisto->GetBinError(j) );
        } //End Loop over bins of hTempHisto
    } //End Loop over inputMapHisto2D
    
    ret_histo2D->GetXaxis()->SetNdivisions(setupHisto_RunHistory_X.iHisto_nBins);
    
    //Setup the TLatex for this iEta sector
    //        TLatex latex_EtaSector;
    //      latex_EtaSector.SetTextSize(0.05);
    //    latex_EtaSector.DrawLatexNDC(0.125, 0.85, ( "i#eta = " + getString(iEta) ).c_str() );
    
    return ret_histo2D;
} //End VisualizeUniformity::getSummarizedRunHistoryHisto2D()
