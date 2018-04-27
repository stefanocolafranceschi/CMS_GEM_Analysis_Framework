//
//  AnalyzeResponseUniformity.cpp
//  
//
//  Created by Brian L Dorney on 28/01/16.
//
//

#include "DetectorMPGD.h"
#include "AnalyzeResponseUniformity.h"
#include <TAxis.h>
#include <Math/DistFunc.h>

using std::cout;
using std::endl;
using std::make_shared;
using std::map;
using std::shared_ptr;
using std::string;
using std::vector;

using QualityControl::Plotter::InfoFit;

using QualityControl::Timing::getString;
using QualityControl::Timing::printROOTFileStatus;
using QualityControl::Timing::HistoSetup;
using QualityControl::Timing::stofSafe;

using namespace QualityControl::Uniformity;

//Default Constructor
AnalyzeResponseUniformity::AnalyzeResponseUniformity(){
    strAnalysisName = "analysis";
} //End Default Constructor

//Set inputs at construction
AnalyzeResponseUniformity::AnalyzeResponseUniformity(AnalysisSetupUniformity inputSetup, DetectorMPGD & inputDet){
    strAnalysisName = "analysis";
    
    //Store Analysis Parameters
    aSetup = inputSetup;
    
    //Store Detector
    detMPGD = inputDet;
} //End Constructor

void AnalyzeResponseUniformity::calcStatistics(SummaryStatistics &inputStatObs, std::multiset<float> &mset_fInputObs, string strObsName){
    //Variable Declaration;
    std::vector<float> vec_fInputObsVariance;
    
    std::multiset<float>::iterator iterQ1 = mset_fInputObs.begin();
    std::multiset<float>::iterator iterQ2 = mset_fInputObs.begin();
    std::multiset<float>::iterator iterQ3 = mset_fInputObs.begin();
    
    //Determine max, min, & mean
    inputStatObs.fMax   = *mset_fInputObs.rbegin();   //Last member of the multiset
    inputStatObs.fMin   = *mset_fInputObs.begin();    //First member of the multiset
    inputStatObs.fMean  = std::accumulate( mset_fInputObs.begin(), mset_fInputObs.end(), 0. ) / mset_fInputObs.size();
    
    //Determine standard deviation
    float fAvg = inputStatObs.fMean;
    vec_fInputObsVariance.resize( mset_fInputObs.size() );
    std::transform(mset_fInputObs.begin(), mset_fInputObs.end(), vec_fInputObsVariance.begin(), [fAvg](float x) { return x - fAvg; });
    inputStatObs.fStdDev = std::sqrt( std::inner_product(vec_fInputObsVariance.begin(), vec_fInputObsVariance.end(), vec_fInputObsVariance.begin(), 0.0) / vec_fInputObsVariance.size() );
    
    //Determine Q1, Q2, & Q3
    std::advance(iterQ1, (int)std::ceil( mset_fInputObs.size() * 0.25 ) );  inputStatObs.fQ1 = *iterQ1;
    std::advance(iterQ2, (int)std::ceil( mset_fInputObs.size() * 0.50 ) );  inputStatObs.fQ2 = *iterQ2;
    std::advance(iterQ3, (int)std::ceil( mset_fInputObs.size() * 0.75 ) );  inputStatObs.fQ3 = *iterQ3;
    
    //Determine IQR
    inputStatObs.fIQR = inputStatObs.fQ3 - inputStatObs.fQ1;
    
    //Make distribution
    inputStatObs.hDist = std::make_shared<TH1F>( TH1F( getNameByIndex(-1, -1, -1, "h", strObsName + "Dataset" ).c_str(), "", 40, inputStatObs.fMean - 5. * inputStatObs.fStdDev, inputStatObs.fMean + 5. * inputStatObs.fStdDev) );
    inputStatObs.hDist->Sumw2();
    inputStatObs.hDist->GetXaxis()->SetTitle( strObsName.c_str() );
    inputStatObs.hDist->GetYaxis()->SetTitle( "N" );
    
    //Fill distribution & determine all outliers
    float fLowerBound = inputStatObs.fQ1 - 1.5 * inputStatObs.fIQR;
    float fUpperBound = inputStatObs.fQ3 + 1.5 * inputStatObs.fIQR;
    inRange inDataRange(fLowerBound, fUpperBound);
    
    inputStatObs.mset_fOutliers = mset_fInputObs;
    for (auto iterSet = inputStatObs.mset_fOutliers.begin(); iterSet != inputStatObs.mset_fOutliers.end(); ++iterSet) { //Loop Over input multiset
        inputStatObs.hDist->Fill( (*iterSet) );
    } //End Loop Over input multiset
    
    
    shared_ptr<TF1> fitDist_Gaus = std::make_shared<TF1>( TF1( getNameByIndex(-1, -1, -1, "fit", strObsName + "Dataset" ).c_str(), "gaus(0)", inputStatObs.fMean - 5. * inputStatObs.fStdDev, inputStatObs.fMean + 5. * inputStatObs.fStdDev) );
    fitDist_Gaus->SetParameter(0, inputStatObs.hDist->GetBinContent( inputStatObs.hDist->GetMaximumBin() ) );
    fitDist_Gaus->SetParameter(1, inputStatObs.hDist->GetMean() );
    fitDist_Gaus->SetParameter(2, inputStatObs.hDist->GetRMS() );
    
    inputStatObs.hDist->Fit(fitDist_Gaus.get(),"QM","", inputStatObs.fMean - 5. * inputStatObs.fStdDev, inputStatObs.fMean + 5. * inputStatObs.fStdDev );
    
    shared_ptr<TF1> fitDist_Landau = std::make_shared<TF1>( TF1( getNameByIndex(-1, -1, -1, "fit", strObsName + "Dataset" ).c_str(), "landau(0)", inputStatObs.fMean - 5. * inputStatObs.fStdDev, inputStatObs.fMean + 5. * inputStatObs.fStdDev) );
    fitDist_Landau->SetParameter(0, inputStatObs.hDist->GetBinContent( inputStatObs.hDist->GetMaximumBin() ) );
    fitDist_Landau->SetParameter(1, inputStatObs.hDist->GetMean() );
    fitDist_Landau->SetParameter(2, inputStatObs.hDist->GetRMS() );
    
    inputStatObs.hDist->Fit(fitDist_Landau.get(),"QM","", inputStatObs.fMean - 5. * inputStatObs.fStdDev, inputStatObs.fMean + 5. * inputStatObs.fStdDev );
    
    // stefano chi2 fit
    //shared_ptr<TF1> fitDist_Landau = std::make_shared<TF1>( TF1( getNameByIndex(-1, -1, -1, "fit", strObsName + "Dataset" ).c_str(), "[0]*pow(x,360/2-1)*TMath::Exp(-x/2*)", inputStatObs.fMean - 5. * inputStatObs.fStdDev, inputStatObs.fMean + 5. * inputStatObs.fStdDev) );
    //fitDist_Landau->SetParameter(0, inputStatObs.hDist->GetBinContent( inputStatObs.hDist->GetMaximumBin() ) );
    //TString formula ="[0]*ROOT::Math:::chisquared_pdf(x,"+inputStatObs.hDist->GetEntries()+"))";
    //TF1 * f1 = new TF1("f1","[0]*[1]*[2}*pow(x,360/2-1)*TMath::Exp(-x/2)");
    //TF1 * f1 = new TF1("f1","[0]+[1]*(pow(x,[2]/2-1)*TMath::Exp(-x/2)) /  (pow(2,[2]/2)*TMath::Factorial(int([2]/2)-1))");
    //TF1 * f1 = new TF1("f1","[0]*[1]*[2]*(pow(x,50/2-1)*TMath::Exp(-x/2) ) /  ( (pow(2,50/2) * TMath::Factorial((50/2)-1) ) )",-2,100);
    //TF1 *f1 = new TF1("f1","[0]*[1]*[2]*ROOT::Math::chisquared_pdf(x,3)");
    //f1->SetParLimits(1,inputStatObs.hDist->GetEntries(),inputStatObs.hDist->GetEntries());
    //inputStatObs.hDist->Fit(fitDist_Landau.get(),"QM","", inputStatObs.fMean - 5. * inputStatObs.fStdDev, inputStatObs.fMean + 5.* inputStatObs.fStdDev );
    //cout << "entries " << inputStatObs.hDist->GetEntries() << endl;
    
    
    //shared_ptr<TF1> fitDist_ChisQuare = std::make_shared<TF1>( TF1( getNameByIndex(-1, -1, -1, "fit", strObsName + "Dataset" ).c_str(), "f1", inputStatObs.fMean - 5. * inputStatObs.fStdDev, inputStatObs.fMean + 5. * inputStatObs.fStdDev) );
    //fitDist_ChisQuare->SetParameter(0, inputStatObs.hDist->GetBinContent( inputStatObs.hDist->GetMaximumBin() ) );
    //fitDist_ChisQuare->SetParameter(0,1);
    //fitDist_ChisQuare->SetParameter(1,1);
    //fitDist_ChisQuare->SetParameter(2,100);
    //fitDist_ChisQuare->SetParLimits(0,1,1);
    //fitDist_ChisQuare->SetParLimits(1,1,200);
    //fitDist_ChisQuare->SetParLimits(2,1,500000);
    //inputStatObs.hDist->Fit(fitDist_ChisQuare.get(),"QM","", inputStatObs.fMean - 5. * inputStatObs.fStdDev, inputStatObs.fMean + 5. * inputStatObs.fStdDev );
    
    
    float fNormChi2_Gaus = fitDist_Gaus->GetChisquare() / fitDist_Gaus->GetNDF();
    float fNormChi2_Landau = fitDist_Landau->GetChisquare() / fitDist_Landau->GetNDF();
    
    
    if (fNormChi2_Gaus < fNormChi2_Landau) {
        inputStatObs.fitDist = fitDist_Gaus;
    }
    else{
        inputStatObs.fitDist = fitDist_Landau;
    }
    
    return;
} //End AnalyzeResponseUniformity::calcStatistics()

bool AnalyzeResponseUniformity::isQualityFit(shared_ptr<TF1> fitInput){
    //Variable Declaration
    bool bIsQuality = true;
    
    for( int i=0; i < fitInput->GetNpar(); i++){
        bIsQuality = isQualityFit(fitInput, i);
        
        if(!bIsQuality) break;
    }
    
    return bIsQuality;
} //End AnalyzeResponseUniformity::isQualityFit()

bool AnalyzeResponseUniformity::isQualityFit(shared_ptr<TF1> fitInput, int iPar){
    //Variable Declaration
    double dPar, dPar_Err;
    double dParLimit_Lower, dParLimit_Upper;
    
    dPar = fitInput->GetParameter(iPar);
    dPar_Err = fitInput->GetParError(iPar);
    fitInput->GetParLimits(iPar,dParLimit_Lower,dParLimit_Upper);
    
    if (dPar == 0) {
        return false;
    }
    else if ( ( fabs(dPar - dParLimit_Lower) / dPar ) < 0.001 ){
        return false;
    }
    else if ( ( fabs(dPar - dParLimit_Upper) / dPar) < 0.001 ){
        return false;
    }
    else if ( ( dPar_Err / dPar ) > 0.1 ) {
        return false;
    }
    
    return true;
} //End AnalyzeResponseUniformity::isQualityFit()

TF1 AnalyzeResponseUniformity::getFit(int iEta, int iPhi, int iSlice, InfoFit & setupFit, shared_ptr<TH1F> hInput, TSpectrum &specInput ){
    //Variable Declaration
    float fLimit_Max = hInput->GetBinLowEdge(hInput->GetNbinsX() + 1 );
    float fLimit_Min = hInput->GetBinLowEdge(1);
    
    vector<string>::const_iterator iterVec_IGuess; //Iterator to use for setting initial guess of fit
    vector<float> vec_fFitRange;
    
    for (auto iterRange = setupFit.m_vec_strFit_Range.begin(); iterRange != setupFit.m_vec_strFit_Range.end(); ++iterRange) { //Loop Over Fit Range
        vec_fFitRange.push_back( getParsedInput( (*iterRange), hInput, specInput ) );
    } //End Loop Over Fit Range
    
    if (vec_fFitRange.size() > 1) {
        fLimit_Min = (*std::min_element(vec_fFitRange.begin(), vec_fFitRange.end() ) );
        fLimit_Max = (*std::max_element(vec_fFitRange.begin(), vec_fFitRange.end() ) );
    }
    
    TF1 ret_Func( getNameByIndex(iEta, iPhi, iSlice, "fit", setupFit.m_strFit_Name).c_str(), setupFit.m_strFit_Formula.c_str(), fLimit_Min, fLimit_Max);
    
    //Check to see if the number of parameters in the TF1 meets the expectation
    if ( ret_Func.GetNpar() < setupFit.m_vec_strFit_ParamIGuess.size() || ret_Func.GetNpar() < setupFit.m_vec_strFit_ParamLimit_Min.size() || ret_Func.GetNpar() < setupFit.m_vec_strFit_ParamLimit_Max.size() ) { //Case: Set points for initial parameters do not meet expectations
        
        printClassMethodMsg("AnalyzeResponseUniformity","getFit","Error! Number of Parameters in Function Less Than Requested Initial Guess Parameters!");
        printClassMethodMsg("AnalyzeResponseUniformity","getFit", ("\tNum Parameter: " + getString( ret_Func.GetNpar() ) ).c_str() );
        printClassMethodMsg("AnalyzeResponseUniformity","getFit", ("\tNum Initial Guesses: " + getString( setupFit.m_vec_strFit_ParamIGuess.size() ) ).c_str() );
        printClassMethodMsg("AnalyzeResponseUniformity","getFit", ("\tNum Initial Guess Limits (Min): " + getString( setupFit.m_vec_strFit_ParamLimit_Min.size() ) ).c_str() );
        printClassMethodMsg("AnalyzeResponseUniformity","getFit", ("\tNum Initial Guess Limits (Max): " + getString( setupFit.m_vec_strFit_ParamLimit_Max.size() ) ).c_str() );
        printClassMethodMsg("AnalyzeResponseUniformity","getFit", "No Initial Parameters Have Been Set! Please Cross-Check Input Analysis Config File" );
        
        return ret_Func;
    } //End Case: Set points for initial parameters do not meet expectations
    
    //Set Fit Parameters - Initial Value
    //------------------------------------------------------
    //Keywords are defined in vec_strSupportedKeywords
    for (int i=0; i<setupFit.m_vec_strFit_ParamIGuess.size(); ++i) { //Loop over parameters - Initial Guess
        ret_Func.SetParameter(i, getParsedInput( setupFit.m_vec_strFit_ParamIGuess[i], hInput, specInput) );
    } //End Loop over parameters - Initial Guess
    
    //Set Fit Parameters - Boundaries
    //------------------------------------------------------
    if (setupFit.m_vec_strFit_ParamLimit_Min.size() == setupFit.m_vec_strFit_ParamLimit_Max.size() ) { //Check: Stored Parameter Limits Match
        
        //Here we use vec_strFit_ParamLimit_Min but we know it has the same number of parameters as vec_strFit_ParamLimit_Max
        //For each fit parameter, set the boundary
        for (int i=0; i<setupFit.m_vec_strFit_ParamLimit_Min.size(); ++i) { //Loop over boundary parameters
            fLimit_Min = getParsedInput(setupFit.m_vec_strFit_ParamLimit_Min[i], hInput, specInput);
            fLimit_Max = getParsedInput(setupFit.m_vec_strFit_ParamLimit_Max[i], hInput, specInput);
            
            (fLimit_Max > fLimit_Min) ? ret_Func.SetParLimits(i, fLimit_Min, fLimit_Max ) : ret_Func.SetParLimits(i, fLimit_Max, fLimit_Min );
        } //End Loop over boundary parameters
    } //End Check: Stored Parameter Limits Match
    
    //Set Fit Parameters - Fixed?
    //------------------------------------------------------
    
    //Placeholder; maybe we add functionality in the future
    
    //Set Other Fit Data Members
    //------------------------------------------------------
    ret_Func.SetLineColor(setupFit.m_iColor);
    ret_Func.SetLineWidth(setupFit.m_fSizeLine);
    
    //Return fit
    //------------------------------------------------------
    return ret_Func;
} //End AnalyzeResponseUniformity::getFit()

float AnalyzeResponseUniformity::getParsedInput(std::string &strInputExp, std::shared_ptr<TH1F> hInput, TSpectrum &specInput){
    //Variable Declaration
    map<string, float> map_key2Val;
    
    //Search the input expression for each of the supported keywords
    //Store these Keywords with their values
    for (int i=0; i < vec_strSupportedKeywords.size(); ++i) { //Loop Through Supported Keywords
        if ( strInputExp.find( vec_strSupportedKeywords[i] ) != std::string::npos ) { //Case: Keyword Found!
            map_key2Val[vec_strSupportedKeywords[i]] = getValByKeyword( vec_strSupportedKeywords[i], hInput, specInput );
        } //End Case: Keyword Found!
    } //End Loop Through Supported Keywords
    
    //Check if map_key2Val has any entries, if so user requested complex expression; parse!
    //If map_key2Val is empty, user has a numeric input; convert to float!
    if (map_key2Val.size() > 0) { //Case: Complex Expression!
        //Setup the expression parser
        symbol_table_t symbol_table;    //Stores the variables in expression and maps them to C++ objects
        expression_t expression;        //Stores the actual expression & the symbol table
        parser_t parser;                //Parses the information for evaluation
        
        //Load all found keywords into the symbol table
        for (auto iterMap = map_key2Val.begin(); iterMap != map_key2Val.end(); ++iterMap) {
            symbol_table.add_variable( (*iterMap).first, (*iterMap).second);
        }
        
        //Give the expression the variables it should have
        expression.register_symbol_table(symbol_table);
        
        //Compile the parsing
        parser.compile(strInputExp, expression);
        
        //Return value to the user
        return expression.value();
    } //End Case: Complex Expression!
    else{ //Case: Numeric Input
        return stofSafe( strInputExp );
    } //End Case: Numeric Input
} //End AnalyzeResponseUniformity::getParsedInput()

TGraphErrors AnalyzeResponseUniformity::getGraph(int iEta, int iPhi, HistoSetup & setupHisto){
    //Variable Declaration
    string strPrefix = "g";
    string strName = getNameByIndex(iEta, iPhi, -1, strPrefix, setupHisto.strHisto_Name);
    
    //Graph Declaration
    TGraphErrors ret_Graph;
    
    //Set Histo Data Members
    ret_Graph.GetXaxis()->SetTitle( setupHisto.strHisto_Title_X.c_str() );
    
    return ret_Graph;
} //End AnalyzeResponseUniformity::getGraph()

//Returns a histogram whose parmeters match those defined in hte input HistoSetup object
TH1F AnalyzeResponseUniformity::getHistogram(int iEta, int iPhi, HistoSetup &setupHisto){
    //Variable Declaration
    string strPrefix = "h";
    string strName = getNameByIndex(iEta, iPhi, -1, strPrefix, setupHisto.strHisto_Name);
    
    //Histo Declaration
    TH1F ret_Histo(strName.c_str(), "", setupHisto.iHisto_nBins, setupHisto.fHisto_xLower, setupHisto.fHisto_xUpper );
    
    //Set Histo Data Members
    ret_Histo.SetXTitle( setupHisto.strHisto_Title_X.c_str() );
    ret_Histo.SetYTitle( setupHisto.strHisto_Title_Y.c_str() );
    
    ret_Histo.Sumw2();
    
    //Set Directory to the global directory
    ret_Histo.SetDirectory(gROOT);
    
    //Return Histogram
    return ret_Histo;
} //End AnalyzeResponseUniformity::getHistogram()

TH2F AnalyzeResponseUniformity::getHistogram2D(int iEta, int iPhi, HistoSetup &setupHisto_X, HistoSetup &setupHisto_Y){
    //Variable Declaration
    string strPrefix = "h";
    string strName = getNameByIndex(iEta, iPhi, -1, strPrefix, setupHisto_Y.strHisto_Name + "_v_" + setupHisto_X.strHisto_Name );
    
    //Histo Declaration
    TH2F ret_Histo(strName.c_str(), "", setupHisto_X.iHisto_nBins, setupHisto_X.fHisto_xLower, setupHisto_X.fHisto_xUpper, setupHisto_Y.iHisto_nBins, setupHisto_Y.fHisto_xLower, setupHisto_Y.fHisto_xUpper);
    
    //Set Histo Data Members
    //Setup the TLatex for this iEta sector
    TLatex latex_EtaSector;
    latex_EtaSector.SetTextSize(0.05);
    latex_EtaSector.DrawLatexNDC(0.125, 0.85, ( "i#eta = " + getString(iEta) ).c_str() );
    
    ret_Histo.SetXTitle( setupHisto_X.strHisto_Title_X.c_str() );
    ret_Histo.SetYTitle( setupHisto_Y.strHisto_Title_X.c_str() );
    
    ret_Histo.Sumw2();
    
    //Set Directory to the global directory
    ret_Histo.SetDirectory(gROOT);
    
    //Return Histogram
    return ret_Histo;
} //End AnalyzeResponseUniformity::getHistogram2D

//Formats a given input string such that it follows the iEta, iPhi, iSlice naming convention
string AnalyzeResponseUniformity::getNameByIndex(int iEta, int iPhi, int iSlice, std::string & strInputPrefix, std::string & strInputName){
    //Variable Declaration
    string ret_Name = "";
    
    if (iSlice > -1) {
        ret_Name = strInputPrefix + "_iEta" + getString(iEta) + "iPhi" + getString(iPhi) + "Slice" + getString(iSlice) + "_" + strInputName;
    }
    else if (iPhi > -1){ //Case: Specific (iEta,iPhi) sector
        ret_Name = strInputPrefix + "_iEta" + getString(iEta) + "iPhi" + getString(iPhi) + "_" + strInputName;
    } //End Case: Specific (iEta,iPhi) sector
    else if (iEta > -1){
        ret_Name = strInputPrefix + "_iEta" + getString(iEta) + "_" + strInputName;
    }
    else{ //Case: iEta Sector, sum over sector's iPhi
        ret_Name = strInputPrefix + "_Summary_" + strInputName;
    } //End Case: iEta Sector, sum over sector's iPhi
    
    return ret_Name;
} //End AnalyzeResponseUniformity::getNameByIndex()

//Formats a given input string such that it follows the iEta, iPhi, iSlice naming convention
//Will take text written in quotes without requiring a variable
string AnalyzeResponseUniformity::getNameByIndex(int iEta, int iPhi, int iSlice, const std::string & strInputPrefix, const std::string & strInputName){
    //Variable Declaration
    string ret_Name = "";
    
    if (iSlice > -1) {
        ret_Name = strInputPrefix + "_iEta" + getString(iEta) + "iPhi" + getString(iPhi) + "Slice" + getString(iSlice) + "_" + strInputName;
    }
    else if (iPhi > -1){ //Case: Specific (iEta,iPhi) sector
        ret_Name = strInputPrefix + "_iEta" + getString(iEta) + "iPhi" + getString(iPhi) + "_" + strInputName;
    } //End Case: Specific (iEta,iPhi) sector
    else if (iEta > -1){
        ret_Name = strInputPrefix + "_iEta" + getString(iEta) + "_" + strInputName;
    }
    else{ //Case: iEta Sector, sum over sector's iPhi
        ret_Name = strInputPrefix + "_Summary_" + strInputName;
    } //End Case: iEta Sector, sum over sector's iPhi
    
    return ret_Name;
} //End AnalyzeResponseUniformity::getNameByIndex()

//Searches the input TF1 for a parameter with meaning given by strParam and stored in HistoSetup
//This parameter is then returned to the user
float AnalyzeResponseUniformity::getParam( shared_ptr<TF1> fitInput, InfoFit & setupFit, std::string strParam ){
    //Variable Declaration
    int iParamPos = -1;
    
    float ret_Val = -1;
    
    vector<string>::iterator iterParamMeaning = std::find(setupFit.m_vec_strFit_ParamMeaning.begin(), setupFit.m_vec_strFit_ParamMeaning.end(), strParam);
    
    if ( iterParamMeaning != setupFit.m_vec_strFit_ParamMeaning.end() ) { //Case: Parameter Found!!!
        
        iParamPos = std::distance(setupFit.m_vec_strFit_ParamMeaning.begin(), iterParamMeaning);
        
        ret_Val = fitInput->GetParameter(iParamPos);
    } //End Case: Parameter Found!!!
    else{ //Case: Parameter NOT Found
        printClassMethodMsg("AnalyzeResponseUniformity","getParam",("Error! - I did not find your requested parameter: " + strParam + "!\n").c_str() );
        printClassMethodMsg("AnalyzeResponseUniformity","getPeakPos","\tPlease Cross-check input analysis config file.\n");
    } //End Case: Parameter NOT Foun
    
    return ret_Val;
} //End AnalyzeResponseUniformity::getParam

//Searches the input TF1 for a parameter with meaning given by strParam and stored in HistoSetup
//The error on this parameter is then returned to the user
float AnalyzeResponseUniformity::getParamError( shared_ptr<TF1> fitInput, Plotter::InfoFit & setupFit, std::string strParam ){
    //Variable Declaration
    int iParamPos = -1;
    
    float ret_Val = -1;
    
    vector<string>::iterator iterParamMeaning = std::find(setupFit.m_vec_strFit_ParamMeaning.begin(), setupFit.m_vec_strFit_ParamMeaning.end(), strParam);
    
    if ( iterParamMeaning != setupFit.m_vec_strFit_ParamMeaning.end() ) { //Case: Parameter Found!!!
        
        iParamPos = std::distance(setupFit.m_vec_strFit_ParamMeaning.begin(), iterParamMeaning);
        
        ret_Val = fitInput->GetParError(iParamPos);
    } //End Case: Parameter Found!!!
    else{ //Case: Parameter NOT Found
        printClassMethodMsg("AnalyzeResponseUniformity","getParamError",("Error! - I did not find your requested parameter: " + strParam + "!\n").c_str() );
        printClassMethodMsg("AnalyzeResponseUniformity","getParamError","\tPlease Cross-check input analysis config file.\n");
    } //End Case: Parameter NOT Foun
    
    return ret_Val;
} //End AnalyzeResponseUniformity::getParamError

//Given an input histogram and TSpectrum returns a numeric value based on the input keyword; supported keywords are "AMPLITUDE,MEAN,PEAK,SIGMA"
float AnalyzeResponseUniformity::getValByKeyword(string strInputKeyword, shared_ptr<TH1F> hInput, TSpectrum &specInput){
    
    //Try to automatically assign a value
    if ( 0 == strInputKeyword.compare("AMPLITUDE") ) { //Case: Histo Amplitude
        return hInput->GetBinContent( hInput->GetMaximumBin() );
    } //End Case: Histo Amplitude
    else if( 0 == strInputKeyword.compare("FWHM") ){ //Case: Full Width Half Max
        //Right now as estimate we just use the histo RMS
        return hInput->GetRMS();
    } //End Case: Full Width Half Max
    else if( 0 == strInputKeyword.compare("HWHM") ){ //Case: Half Width Half Max
        //Right now as estimate we just use half the histo RMS
        return 0.5 * hInput->GetRMS();
    } //End Case: Half Width Half Max
    else if (0 == strInputKeyword.compare("MEAN") ) { //Case: Histo Mean
        return hInput->GetMean();
    } //End Case: Histo Mean
    else if ( 0 == strInputKeyword.compare("PEAK") ){ //Case: Histo Peak
        Double_t *dPeakPos = specInput.GetPositionX();
        
        return dPeakPos[0];
    } //End Case: Histo Peak
    else if (0 == strInputKeyword.compare("SIGMA") ) { //Case: Histo RMS
        return hInput->GetRMS();
    } //End Case: Histo RMS
    else{ //Case: manual assignment
        printClassMethodMsg("AnalyzeResponseUniformity","getValByKeyword","Error! Input Keyword Not Recognized");
        printClassMethodMsg("AnalyzeResponseUniformity","getValByKeyword", ("\tGiven: " + strInputKeyword ).c_str() );
        printClassMethodMsg("AnalyzeResponseUniformity","getValByKeyword","\tRecognized Keywords:\n");
        
        for (int i=0; i < vec_strSupportedKeywords.size(); ++i) {
            printClassMethodMsg("AnalyzeResponseUniformity","getValByKeyword", vec_strSupportedKeywords[i].c_str() );
        }
        
        printClassMethodMsg("AnalyzeResponseUniformity","getValByKeyword","\tUndefined Behavior May Occur");
        
        return -1e12;
    } //End Case: manual assignment
} //End AnalyzeResponseUniformity::getValByKeyword()
