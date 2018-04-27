//
//  AnalyzeResponseUniformity.h
//  
//
//  Created by Brian L Dorney on 28/01/16.
//
//

#ifndef ____AnalyzeResponseUniformity__
#define ____AnalyzeResponseUniformity__

//C++ Includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <stdio.h>
#include <string>
#include <vector>
#include "TLatex.h"

//Framework Includes
//#include "DetectorMPGD.h"
#include "FrameworkBase.h"
//#include "ParameterLoaderDetector.h"
#include "PlotterUtilityTypes.h"
#include "TimingUtilityFunctions.h"
#include "TimingUtilityTypes.h"
#include "UniformityUtilityOperators.h"
#include "UniformityUtilityTypes.h"

//ROOT Includes
#include "TDirectory.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
//#include "TPolyMarker.h"
#include "TROOT.h"
#include "TSpectrum.h"

//C++ Mathematical Expression Toolkit Library
//By Arash Partow; http://www.partow.net/programming/exprtk/index.html
#include "exprtk.hpp"

namespace QualityControl {
    namespace Uniformity {
        //Used for evaluating fit limits given by user
        typedef exprtk::symbol_table<float> symbol_table_t;
        typedef exprtk::expression<float> expression_t;
        typedef exprtk::parser<float> parser_t;
        
        class AnalyzeResponseUniformity : public FrameworkBase {
            
        public:
            //Constructors
            //------------------------------------------------------------------------------------------------------------------------------------------
            //Default
            AnalyzeResponseUniformity();
            
            //Set the detector and anlysis parameters at construction
            AnalyzeResponseUniformity(AnalysisSetupUniformity inputSetup, DetectorMPGD & inputDet);
            
            //Actions - Methods that Do Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            
            //Getters - Methods that Get (i.e. Return) Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            //virtual DetectorMPGD getDetector(){ return detMPGD; };
            
            //Printers - Methods that Print Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            
            //Setters - Methods that Set Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            //Sets the Analysis Setup
            //virtual void setAnalysisParameters(Uniformity::AnalysisSetupUniformity inputSetup){ aSetup = inputSetup; return; };
            
            //Sets the Detector
            //virtual void setDetector(Uniformity::DetectorMPGD & inputDet){ detMPGD = inputDet; return; };
            
        protected:
            //Actions - Methods that Do Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            //Calculates Summary Statistics
            void calcStatistics(SummaryStatistics &inputStatObs, std::multiset<float> &mset_fInputObs, std::string strObsName);
            
            //Determines if a fit is "good"
            //Good fits have:
            //		fit parameter != 0 (hold that thought...better implementation needed?)
            //		fit parameter != parameter limit
            //		Percent Error on fit parameter (sigma/value) <= 0.1
            bool isQualityFit(std::shared_ptr<TF1> fitInput);

            //As above but for a specific parameter iPar
            bool isQualityFit(std::shared_ptr<TF1> fitInput, int iPar);

            //Getters - Methods that Get (i.e. Return) Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            
            //Returns a fit whose parameters match those defined in the AnalysisSetupUniformity
            //TF1 getFit(int iEta, int iPhi, int iSlice, Timing::HistoSetup & setupHisto, std::shared_ptr<TH1F> hInput, TSpectrum &specInput );
            TF1 getFit(int iEta, int iPhi, int iSlice, Plotter::InfoFit & setupFit, std::shared_ptr<TH1F> hInput, TSpectrum &specInput );
            
            //strInput is understood to be a number or an algebraic expression
            //If strInput contains a substr matching to one or more elements of vec_strSupportedKeywords it is treated as an algebraic expression
            //If so it parses the expression and returns a numeric value
            //If not it is treated as a number and attempts to convert strInput directly
            float getParsedInput(std::string &strInputExp, std::shared_ptr<TH1F> hInput, TSpectrum &specInput);
            
            //Returns a TGraph Errors whose parameters match those defined in the input HistoSetup object
            TGraphErrors getGraph(int iEta, int iPhi, Timing::HistoSetup &setupHisto);
            
            //Returns a histogram whose parmeters match those defined in the input HistoSetup object
            TH1F getHistogram(int iEta, int iPhi, Timing::HistoSetup &setupHisto);
            
            //Returns a histogram whose parmeters match those defined in the input HistoSetup object
            TH2F getHistogram2D(int iEta, int iPhi, Timing::HistoSetup &setupHisto_X, Timing::HistoSetup &setupHisto_Y);
            
            //Formats a given input string such that it follows the iEta, iPhi, iSlice naming convention
            std::string getNameByIndex(int iEta, int iPhi, int iSlice, std::string & strInputPrefix, std::string & strInputName);
            std::string getNameByIndex(int iEta, int iPhi, int iSlice, const std::string & strInputPrefix, const std::string & strInputName);
            
            //Searchs the input fit for the given variable (strParam); returns it
            float getParam( std::shared_ptr<TF1> fitInput, Plotter::InfoFit & setupFit, std::string strParam );
            
            //Searchs the input fit for the error on the given variable (strParam); returns it
            float getParamError( std::shared_ptr<TF1> fitInput, Plotter::InfoFit & setupFit, std::string strParam );
            
            //Given an input histogram and TSpectrum returns a numeric value based on the input keyword; supported keywords are "AMPLITUDE,MEAN,PEAK,SIGMA"
            float getValByKeyword(std::string strInputKeyword, std::shared_ptr<TH1F> hInput, TSpectrum &specInput);
            
            //Data Members
            //------------------------------------------------------------------------------------------------------------------------------------------
            std::string strAnalysisName;
            
            const std::vector<std::string> vec_strSupportedKeywords = {"AMPLITUDE", "FWHM", "HWHM", "MEAN","PEAK","SIGMA"}; //Supported Keywords for fit setup
            
            //AnalysisSetupUniformity aSetup; //Container to define the analysis setup
            
            //DetectorMPGD detMPGD; //Link to header file in AnalyzeResponseUniformity.cpp; Detector object
        }; //End class AnalyzeResponseUniformity
    } //End namespace Uniformity
} //End namespace QualityControl

#endif /* defined(____AnalyzeResponseUniformity__) */
