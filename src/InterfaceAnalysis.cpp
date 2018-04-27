//
//  InterfaceAnalysis.cpp
//  
//
//  Created by Brian L Dorney on 03/05/16.
//
//

#include "InterfaceAnalysis.h"

using std::cout;
using std::endl;
using std::map;
using std::string;

using namespace QualityControl::Uniformity;

//Default Constructor
InterfaceAnalysis::InterfaceAnalysis(){
    //bInitialized = false;
    bVerboseMode = false;
} //End Default Constructor

//Performs the analysis on the detMPGD object defined by rSetup and aSetup on the input files stored in vec_pairedRunList
void InterfaceAnalysis::analyzeInput(){
    //Variable Declaration
    
    //Evaluate the two possible usage cases
    //  Case 1: input created by amoreSRS (rSetup.bInputFromFrmwrk == false)
    //  Case 2: input created by framework (rSetup.bInputFromFrmwrk == true)
    if (rSetup.bInputFromFrmwrk) { //Case: Input created by CMS_GEM_Analysis_Framework
        analyzeInputFrmwrk();
    } //End Case: Input created by CMS_GEM_Analysis_Framework
    else{ //Case: Input created by amoreSRS
        analyzeInputAmoreSRS();
    } //Case: Input created by amoreSRS
    
    return;
} //End InterfaceAnalysis::analyzeInput()

//Runs the analysis framework on input created by amoreSRS
void InterfaceAnalysis::analyzeInputAmoreSRS(){
    //TFile does not automatically own histograms
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    
    ReadoutSectorEta etaSector;
    
    string strTempRunName;
    
    TFile *file_ROOTInput, *file_ROOTOutput_All, *file_ROOTOutput_Single;
    
    //Debugging
    //cout<<"InterfaceAnalysis::analyzeInputAmoreSRS(): detMPGD.getName() = " << detMPGD.getName() << endl;
    
    //Print the detector info & analysis parameters to the user
    //------------------------------------------------------
    if (bVerboseMode) { //Case: Verbose Printing
        
        //Detector info
        cout<<"iEta\tY_Pos\tWidth\tX_Low\tX_High\tX_Low\tX_High\tX_Low\tX_High\n";
        for (int i=1; i <= detMPGD.getNumEtaSectors(); ++i) { //Loop through Detector's ReadoutSectorEta objects
            cout<<i<<"\t";
            
            etaSector = detMPGD.getEtaSector(i);
            
            cout<<etaSector.fPos_Y<<"\t"<<etaSector.fWidth<<"\t";
            
            for (auto iterPhi = etaSector.map_sectorsPhi.begin(); iterPhi != etaSector.map_sectorsPhi.end(); ++iterPhi) { //Loop through ReadoutSectorEta's ReadoutSectorPhi objects
                cout<<(*iterPhi).second.fPos_Xlow<<"\t"<<(*iterPhi).second.fPos_Xhigh<<"\t";
            } //End Loop through ReadoutSectorEta's ReadoutSectorPhi objects
            
            cout<<endl;
        } //End Loop through Detector's ReadoutSectorEta objects
        
        //Analysis parameters
        cout<<"-----------------Hit Parameters-----------------\n";
        cout<<"Hit ADC, Min = " << aSetup.selHit.iCut_ADCNoise << endl;
        cout<<"Hit ADC, Max = " << aSetup.selHit.iCut_ADCSat << endl;
        cout<<"Hit Multi, Min = " << aSetup.selHit.iCut_MultiMin << endl;
        cout<<"Hit Multi, Max = " << aSetup.selHit.iCut_MultiMax << endl;
        cout<<"Hit Time, Min = " << aSetup.selHit.iCut_TimeMin << endl;
        cout<<"Hit Time, Max = " << aSetup.selHit.iCut_TimeMax << endl;
        cout<<"-----------------Cluster Parameters-----------------\n";
        cout<<"Clust ADC, Min = " << aSetup.selClust.iCut_ADCNoise << endl;
        cout<<"Clust Size, Min = " << aSetup.selClust.iCut_SizeMin << endl;
        cout<<"Clust Size, Max = " << aSetup.selClust.iCut_SizeMax << endl;
        cout<<"Clust Time, Min = " << aSetup.selClust.iCut_TimeMin << endl;
        cout<<"Clust Time, Max = " << aSetup.selClust.iCut_TimeMax << endl;
    } //End Case: Verbose Printing
    
    //Debugging
    //cout<<"InterfaceAnalysis::analyzeInputAmoreSRS(): detMPGD.getName() = " << detMPGD.getName() << endl;
    
    //Loop over input files
    //------------------------------------------------------
    for (int i=0; i < vec_pairedRunList.size(); ++i) { //Loop over vec_pairedRunList
        //Wipe physics objects from previous file (prevent double counting)
        detMPGD.resetPhysObj();
        
        //Open this run's root file & check to see if data file opened successfully
        //------------------------------------------------------
        file_ROOTInput = new TFile(vec_pairedRunList[i].second.c_str(),"READ","",1);
        
        if ( !file_ROOTInput->IsOpen() || file_ROOTInput->IsZombie() ) { //Case: failed to load ROOT file
            perror( ("InterfaceAnalysis::analyzeInputAmoreSRS() - error while opening file: " + vec_pairedRunList[i].second ).c_str() );
            Timing::printROOTFileStatus(file_ROOTInput);
            std::cout << "Skipping!!!\n";
            
            continue;
        } //End Case: failed to load ROOT file
        
        //Debugging
        //cout<<"InterfaceAnalysis::analyzeInputAmoreSRS(): Pre Hit Ana detMPGD.getName() = " << detMPGD.getName() << endl;
        
        //Hit Analysis
        //------------------------------------------------------
        //Force the hit analysis if the user requested cluster reconstruction
        if ( rSetup.bAnaStep_Hits ) { //Case: Hit Analysis
            //Hit Selection
            hitSelector.setRunNum(vec_pairedRunList[i].first);
            hitSelector.setHits(file_ROOTInput, detMPGD, aSetup);
            
            if (bVerboseMode) { //Print Number of Selected Hits to User
                cout<<vec_pairedRunList[i].second << " has " << detMPGD.getHits().size() << " hits passing selection" << endl;
            } //End Print Number of Selected Hits to User
            
            //Load the required input parameters
            if (i == 0) { hitAnalyzer.setAnalysisParameters(aSetup); } //Fixed for all runs
            
            if (i == 0 || rSetup.bMultiOutput) { hitAnalyzer.initHistosHits(detMPGD); }
            
            //Hit Analysis
            hitAnalyzer.setRunNum(vec_pairedRunList[i].first);
            hitAnalyzer.fillHistos(detMPGD);
        } //End Case: Hit Analysis
        
        //Cluster Analysis
        //------------------------------------------------------
        if ( rSetup.bAnaStep_Clusters ) { //Case: Cluster Analysis
            //Cluster Selection
            clustSelector.setRunNum(vec_pairedRunList[i].first);
            clustSelector.setClusters(file_ROOTInput, detMPGD, aSetup);
            
            if (bVerboseMode) { //Print Number of Selected Clusters to User
                cout<<vec_pairedRunList[i].second << " has " << detMPGD.getClusters().size() << " clusters passing selection" << endl;
            } //End Print Number of Selected Clusters to User
            
            //Load the required input parameters
            if (i == 0) { clustAnalyzer.setAnalysisParameters(aSetup); } //Fixed for all runs
            
            //Initialize the cluster histograms for all runs; if this is the first run
            if (i == 0 || rSetup.bMultiOutput) {
                clustAnalyzer.initGraphsClusters(detMPGD);
                clustAnalyzer.initHistosClusters(detMPGD);
            }
            
            //Initialize the cluster histograms specific to this run
            clustAnalyzer.initHistosClustersByRun(vec_pairedRunList[i].first, detMPGD);
            
            //Cluster Analysis
            clustAnalyzer.setRunNum(vec_pairedRunList[i].first);
            clustAnalyzer.fillHistos(detMPGD);
        } //End Case: Cluster Analysis
        
        //User requests multiple output files?
        //------------------------------------------------------
        if (rSetup.bMultiOutput) { //Case: User wants one output file per input file
            //Setup the name of the output file
            strTempRunName = vec_pairedRunList[i].second;
            
            if ( strTempRunName.find("dataTree.root") != string::npos ){
                strTempRunName.erase(strTempRunName.find("dataTree.root"), strTempRunName.length() - strTempRunName.find("dataTree.root") );
                strTempRunName = strTempRunName + "Ana.root";
                //replaceSubStr1WithSubStr2(strTempRunName, "dataTree.root", "Ana.root");
            } //End Case: Input Tree File
            else if ( strTempRunName.find(".root") != string::npos ) {
                strTempRunName.erase(strTempRunName.find(".root"), strTempRunName.length() - strTempRunName.find(".root") );
                strTempRunName = strTempRunName + "Ana.root";
                //replaceSubStr1WithSubStr2(strTempRunName, ".root", "Ana.root");
            } //End Case: Other ROOT file
            
            //Create the file
            file_ROOTOutput_Single = new TFile(strTempRunName.c_str(), rSetup.strFile_Output_Option.c_str(),"",1);
            
            if ( !file_ROOTOutput_Single->IsOpen() || file_ROOTOutput_Single->IsZombie() ) { //Case: failed to load ROOT file
                perror( ("InterfaceAnalysis::analyzeInputAmoreSRS() - error while opening file: " + strTempRunName ).c_str() );
                Timing::printROOTFileStatus(file_ROOTOutput_Single);
                std::cout << "Skipping!!!\n";
                
                //Close the file & delete pointer before the next iter
                file_ROOTInput->Close();
                delete file_ROOTInput;
                
                //Move to next iteration
                continue;
            } //End Case: failed to load ROOT file
            
            //Store the results
            storeResults(file_ROOTOutput_Single, strTempRunName);
            
            //Close the file before the next iteration
            file_ROOTOutput_Single->Close();
            delete file_ROOTOutput_Single;
        } //End Case: User wants one output file per input file
        
        //Close the file & delete pointer before the next iter
        //------------------------------------------------------
        file_ROOTInput->Close();
        delete file_ROOTInput;
    } //End Loop over vec_pairedRunList
    
    //Debugging
    //cout<<"InterfaceAnalysis::analyzeInputAmoreSRS(): Pre Storing detMPGD.getName() = " << detMPGD.getName() << endl;
    
    //Create the summary TFile (only if multi file output is NOT requested)
    //------------------------------------------------------
    if(!rSetup.bMultiOutput){ //Case: Create Summary TFile
        file_ROOTOutput_All = new TFile(rSetup.strFile_Output_Name.c_str(), rSetup.strFile_Output_Option.c_str(),"",1);
        
        if ( !file_ROOTOutput_All->IsOpen() || file_ROOTOutput_All->IsZombie() ) { //Case: failed to load ROOT file
            perror( ("InterfaceAnalysis::analyzeInputAmoreSRS() - error while opening file: " + rSetup.strFile_Output_Name ).c_str() );
            Timing::printROOTFileStatus(file_ROOTOutput_All);
            std::cout << "Exiting!!!\n";
            
            return;
        } //End Case: failed to load ROOT file
        
        //Store the results
        storeResults(file_ROOTOutput_All, rSetup.strFile_Output_Name);
        
        //Close the file
        file_ROOTOutput_All->Close();
    } //End Case: Create Summary TFile
    
    return;
} //End InterfaceAnalysis::analyzeInputAmoreSRS()

//Runs the analysis framework on input created by the CMS_GEM_AnalysisFramework
void InterfaceAnalysis::analyzeInputFrmwrk(){
    //TFile does not automatically own histograms
    TH1::AddDirectory(kFALSE);
    
    //Variable Declaration
    TFile *file_ROOTInput, *file_ROOTOutput_All, *file_ROOTOutput_Single;
    
    //Loop over input files
    for (int i=0; i < vec_pairedRunList.size(); ++i) { //Loop over input files
        //Open this run's root file & check to see if data file opened successfully
        //------------------------------------------------------
        file_ROOTInput = new TFile(vec_pairedRunList[i].second.c_str(),"READ","",1);
        
        if ( !file_ROOTInput->IsOpen() || file_ROOTInput->IsZombie() ) { //Case: failed to load ROOT file
            perror( ("InterfaceAnalysis::analyzeInputFrmwrk() - error while opening file: " + vec_pairedRunList[i].second ).c_str() );
            Timing::printROOTFileStatus(file_ROOTInput);
            std::cout << "Skipping!!!\n";
            
            continue;
        } //End Case: failed to load ROOT file
        
        //Hit Analysis
        //------------------------------------------------------
        //Force the hit analysis if the user requested cluster reconstruction
        if ( rSetup.bAnaStep_Hits ) { //Case: Hit Analysis
            
            //Placeholder
            
        } //End Case: Hit Analysis
        
        //Cluster Analysis
        //------------------------------------------------------
        if ( rSetup.bAnaStep_Clusters ) { //Case: Cluster Analysis
            //Load the required input parameters
            if (i == 0) { clustAnalyzer.setAnalysisParameters(aSetup); } //Fixed for all runs
            
            //Load previous cluster histograms & setup the detector
            clustAnalyzer.loadHistosFromFile(rSetup.strFile_Config_Map, file_ROOTInput);
            detMPGD = clustAnalyzer.getDetector();
            
            //Initialize Graphs
            clustAnalyzer.initGraphsClusters(detMPGD);
        } //End Case: Cluster Analysis
        
        //Store the Output
        //When running over framework input only the one output file per input file mode is supported
        //  e.g. there is no "aggregate file"
        //------------------------------------------------------
        //Setup the name of the output file
        string strTempRunName = vec_pairedRunList[i].second;
        
        if ( strTempRunName.find("Ana.root") != string::npos ) {
            strTempRunName.erase(strTempRunName.find("Ana.root"), strTempRunName.length() - strTempRunName.find("Ana.root") );
            strTempRunName = strTempRunName + "NewAna.root";
            //replaceSubStr1WithSubStr2(strTempRunName,"Ana.root","NewAna.root");
        } //End Case: Other ROOT file
        
        //Create the file
        file_ROOTOutput_Single = new TFile(strTempRunName.c_str(), rSetup.strFile_Output_Option.c_str(),"",1);
        
        if ( !file_ROOTOutput_Single->IsOpen() || file_ROOTOutput_Single->IsZombie() ) { //Case: failed to load ROOT file
            perror( ("InterfaceAnalysis::analyzeInputFrmwrk() - error while opening file: " + strTempRunName ).c_str() );
            Timing::printROOTFileStatus(file_ROOTOutput_Single);
            std::cout << "Skipping!!!\n";
            
            //Close the file & delete pointer before the next iter
            file_ROOTInput->Close();
            delete file_ROOTInput;
            
            //Move to next iteration
            continue;
        } //End Case: failed to load ROOT file
        
        //Store the results
        //cout<<"InterfaceAnalysis::analyzeInputFrmwrk(): file_ROOTOutput_Single = " << file_ROOTOutput_Single << endl;
        storeResults(file_ROOTOutput_Single, strTempRunName);
        
        //Close the files & delete pointers before the next iter
        //------------------------------------------------------
        file_ROOTInput->Close();
        file_ROOTOutput_Single->Close();
        
        delete file_ROOTInput;
        delete file_ROOTOutput_Single;
    } //End Loop over input files
    
    return;
} //End InterfaceAnalysis::analyzeInputFrmwrk()

void InterfaceAnalysis::storeResults(TFile * file_Results, string strFileName){
    //Variable Declaration
    map<string,string> map_clust_ObsAndDrawOpt; //Cluster observables & draw option
    map<string,string> map_hit_ObsAndDrawOpt;   //as above but for hits
    map<string,string> map_res_ObsAndDrawOpt;   //as above but for results (e.g. fits)
    
    string strTempRunName;
    
    //cout<<"InterfaceAnalysis::storeResults(): file_Results = " << file_Results << endl;
    //cout<<"InterfaceAnalysis::storeResults(): detMPGD.getName() = " << detMPGD.getName() << endl;
    
    //Store Histograms After Analyzing all input files
    //------------------------------------------------------
    cout<<"<<<<<<<< Storing Hit Histograms >>>>>>>>\n";
    if ( rSetup.bAnaStep_Hits){ //Case: Hits
        //Store Hit Histograms
        hitAnalyzer.storeHistos(file_Results, detMPGD);
        
        //Setup name for strip list
        strTempRunName = strFileName;
        if ( strTempRunName.find("dataTree.root") != string::npos ){
            strTempRunName.erase(strTempRunName.find("dataTree.root"), strTempRunName.length() - strTempRunName.find("dataTree.root") );
            strTempRunName = strTempRunName + "_DeadStripList.txt";
            //replaceSubStr1WithSubStr2(strTempRunName, "dataTree.root", "Ana.root");
        } //End Case: Input Tree File
        else if ( strTempRunName.find(".root") != string::npos ) {
            strTempRunName.erase(strTempRunName.find(".root"), strTempRunName.length() - strTempRunName.find(".root") );
            strTempRunName = strTempRunName + "_DeadStripList.txt";
            //replaceSubStr1WithSubStr2(strTempRunName, ".root", "Ana.root");
        } //End Case: Other ROOT file
        
        //Store Dead Strip List
        cout<<"<<<<<<<< Making Dead Strip List >>>>>>>>\n";
        hitAnalyzer.findDeadStrips(detMPGD, strTempRunName);
    } //End Case: Hits
    
    cout<<"<<<<<<<< Storing Cluster Histograms >>>>>>>>\n";
    if ( rSetup.bAnaStep_Clusters) clustAnalyzer.storeHistos(file_Results, detMPGD);
    
    //Fit Histograms After Analyzing all input files
    //------------------------------------------------------
    if ( rSetup.bAnaStep_Fitting){ //Case: Fitting Stored Distributions
        if ( rSetup.bAnaStep_Clusters){ //Case: Cluster Analysis
            cout<<"<<<<<<<< Fitting Cluster Histograms >>>>>>>>\n";
            clustAnalyzer.fitHistos(detMPGD);
            
            cout<<"<<<<<<<< Storing Cluster Fits >>>>>>>>\n";
            clustAnalyzer.storeFits(file_Results, detMPGD);
        } //End Case: Cluster Analysis
    } //End Case: Fitting Stored Distributions
    
    //Visualize Results
    //------------------------------------------------------
    if ( rSetup.bAnaStep_Visualize ) { //Case: Visualize Output
        visualizeUni.setAnalysisParameters(aSetup);
        visualizeUni.setAutoSaveCanvas( rSetup.bVisPlots_AutoSaving );
        visualizeUni.setDetector(detMPGD);
        
        if (rSetup.bAnaStep_Hits) { //Case: Hit Analysis
            cout<<"<<<<<<<< Making Hit Summary Plots >>>>>>>>\n";
            
            map_hit_ObsAndDrawOpt["HitADC"]="E1";
            map_hit_ObsAndDrawOpt["HitMulti"]="E1";
            map_hit_ObsAndDrawOpt["HitTime"]="E1";
            
            visualizeUni.storeCanvasHistoSegmented(file_Results, "HitPos", "E1", rSetup.bVisPlots_PhiLines);
            visualizeUni.storeListOfCanvasesHistoSegmented(file_Results, map_hit_ObsAndDrawOpt, false);
        } //End Case: Hit Analysis
        
        if (rSetup.bAnaStep_Clusters) { //Case: Cluster Analysis
            cout<<"<<<<<<<< Making Cluster Summary Plots >>>>>>>>\n";
            
            map_clust_ObsAndDrawOpt["ClustADC"]="E1";
            map_clust_ObsAndDrawOpt["ClustMulti"]="E1";
            map_clust_ObsAndDrawOpt["ClustSize"]="E1";
            map_clust_ObsAndDrawOpt["ClustTime"]="E1";
            
            visualizeUni.storeCanvasHistoSegmented(file_Results, "ClustPos", "E1", rSetup.bVisPlots_PhiLines);
            visualizeUni.storeListOfCanvasesHistoSegmented(file_Results, map_clust_ObsAndDrawOpt, false);
            
            if( !rSetup.bInputFromFrmwrk ) {
                cout<<"<<<<<<<< Making Cluster Run History Summary Plots >>>>>>>>\n";
                
                visualizeUni.storeCanvasHisto2DHistorySegmented(file_Results, "HistoryClustADC", "COLZ", true); //Eta Level
                visualizeUni.storeCanvasHisto2DHistorySegmented(file_Results, "HistoryClustADC", "COLZ", false); //Phi Level
                visualizeUni.storeCanvasHisto2DHistorySegmented(file_Results, "HistoryClustTime", "COLZ", true); //Eta Level
                visualizeUni.storeCanvasHisto2DHistorySegmented(file_Results, "HistoryClustTime", "COLZ", false); //Phi Level
            }
            
            if (rSetup.bAnaStep_Fitting) { //Case: Fitting
                cout<<"<<<<<<<< Making Uniformity Summary Plots >>>>>>>>\n";
                
                map_res_ObsAndDrawOpt["ResponseFitChi2"]="APE1";
                map_res_ObsAndDrawOpt["ResponseFitPkPos"]="APE1";
                map_res_ObsAndDrawOpt["ResponseFitPkRes"]="APE1";
                
                visualizeUni.storeCanvasData(file_Results, "ResponseFitChi2", "E1",false);
                visualizeUni.storeCanvasData(file_Results, "ResponseFitPkPos", "E1",false);
                visualizeUni.storeCanvasData(file_Results, "ResponseFitPkPos", "E1",true);
                visualizeUni.storeCanvasData(file_Results, "ResponseFitPkRes", "E1",false);
                visualizeUni.storeCanvasFits(file_Results, "COLZTEXT");
                
                visualizeUni.storeListOfCanvasesGraph(file_Results,map_res_ObsAndDrawOpt, rSetup.bVisPlots_PhiLines);
                visualizeUni.storeCanvasGraph2D(file_Results,"ResponseFitPkPos","TRI2Z",false);
                visualizeUni.storeCanvasGraph2D(file_Results,"ResponseFitPkPos","TRI2Z",true);	//Normalized version
                visualizeUni.storeCanvasGraph2D(file_Results,"ResponseFitPkRes","TRI2Z",false);
                visualizeUni.storeCanvasGraph2D(file_Results,"ResponseFitPkRes","TRI2Z",true);	//Normalized version
            } //End Case: Fitting
        } //End Case: Cluster Analysis
    } //End Case: Visualize Output
    
    return;
} //End InterfaceAnalysis::storeResults()
