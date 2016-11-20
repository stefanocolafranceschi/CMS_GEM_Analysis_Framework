//
//  cosmicStandTest.cpp
//  
//
//  Created by Brian L Dorney on 27/10/16.
//
//

//C++ Includes
#include <memory>
#include <stdio.h>
#include <string>
#include <vector>

//Framework Includes
#include "DetectorTiming.h"
#include "SelectorTiming.h"
#include "TestStandVME.h"
#include "TimingEvent.h"
#include "TimingRunSetup.h"
#include "TimingUtilityTypes.h"

//ROOT Includes

//namespaces
using std::cout;
using std::endl;
using std::shared_ptr;
using std::string;
using std::vector;

using QualityControl::Timing::DetectorTiming;
using QualityControl::Timing::DetectorMatrix;
using QualityControl::Timing::EventDigi;
using QualityControl::Timing::EventReco;
using QualityControl::Timing::HwVMEBoard;
using QualityControl::Timing::DetType;
using QualityControl::Timing::RunSetup;
using QualityControl::Timing::SelectorTiming;
using QualityControl::Timing::TDCDataDigi;
using QualityControl::Timing::TestStandVME;
using QualityControl::Timing::VMETypes;

using QualityControl::Timing::kDUT;
using QualityControl::Timing::kTrigger;
using QualityControl::Timing::kVMETDC;

void printEventDetMatrix(DetectorMatrix & inputDetMatrix){
    for (auto iterMatrix = inputDetMatrix.m_map_detectors.begin(); iterMatrix != inputDetMatrix.m_map_detectors.end(); ++iterMatrix) {
        for (auto iterDet = (*iterMatrix).second->getDetectorPtrBegin(); iterDet != (*iterMatrix).second->getDetectorPtrEnd(); ++iterDet) {
            cout<<"----" << (*iterDet).second->getName() << "----\n";
            cout<<"Ch\tData\n";
            for (int i=0; i < (*iterDet).second->getNChan() ; ++i) {
                cout<<i<<"\t"<<(*iterDet).second->getChanData(i)<<endl;
            }
        } //End Loop Over Detectors in Trigger
    } //End Loop Over Trigger Detectors
} //End printEventDetMatrix

void printEvent(EventReco & evtReco){
    cout<<"Evt = " << evtReco.m_uiEvtCount << endl;
    
    cout<<"=======Trigger Info=======\n";
    printEventDetMatrix(evtReco.m_detMatrix_Trigger);
    
    cout<<"=======DUT Info=======\n";
    printEventDetMatrix(evtReco.m_detMatrix_DUTs);
    
    return;
} //End printEvent()

int main(){
    TestStandVME cosmicStand;
    
    //Top PMT's
    shared_ptr<DetectorTiming> det_PMT1(new DetectorTiming);
    shared_ptr<DetectorTiming> det_PMT2(new DetectorTiming);
    shared_ptr<DetectorTiming> det_PMT3(new DetectorTiming);
    shared_ptr<DetectorTiming> det_PMT4(new DetectorTiming);
    shared_ptr<DetectorTiming> det_PMT5(new DetectorTiming);

    //Bottom PMT's
    shared_ptr<DetectorTiming> det_PMT6(new DetectorTiming);
    shared_ptr<DetectorTiming> det_PMT8(new DetectorTiming);
    shared_ptr<DetectorTiming> det_PMT9(new DetectorTiming);
    shared_ptr<DetectorTiming> det_PMT11(new DetectorTiming);
    shared_ptr<DetectorTiming> det_PMT12(new DetectorTiming);
    
    //Detectors
    shared_ptr<DetectorTiming> det_GE11Long2(new DetectorTiming);
    shared_ptr<DetectorTiming> det_GE11Long4(new DetectorTiming);
    
    det_PMT1->setName("PMT1");  det_PMT1->setBaseAddress( "400C" ); det_PMT1->setChannelMapVME2Det(0,0);
    det_PMT2->setName("PMT2");  det_PMT2->setBaseAddress( "400C" ); det_PMT2->setChannelMapVME2Det(1,0);
    det_PMT3->setName("PMT3");  det_PMT3->setBaseAddress( "400C" ); det_PMT3->setChannelMapVME2Det(2,0);
    det_PMT4->setName("PMT4");  det_PMT4->setBaseAddress( "400C" ); det_PMT4->setChannelMapVME2Det(3,0);
    det_PMT5->setName("PMT5");  det_PMT5->setBaseAddress( "400C" ); det_PMT5->setChannelMapVME2Det(4,0);
    
    det_PMT6->setName("PMT6");  det_PMT6->setBaseAddress( "400C" ); det_PMT6->setChannelMapVME2Det(5,0);
    det_PMT8->setName("PMT8");  det_PMT8->setBaseAddress( "400C" ); det_PMT8->setChannelMapVME2Det(7,0);
    det_PMT9->setName("PMT9");  det_PMT9->setBaseAddress( "400C" ); det_PMT9->setChannelMapVME2Det(8,0);
    det_PMT11->setName("PMT11"); det_PMT11->setBaseAddress( "400C" );det_PMT11->setChannelMapVME2Det(10,0);
    det_PMT12->setName("PMT12"); det_PMT12->setBaseAddress( "400C" );det_PMT12->setChannelMapVME2Det(12,0);
    
    det_GE11Long2->setName("GE11-VII-L-CERN-0002"); //3/2/T
    det_GE11Long2->setBaseAddress("0011");
    det_GE11Long2->setChannelMapVME2Det(8, 1);
    det_GE11Long2->setChannelMapVME2Det(9, 1);
    det_GE11Long2->setChannelMapVME2Det(10, 1);
    det_GE11Long2->setChannelMapVME2Det(11, 2);
    det_GE11Long2->setChannelMapVME2Det(12, 3);
    det_GE11Long2->setChannelMapVME2Det(13, 4);

    det_GE11Long4->setName("GE11-VII-L-CERN-0004");  //3/2/B
    det_GE11Long4->setBaseAddress("0011");
    det_GE11Long4->setChannelMapVME2Det(16, 1);
    det_GE11Long4->setChannelMapVME2Det(17, 2);
    det_GE11Long4->setChannelMapVME2Det(18, 3);
    det_GE11Long4->setChannelMapVME2Det(19, 4);
    det_GE11Long4->setChannelMapVME2Det(20, 5);
    det_GE11Long4->setChannelMapVME2Det(21, 6);
    
    cosmicStand.setDetector(0,1, det_PMT1, kTrigger);
    cosmicStand.setDetector(0,2, det_PMT2, kTrigger);
    cosmicStand.setDetector(0,3, det_PMT3, kTrigger);
    cosmicStand.setDetector(0,4, det_PMT4, kTrigger);
    cosmicStand.setDetector(0,5, det_PMT5, kTrigger);
    
    cosmicStand.setDetector(3,2, det_GE11Long2, kDUT);
    cosmicStand.setDetector(3,2, det_GE11Long4, kDUT);
    
    cosmicStand.setDetector(7,1, det_PMT6, kTrigger);
    cosmicStand.setDetector(7,2, det_PMT8, kTrigger);
    cosmicStand.setDetector(7,3, det_PMT9, kTrigger);
    cosmicStand.setDetector(7,4, det_PMT11, kTrigger);
    cosmicStand.setDetector(7,5, det_PMT12, kTrigger);
    
    auto standMatrix = cosmicStand.getMatrix();
    for (auto iterSuperChamber = standMatrix.m_map_detectors.begin(); iterSuperChamber != standMatrix.m_map_detectors.end(); ++iterSuperChamber) { //Loop Over Superchambers
        for (auto iterDet = (*iterSuperChamber).second->getDetectorPtrBegin(); iterDet != (*iterSuperChamber).second->getDetectorPtrEnd(); ++iterDet) {
            cout<<((*iterSuperChamber).first).first<<"\t";
            cout<<((*iterSuperChamber).first).second<<"\t";
            cout<<std::distance((*iterSuperChamber).second->getDetectorPtrBegin(), iterDet)<<"\t";
            cout<<(*iterDet).second->getName()<<endl;
        }
    }
    
    /*EventDigi evtDigi;
    
    TDCDataDigi digiPMT, digiDet;
    
    digiPMT.m_uiEvtCount = 1;
    digiPMT.m_strBaseAddress = "400C";
    digiPMT.m_map_fTime[1]=13.3;
    digiPMT.m_map_fTime[7]=17.0;
    
    digiDet.m_uiEvtCount = 1;
    digiDet.m_strBaseAddress = "0011";
    digiDet.m_map_fTime[11]=15.1;
    digiDet.m_map_fTime[13]=15.2;
    digiDet.m_map_fTime[15]=15.3;
    digiDet.m_map_fTime[17]=15.4;
    
    evtDigi.m_uiEvtCount = 1;
    evtDigi.m_map_TDCData[digiPMT.m_strBaseAddress]=digiPMT;
    evtDigi.m_map_TDCData[digiDet.m_strBaseAddress]=digiDet;
    
    EventReco evtReco = cosmicStand.getEventReco(evtDigi);
    
    printEvent(evtReco);*/
    
    HwVMEBoard vmeBoard_A;
    vmeBoard_A.m_strBaseAddress = "0x400C0000";
    vmeBoard_A.m_vme_type = kVMETDC;
    
    HwVMEBoard vmeBoard_B;
    vmeBoard_B.m_strBaseAddress = "0x00110000";
    vmeBoard_B.m_vme_type = kVMETDC;
    
    RunSetup daqSetup;
    daqSetup.m_map_vmeBoards[vmeBoard_A.m_strBaseAddress]=vmeBoard_A;
    daqSetup.m_map_vmeBoards[vmeBoard_B.m_strBaseAddress]=vmeBoard_B;
    
    string strInputFileName = "/Users/dorney/Desktop/MyResearch/GitRepos/CMS_GEM_Analysis_Framework/data/cosmicTesting/GE11-VII-L-CERN-0004_Run010_Timing_700uA_SBITModeiPhi_100kEvt.root";
    SelectorTiming selector(daqSetup);
    
    cout<<"main() - No Fault 1\n";
    
    vector<EventReco> vec_recoEvts = selector.getEventsReco(strInputFileName, cosmicStand);
    
    cout<<"main() - No Fault 2\n";
    
    cout<<"vec_recoEvts.size() = " << vec_recoEvts.size() << endl;
    
    for (int i=0; i < vec_recoEvts.size(); ++i) {
        if ( (i % 1000) == 0) {
            printEvent(vec_recoEvts[i]);
        }
    }
    
    return 1;
} //End main()