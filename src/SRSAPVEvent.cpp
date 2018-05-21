#include <cmath>
#include "SRSAPVEvent.h"


//=====================================================
SRSAPVEvent::SRSAPVEvent(int fec_no, int fec_channel, int apv_id, int zeroSupCut, int packetSize) {
  this->Clear() ;
  SRSMapping * mapping = SRSMapping::GetInstance();
  
  fPlane        = "GEM1X" ;
  fDetector     = "GEM1" ;
  fReadoutBoard = "CARTESIAN" ;
  fDetectorType = "STANDARD" ;
  
  fPlaneSize       = 102.4 ;
  fAPVIndexOnPlane = 2;
  fAPVOrientation  = 0 ;
  fAPVstripmapping = 0 ;
  fAPVstripmapping = 0 ;
  fAPVHeaderLevel  = 1300 ;
  fEtaSectorPos = 0 ;
  
  fTrapezoidDetOuterRadius = 430.0 ;
  fTrapezoidDetInnerRadius = 220.0 ;
  
  fPedSubFlag           = false ;
  fCommonModeFlag       = false ;
  fIsCosmicRunFlag      = false ;
  fIsPedestalRunFlag    = false ;
  fIsRawPedestalRunFlag = false ;
  
  fFECNo = fec_no ;
  fAPVID = apv_id ;
  
  fADCChannel = fec_channel ;
  fPacketSize = packetSize;
  
  fAPVGain = 1.0 ;
  
  fZeroSupCut = zeroSupCut ;
  fMeanAPVnoise = 0;
  
  fAPVKey          = mapping->GetAPVNoFromID(apv_id);
  fAPV             = mapping->GetAPVFromID(apv_id);
  fAPVHeaderLevel  = mapping->GetAPVHeaderLevelFromID(apv_id);
  fAPVIndexOnPlane = mapping->GetAPVIndexOnPlane(apv_id);
  fAPVOrientation  = mapping->GetAPVOrientation(apv_id);
  fAPVstripmapping = mapping->GetAPVstripmapping(apv_id);

    
  fPlane        = mapping->GetDetectorPlaneFromAPVID(apv_id);
  fDetector     = mapping->GetDetectorFromPlane(fPlane) ;
  fReadoutBoard = mapping->GetReadoutBoardFromDetector(fDetector) ;
  fDetectorType = mapping->GetDetectorTypeFromDetector(fDetector) ;
  
    
  if(fReadoutBoard == "CARTESIAN") {
    fPlaneSize         = (mapping->GetCartesianReadoutMap(fPlane))[1];
    fNbOfAPVsFromPlane = (int) ((mapping->GetCartesianReadoutMap(fPlane)) [2]) ;
    fPadDetectorMap.resize(5) ;
        //    printf("SRSAPVEvent::SRSAPVEvent  detType=%s, plane=%s, planeSize=%f, fNbOfAPVsFromPlane=%d\n", fDetectorType.Data(), fPlane.Data(), fPlaneSize, fNbOfAPVsFromPlane) ;
    }
    
    if(fReadoutBoard == "UV_ANGLE") {
        fTrapezoidDetLength      = (mapping->GetUVangleReadoutMap(fDetector)) [0] ;
        fTrapezoidDetInnerRadius = (mapping->GetUVangleReadoutMap(fDetector)) [1] ;
        fTrapezoidDetOuterRadius = (mapping->GetUVangleReadoutMap(fDetector)) [2] ;
        fNbOfAPVsFromPlane = (int) ((mapping->GetUVangleReadoutMap(fPlane)) [1]) ;
        fPadDetectorMap.resize(5) ;
        //    printf("SRSAPVEvent::SRSAPVEvent => fDetector =%s,  fPlane=%s, detLength=%f, innerRadius=%f, outerRadius=%f, fNbOfAPVsFromPlane=%d \n",fDetector.Data(), fPlane.Data(), fTrapezoidDetLength, fTrapezoidDetInnerRadius, fTrapezoidDetOuterRadius, fNbOfAPVsFromPlane);
    }
    
    if(fReadoutBoard == "1DSTRIPS") {
        fPlaneSize = (mapping->Get1DStripsReadoutMap(fPlane)) [1] ;
        fNbOfAPVsFromPlane = (int) ((mapping->Get1DStripsReadoutMap(fPlane)) [2]) ;
        fPadDetectorMap.resize(5) ;
        //    printf("SRSAPVEvent::SRSAPVEvent  fPlane=%s, planeSize=%f, fNbOfAPVsFromPlane=%d \n", fPlane.Data(), fPlaneSize, fNbOfAPVsFromPlane);
    }
    
    if(fReadoutBoard == "CMSGEM") {
        fEtaSectorPos = (mapping->GetCMSGEMReadoutMap(fPlane)) [0] ;
        fPlaneSize    = (mapping->GetCMSGEMReadoutMap(fPlane)) [1] ;
        fNbOfAPVsFromPlane = (int) ((mapping->GetCMSGEMReadoutMap(fPlane)) [2]) ;
        fPadDetectorMap.resize(5) ;
	//printf("SRSAPVEvent::SRSAPVEvent   fPlane=%s, etaSectorPos=%f, planeSize=%f, fNbOfAPVsFromPlane=%d, \n", fPlane.c_str(), fEtaSectorPos, fPlaneSize, fNbOfAPVsFromPlane);
    }
    
    else if(fReadoutBoard == "PADPLANE") {
      std::vector< int> apvChPadCh = mapping->GetPadChannelsMapping(apv_id) ;
      std::vector< int>::const_iterator apvChPadCh_itr ;
      for (apvChPadCh_itr = apvChPadCh.begin(); apvChPadCh_itr != apvChPadCh.end(); ++apvChPadCh_itr) {
	int apvChPadCh = * apvChPadCh_itr ;
	int padCh =  (apvChPadCh >> 8) & 0xffff ;
	int apvCh = apvChPadCh & 0xff ;
	fapvChToPadChMap[apvCh] = padCh ;
	fPadDetectorMap = mapping->GetPadDetectorMap(fDetector) ;
	fPadDetectorMap.resize(5) ;
      }
    }
    //  printf("SRSAPVEvent::SRSAPVEvent   fPlane=%s, etaSectorPos=%f, planeSize=%f, fNbOfAPVsFromPlane=%d, \n", fPlane.Data(), fEtaSectorPos, fPlaneSize, fNbOfAPVsFromPlane) ;
}

//=====================================================
void SRSAPVEvent::Clear() {
    fRawData32bits.clear() ;
    fRawData16bits.clear() ;
    fapvTimeBinDataMap.clear() ;
    fPedestalData.clear() ;
    fRawPedestalData.clear() ;
    fMaskedChannels.clear() ;
    fPedestalNoises.clear() ;
    fRawPedestalNoises.clear() ;
    
    fCommonModeOffsets.clear() ;
    fCommonModeOffsets_odd.clear() ;
    fCommonModeOffsets_even.clear() ;
    fPadDetectorMap.clear() ;
    //  printf(" === SRSAPVEvent::Clear()\n") ;
}

//=====================================================
SRSAPVEvent::~SRSAPVEvent(){
    Clear() ;
}

//=====================================================

void SRSAPVEvent::Add32BitsRawData(unsigned int rawData32bits) {
    //  printf("SRSAPVEvent::Add32BitsRawData()\n");
    fRawData32bits.push_back(rawData32bits) ;
}

//=====================================================
void SRSAPVEvent::Set32BitsRawData(std::vector<unsigned int> rawData32bits) {
    //  printf("SRSAPVEvent::Set32BitsRawData()\n");
    fRawData32bits.clear() ;
    fRawData32bits = rawData32bits ;
}

//=====================================================
int  SRSAPVEvent::APVchannelCorrection(int chNo) {
  chNo = (32 * (chNo%4)) + (8 * (int)(chNo/4)) - (31 * (int)(chNo/16)) ;
  return chNo ;

}


//=====================================================
int SRSAPVEvent::NS2StripMapping(int chNo) {
    if((chNo%2)==1){
        chNo=((chNo-1)/2)+64;
    }
    else{
        chNo=(chNo/2);
    }
    return chNo ;
}

//=====================================================
int SRSAPVEvent::CMSStripMapping(int chNo) {
  if((chNo%2)==1){ //Case: Input is Odd
    //return (128 - ( (chNo - 1) / 2 ) );//ODD
    return ( ( chNo - 1) / 2 );//ODD input va da 1 a 127 ->0 63
  } //End Case: Input is Odd
  else{ //Case: Input is Even
    //return ( ( chNo / 2 ) - 1 );//EVEN
    return ( 127 - ( chNo / 2 ) );//EVEN  126 a 0                  ->64...65.. 127
  } //End Case: Input is Even
}


//=====================================================
//For GE1/1-VII-L/S design
//Use for specific sectors from the set {(7,1),(7,3),(8,1),(8,3)}
int SRSAPVEvent::CMSStripMapping1(int chNo) {
  if((chNo%2)==1){ //Case: Input is Odd
    //return ( (chNo + 125) / 2 );//ODD
    return ( ( 127 - chNo ) / 2 );//ODD
  } //End Case: Input is Odd
  else{ //Case: Input is Even
    //return ( ( 128 - chNo) / 2 );//EVEN
    return (1+ ( chNo + 126 ) / 2 );//EVEN            127...64
  } //End Case: Input is Even
    
    //return chNo ;
} //End SRSAPVEvent::CMSStripMapping1()

//=====================================================
//For GE1/1-VII-L/S design
//Use for specific sectors from the set {(5,1),(6,1)}
int SRSAPVEvent::CMSStripMapping2(int chNo) {
  if((chNo%2)==1){ //Case: Input is Odd
    //return ( ( 127 - chNo ) / 2 );//ODD
    return (1+ ( chNo + 125 ) / 2 );//ODD
  } //End Case: Input is Odd
  else{ //Case: Input is Even
    //return ( ( chNo + 126 ) / 2 );//EVEN
    return ( ( 128 - chNo ) / 2 )-1;//EVEN
  } //End Case: Input is Even

    //return chNo ;
} //End SRSAPVEvent::CMSStripMapping2()

//=====================================================
//For GE1/1-VII-S design
//Use for specific sectors from ieta rows 2 & 4
int SRSAPVEvent::CMSStripMapping3(int chNo) {
  if((chNo%2)==1){ //Case: Input is Odd
    return ( ( 127 - chNo ) / 2 );//ODD
  } //End Case: Input is Odd
  else{ //Case: Input is Even
    return ( 1+( chNo + 126 ) / 2 );//EVEN
  } //End Case: Input is Even

    //return chNo ;
} //End SRSAPVEvent::CMSStripMapping3()

//=====================================================
int SRSAPVEvent::CMSStripMapping4(int chNo) {
  //chNo=4;
  return chNo ;
}


//=====================================================
int SRSAPVEvent::MMStripMappingAPV1(int chNo){
    if((chNo%2)==1){
        chNo=((chNo-1)/2) + 32;
    }
    else{
        chNo=(chNo/2);
        if (chNo < 32) chNo = 31 - chNo;
        else if (chNo > 37) chNo = 159 - chNo ;
        else chNo += 90;
    }
    return chNo;
}

//=====================================================
int SRSAPVEvent::MMStripMappingAPV2(int chNo){
    if((chNo%2)==1){
        chNo=((chNo-1)/2) + 27;
    }
    else{
        chNo=(chNo/2);
        if (chNo < 27) chNo = 26 - chNo;
        else if (chNo > 38) chNo = 154 - chNo;
        else chNo += 89;
    }
    return chNo;
}

//=====================================================
int SRSAPVEvent::MMStripMappingAPV3(int chNo){
    if((chNo%2)==1){
        chNo=((chNo-1)/2) + 26;
    }
    else{
        chNo=(chNo/2);
        if (chNo < 26) chNo = 25 - chNo;
        else if (chNo > 31) chNo = 153 - chNo;
        else chNo += 96;
    }
    return chNo;
}

//=====================================================
//int SRSAPVEvent::PRadStripsMapping(int chNo) {
//  int chno = chNo / 2 ;
//  chNo = chNo /2 ;
//  if (chNo % 2 == 0) chNo = 31 + (chNo / 2) ;

//  if (chNo % 2 == 0) chNo = 31 + (chNo / 2) ;
//  else {
//    if (chNo < 64 ) chNo = 31 -  ( (chNo + 1) / 2) ;
//    else            chNo = 127 - ( (chNo - 65) / 2 ) ;
//  }
//  printf("SRSAPVEvent::PRadStripsMapping ==>  APVID=%d, chNo=%d, stripNo=%d, \n",fAPVID, chno, chNo) ;

//  return chNo ;
//}

//=====================================================
int SRSAPVEvent::StandardMapping(int chNo) {
    return chNo ;
}

//=====================================================
int SRSAPVEvent::EICStripMapping(int chNo) {
    if(chNo % 2 == 0) chNo = chNo / 2 ;
    else              chNo = 64 + (chNo - 1) / 2 ;
    return chNo ;
}

//=====================================================
int SRSAPVEvent::HMSStripMapping(int chNo) {
    if(chNo % 4 == 0)      chNo = chNo + 2 ;
    else if(chNo % 4 == 1) chNo = chNo - 1 ;
    else if(chNo % 4 == 2) chNo = chNo + 1 ;
    else if(chNo % 4 == 3) chNo = chNo - 2 ;
    else chNo = chNo ;
    return chNo ;
}

//=====================================================
int SRSAPVEvent::StripMapping(int chNo) {

    chNo = this->APVchannelCorrection(chNo) ;
    //    std::cout <<" ++++ detector type "<<fDetectorType<<" APVstrip mapping "<<fAPVstripmapping<<std::endl;
    if (fDetectorType == "CMSGEM") {
        if (fAPVstripmapping==0) chNo = CMSStripMapping(chNo) ;
        if (fAPVstripmapping==1) chNo = CMSStripMapping1(chNo) ;
        if (fAPVstripmapping==2) chNo = CMSStripMapping2(chNo) ;
        if (fAPVstripmapping==3) chNo = CMSStripMapping3(chNo) ;
        if (fAPVstripmapping==4) chNo = CMSStripMapping4(chNo) ;
	//	std::cout << " ++++ strip = "<<chNo<<std::endl;
        //cout << "AA " << chNo << endl;
    }

    //    else if (fDetectorType == "ZIGZAG")    chNo = ZigZagStripMapping(chNo) ;
    else if (fDetectorType == "NS2")       chNo = NS2StripMapping(chNo) ;
    else if (fDetectorType == "EICPROTO1") chNo = EICStripMapping(chNo) ;
    else if (fDetectorType == "HMSGEM")    chNo = HMSStripMapping(chNo) ;
    //else if (fDetectorType == "PRADGEM")   chNo = PRadStripsMapping(chNo) ;
    else                                   chNo = StandardMapping(chNo) ;
    //cout << "BB " << chNo << endl;

    return chNo;
}

//========================================================================================================================
static bool usingGreaterThan(float u, float v) {
    return u > v  ;
}

//=====================================================
void SRSAPVEvent::Print() {
    //  printf("SRSAPVEvent::Print() => Printing APV Data: -1 == No Data, -2 == Underflow, -3 == Overflow\n");
    int Size = fRawData32bits.size() ;
    int Capacity = fRawData32bits.capacity() ;
    int MaxSize =  fRawData32bits.max_size() ;
    std::vector<unsigned int>::const_iterator rawdata_itr ;
    for (rawdata_itr = fRawData32bits.begin(); rawdata_itr != fRawData32bits.end(); ++rawdata_itr) {
        int rawdata = * rawdata_itr ;
        printf("SRSAPVEvent::Print() ====> 32 bits raw data [0x%x]\n",rawdata) ;
    }
}

void SRSAPVEvent::ComputeRawData16bits() {
    
    //printf("SRSAPVEvent::ComputeRawDataZS()==> fRawData32bits.size() = %d \n",(int) (fRawData32bits.size()) ) ;
    fRawData16bits.clear() ;
    std::vector<unsigned int>::const_iterator rawData_itr ;
    for (rawData_itr = fRawData32bits.begin(); rawData_itr != fRawData32bits.end(); ++rawData_itr) {
        unsigned int word32bit = * rawData_itr ;
        //cout << "Word32bit is = " << word32bit << endl;
        
        if (((word32bit >> 8) & 0xffffff) != 0x414443) {
            unsigned int data1 = (word32bit>> 24) & 0xff ;
            unsigned int data2 = (word32bit >> 16) & 0xff ;
            unsigned int data3 = (word32bit >> 8)  & 0xff ;
            unsigned int data4 = word32bit  & 0xff ;
            fRawData16bits.push_back(((data2 << 8) | data1)) ;
            fRawData16bits.push_back(((data4 << 8) | data3)) ;
        }
    }
}


void SRSAPVEvent::ComputeRawData16bitsZS() {
    
    //printf("SRSAPVEvent::ComputeRawDataZS()==> fRawData32bits.size() = %d \n",(int) (fRawData32bits.size()) ) ;
    fRawData16bits.clear() ;
    std::vector<unsigned int>::const_iterator rawData_itr ;
    for (rawData_itr = fRawData32bits.begin(); rawData_itr != fRawData32bits.end(); ++rawData_itr) {
        unsigned int word32bit = * rawData_itr ;
        //cout << "Word32bit is = " << word32bit << endl;
        if (((word32bit >> 8) & 0xffffff) != 0x41505a) {
            unsigned int data1 = (word32bit>> 24) & 0xff ;
            unsigned int data2 = (word32bit >> 16) & 0xff ;
            unsigned int data3 = (word32bit >> 8)  & 0xff ;
            unsigned int data4 = word32bit  & 0xff ;
            fRawData16bits.push_back(((data2 << 8) | data1)) ;
            fRawData16bits.push_back(((data4 << 8) | data3)) ;
        }
    }
}

//=====================================================
void SRSAPVEvent::ComputeTimeBinCommonMode() {
    // printf("  SRSAPVEvent::ComputeTimeBinCommonMode()==>enter \n") ;
    
    int idata = 0  ;
    bool startDataFlag = false ;
    
    fapvTimeBinDataMap.clear() ;
    unsigned int apvheaderlevel = (unsigned int)  fAPVHeaderLevel ;
    int size = fRawData16bits.size() ;
    
    if (size != fPacketSize) {
        printf("SRSAPVEvent::ComputeTimeBinCommonMode() XXXX ERROR XXXXX ==> Packet size %d different from expected %d, header=%d \n",size,fPacketSize, fAPVHeaderLevel ) ;
    }
    
    std::list<float> commonModeOffset,  commonModeOffset_odd, commonModeOffset_even ;

    float vtot = 0;
    for (std::vector<float>::iterator i=fPedestalOffsets.begin();i<fPedestalOffsets.end();i++){
      vtot+=*i;
    }
    //    float apvBaseline = 4096 - TMath::Mean(fPedestalOffsets.begin(), fPedestalOffsets.end()) ;
    float apvBaseline = 4096 - vtot/float(fPedestalOffsets.size());
    
    while(idata < size) {
        //    printf("    SRSAPVEvent::ComputeTimeBinCommonMode()  ==> idata=%d, apvBaseline=%f \n",idata, apvBaseline) ;
        //===============================================================//
        // If 3 consecutive words satisfy this condition below => it is  //
        // an APV header so we could take meaninfull data                //
        //===============================================================//
        if (fRawData16bits[idata] < apvheaderlevel ) {
            idata++ ;
            if (fRawData16bits[idata] < apvheaderlevel ) {
                idata++ ;
                if (fRawData16bits[idata] < apvheaderlevel ) {
                    idata += 10;
                    startDataFlag = true ;
                    continue ;
                }
            }
        }
        
        //===============================================================//
        // That's where the meaninfull data are taken                    //
        // 128 analog word for each apv strip and each timebin           //
        //===============================================================//
        if (startDataFlag == true) {
            commonModeOffset.clear(), commonModeOffset_odd.clear(), commonModeOffset_even.clear() ;
            
            float commMode      = 0 ;
            float commMode_odd  = 0 ;
            float commMode_even = 0 ;
            
            for(int chNo = 0; chNo < NCH; ++chNo) {
                
                int stripNo = StripMapping(chNo) ;
                
                float rawdata = ((float) fRawData16bits[idata]) ;
                float comMode = rawdata ;
                
                float thresohld = 500 ;
                if ((fReadoutBoard == "UV_ANGLE") &&  (fDetectorType == "EICPROTO1") ) thresohld = 500 	;
                
                if (fabs(comMode - apvBaseline) > thresohld) {
                    comMode = apvBaseline ;
                }
                
                rawdata = 4096 - rawdata ;
                comMode = 4096 - comMode ;
                
                if (fIsCosmicRunFlag) {
                    if ((fReadoutBoard == "UV_ANGLE") &&  (fDetectorType == "EICPROTO1") ) {
                        if(stripNo < 64) commonModeOffset_even.push_back(comMode - fPedestalOffsets[stripNo]) ;
                        else             commonModeOffset_odd.push_back(comMode - fPedestalOffsets[stripNo]) ;
                    }
                    else {
                        commonModeOffset.push_back(comMode - fPedestalOffsets[stripNo]) ;
                    }
                }
                
                if (fIsPedestalRunFlag) {
                    if ((fReadoutBoard == "UV_ANGLE") &&  (fDetectorType == "EICPROTO1") ) {
                        if(stripNo < 64) commonModeOffset_even.push_back(comMode - fRawPedestalOffsets[stripNo]) ;
                        else             commonModeOffset_odd.push_back(comMode - fRawPedestalOffsets[stripNo]) ;
                    }
                    else {
                        commonModeOffset.push_back(comMode - fRawPedestalOffsets[stripNo]) ;
                    }
                }
                fapvTimeBinDataMap.insert(std::pair<int, float>(stripNo, rawdata)) ;
                idata++ ;
            }
            
            if ((fReadoutBoard == "UV_ANGLE") &&  (fDetectorType == "EICPROTO1") ) {
                commonModeOffset_odd.sort() ;
                commonModeOffset_even.sort() ;
                
		float ovtot = 0;
		for (std::list<float>::iterator i=commonModeOffset_odd.begin();i!=commonModeOffset_odd.end();i++){
		  ovtot+=*i;
		}
                commMode_odd  = ovtot/float(commonModeOffset_odd.size());
		  //TMath::Mean(commonModeOffset_odd.begin(),commonModeOffset_odd.end()) ;
		float evtot = 0;
		for (std::list<float>::iterator i=commonModeOffset_even.begin();i!=commonModeOffset_even.end();i++){
		  evtot+=*i;
		}
                commMode_even  = evtot/float(commonModeOffset_even.size());
		//                commMode_even = TMath::Mean(commonModeOffset_even.begin(),commonModeOffset_even.end()) ;
                commonModeOffset_odd.clear() ;
                commonModeOffset_even.clear() ;
            }
            else {
                commonModeOffset.sort() ;
		float vtot = 0;
		for (std::list<float>::iterator i=commonModeOffset.begin();i!=commonModeOffset.end();i++){
		  vtot+=*i;
		}
                commMode = vtot/float(commonModeOffset.size());
		//                commMode = TMath::Mean(commonModeOffset.begin(),commonModeOffset.end()) ;
                commonModeOffset.clear() ;
            }
            
            if (fIsRawPedestalRunFlag) {
                if ((fReadoutBoard == "UV_ANGLE") &&  (fDetectorType == "EICPROTO1") ) {
                    fCommonModeOffsets_odd.push_back(0) ;
                    fCommonModeOffsets_even.push_back(0) ;
                }
                else {
                    fCommonModeOffsets.push_back(0) ;
                }
            }
            else {
                if ((fReadoutBoard == "UV_ANGLE") &&  (fDetectorType == "EICPROTO1") ) {
                    fCommonModeOffsets_odd.push_back(commMode_odd) ;
                    fCommonModeOffsets_even.push_back(commMode_even) ;
                }
                else {
                    fCommonModeOffsets.push_back(commMode) ;
                }
            }
            startDataFlag = false ;
            continue ;
        }
        idata++ ;
    }
    //  printf("  SRSAPVEvent::ComputeTimeBinCommonMode()==>exit \n") ;
}

//===================================================================
std::list <SRSHit * >  SRSAPVEvent::ComputeListOfAPVHits() {
    //  printf("SRSAPVEvent::ComputeListOfAPVHits()==> enter \n") ;
    
    fIsCosmicRunFlag = true ;
    fIsPedestalRunFlag = false ;
    fIsRawPedestalRunFlag = false ;
    
    ComputeRawData16bits() ;
    ComputeTimeBinCommonMode() ;
    std::list <SRSHit * > listOfHits ;
    
    int apvTimeBinDataMapSize = fapvTimeBinDataMap.size() ;
    int padNo = 0 ;
    std::string plane = fPlane ;
    
    if (apvTimeBinDataMapSize != 0) {
        std::vector <float> stripPedestalNoise ;
        
	std::pair<std::multimap<int, float>::iterator, std::multimap<int, float>::iterator> stripSetOfTimeBinRawData ;
        for (int stripNo = 0;  stripNo < NCH; stripNo++) {
            //     printf("SRSAPVEvent::ComputeListOfAPVHits()==> stripNo=%d \n", stripNo) ;
            
	  std::vector<float> timeBinADCs;
            int timebin = 0 ;
            
            stripSetOfTimeBinRawData = fapvTimeBinDataMap.equal_range(stripNo) ;
            
	    std::multimap <int, float>::iterator timebin_it ;
            for (timebin_it = stripSetOfTimeBinRawData.first; timebin_it != stripSetOfTimeBinRawData.second; ++timebin_it) {
                float rawdata =  timebin_it->second ;
                //	printf("SRSAPVEvent::ComputeListOfAPVHits()==> stripNo=%d, rawdata=%f \n", stripNo, rawdata) ;
                
                // BASELINE CORRECTION
                if (fCommonModeFlag) {
                    if ((fReadoutBoard == "UV_ANGLE") &&  (fDetectorType == "EICPROTO1") ) {
                        if(stripNo < 64)  rawdata -= fCommonModeOffsets_even[timebin] ;
                        else              rawdata -= fCommonModeOffsets_odd[timebin] ;
                    }
                    else  {
                        rawdata -= fCommonModeOffsets[timebin] ;
                    }
                }
                
                // Pedestal Offset Suppression
                if(fPedSubFlag) rawdata -= fPedestalOffsets[stripNo] ;
                
                // Masked Channels
                rawdata = (1 - fMaskedChannels[stripNo]) * rawdata ;
                
                // APV GAIN CORRECTION DEFAULT VALUE GAIN = 1
                rawdata =  rawdata / fAPVGain ;
                timeBinADCs.push_back(rawdata) ;
                timebin++ ;
            }
            
            if(fReadoutBoard == "PADPLANE") {
                padNo = fapvChToPadChMap[stripNo] ;
                if (padNo == 65535) continue ;
            }
            
            int stripNb = stripNo ;
            if ((fReadoutBoard == "UV_ANGLE") &&  (fDetectorType == "EICPROTO1") ) {
                if(stripNo > 63) {
                    stripNb = stripNo - 64 ;
                    plane = fDetector + "BOT" ;
                }
                else    {
                    plane = fDetector + "TOP" ;
                }
            }
            
            // ZERO SUPPRESSION
            if (fZeroSupCut > 0 ) {
	      float tval= 0;
	      float tmax= 0;
	      for (std::vector<float>::iterator i=timeBinADCs.begin();i<timeBinADCs.end();i++){
		  tval+=*i;
		  if (tmax < *i) tmax = *i;
	      }
	      tval = tval / float(timeBinADCs.size());
	      //	      if( TMath::Mean(timeBinADCs.begin(), timeBinADCs.end()) < fZeroSupCut * fPedestalNoises[stripNo] ) {
	      if( tval < fZeroSupCut * fPedestalNoises[stripNo] ) {
		stripPedestalNoise.push_back(tval) ;
		timeBinADCs.clear() ;
		continue ;
	      }
	      
	      if ( (int) (timeBinADCs.size()) != 0 ) {
		//		float adcs = * (TMath::LocMax(timeBinADCs.begin(), timeBinADCs.end())) ;
		float adcs = tmax;
		if(adcs < 0) adcs = 0 ;
		SRSHit * apvHit = new SRSHit() ;
		apvHit->SetAPVID(fAPVID) ;
		apvHit->IsHitFlag(true) ;
		apvHit->SetDetector(fDetector) ;
		apvHit->SetDetectorType(fDetectorType) ;
		apvHit->SetReadoutBoard(fReadoutBoard) ;
		apvHit->SetPadDetectorMap(fPadDetectorMap) ;
		apvHit->SetPlane(plane) ;
		apvHit->SetPlaneSize(fPlaneSize) ;
		apvHit->SetTrapezoidDetRadius(fTrapezoidDetInnerRadius, fTrapezoidDetOuterRadius) ;
		apvHit->SetAPVOrientation(fAPVOrientation) ;
		apvHit->SetAPVstripmapping(fAPVstripmapping) ;
		apvHit->SetAPVIndexOnPlane(fAPVIndexOnPlane) ;
		apvHit->SetNbAPVsFromPlane(fNbOfAPVsFromPlane) ;
		apvHit->SetTimeBinADCs(timeBinADCs) ;
		apvHit->SetHitADCs(fZeroSupCut, adcs, fIsHitMaxOrTotalADCs) ;
		apvHit->SetPadNo(padNo) ;
		apvHit->SetStripNo(stripNb) ;
		listOfHits.push_back(apvHit) ;
	      }
            }
            
            else {
	      SRSHit * apvHit = new SRSHit() ;
	      apvHit->SetAPVID(fAPVID) ;
	      
	      float tval= 0;
	      float tmax= 0;
	      for (std::vector<float>::iterator i=timeBinADCs.begin();i<timeBinADCs.end();i++){
		tval+=*i;
		if (tmax < *i) tmax = *i;
	      }
	      tval = tval / float(timeBinADCs.size());
	      float adcs = tmax;
	      //float adcs = * (TMath::LocMax(timeBinADCs.begin(), timeBinADCs.end())) ;
	      if ( tval < 5 * fPedestalNoises[stripNo] ) {
		//                if( TMath::Mean(timeBinADCs.begin(), timeBinADCs.end()) < 5 * fPedestalNoises[stripNo] ) {
		apvHit->IsHitFlag(false) ;
		//		adcs = TMath::Mean(timeBinADCs.begin(), timeBinADCs.end());		
		adcs = tval;
	      }
	      
	      apvHit->SetDetector(fDetector) ;
	      apvHit->SetDetectorType(fDetectorType) ;
	      apvHit->SetReadoutBoard(fReadoutBoard) ;
	      apvHit->SetPadDetectorMap(fPadDetectorMap) ;
	      apvHit->SetPlane(plane) ;
	      apvHit->SetPlaneSize(fPlaneSize) ;
	      apvHit->SetTrapezoidDetRadius(fTrapezoidDetInnerRadius, fTrapezoidDetOuterRadius) ;
	      apvHit->SetAPVOrientation(fAPVOrientation) ;
	      apvHit->SetAPVstripmapping(fAPVstripmapping) ;
	      apvHit->SetAPVIndexOnPlane(fAPVIndexOnPlane) ;
	      apvHit->SetNbAPVsFromPlane(fNbOfAPVsFromPlane) ;
	      apvHit->SetTimeBinADCs(timeBinADCs) ;
	      apvHit->SetHitADCs(fZeroSupCut, adcs, fIsHitMaxOrTotalADCs) ;
	      apvHit->SetPadNo(padNo) ;
	      apvHit->SetStripNo(stripNb) ;
	      listOfHits.push_back(apvHit) ;
            }
            timeBinADCs.clear() ;
        }
	float nav1 = 0;
	float nav2 = 0;
	for (std::vector<float>::iterator i=stripPedestalNoise.begin();i<stripPedestalNoise.end();i++){
	  nav1+=*i;
	  nav2+=(*i)*(*i);
	}
	float rms2 = (nav2-nav1*nav1/float(stripPedestalNoise.size()))/(float(stripPedestalNoise.size()-1));
	float rmsnoise = sqrt(rms2);
        fMeanAPVnoise  = rmsnoise;
	//        fMeanAPVnoise  = TMath::RMS(stripPedestalNoise.begin(), stripPedestalNoise.end()) ;
        stripPedestalNoise.clear() ;
    }
    //  printf("SRSAPVEvent::ComputeListOfAPVHits()==> exit \n") ;
    
    return listOfHits ;
}

std::list <SRSHit * >  SRSAPVEvent::ComputeListOfAPVHitsZS() {
    //printf("SRSAPVEvent::ComputeListOfAPVHits()==> enter \n") ;
  
    fIsCosmicRunFlag = true ;
    fIsPedestalRunFlag = false ;
    fIsRawPedestalRunFlag = false ;
    
    this->ComputeRawData16bitsZS() ;
    
    int idata = 0  ;
    std::list <SRSHit * > listOfHits ;
    std::string plane = fPlane ;
    int padNo = 0 ;
    std::vector<float> timeBinADCs;
    timeBinADCs.clear();
    
    bool startDataFlag = false ;
    
    fapvTimeBinDataMap.clear() ;
    unsigned int apvheaderlevel = (unsigned int)  fAPVHeaderLevel ;
    int size = fRawData16bits.size() ;
    
    std::list<float> commonModeOffset,  commonModeOffset_odd, commonModeOffset_even ;
    //float apvBaseline = 4096 - TMath::Mean(fPedestalOffsets.begin(), fPedestalOffsets.end()) ;
    
    //    std::cout << "fRawData16bits.size = " << size << std::endl;
    while(idata < (size-1)) {
        //printf("    SRSAPVEvent::ComputeTimeBinCommonMode()  ==> idata=%d, apvBaseline=%f \n",idata, apvBaseline) ;
        
        //===============================================================//
        // Stefano Colafranceschi ZS code                                //
        //===============================================================//
        
        int adcbin, chNo;
        
        idata+=1;
        // Global Event Dump
        /*
	  cout << "GLOBAL DUMP fRawData16bits.size = " << size << " rawdata = " << fRawData16bits[idata] << " idata = " << idata << " hex is : " << "0x" << std::hex << fRawData16bits[idata] << std::dec << endl;
         cin.get();
         */
        
        // Reading out number of time bins from the header
        std::stringstream ss;
        ss << std::hex << (fRawData16bits[1]>>8);
        ss >> adcbin;
	//	std::cout <<" ++++ NUmber of adc bin"<<adcbin<<std::endl;

        // VALID DATA AFTER idata=4
        if ( idata>=4 ) {
            if ( size>4 ) {
	      //	      std::cout <<" ++++ Starting flag "<<std::endl;
	      startDataFlag = true ;
            }
        }
        
        if (startDataFlag == true) {
	  if ( ((idata-4)%(adcbin+1))==0 ) {
	    chNo = fRawData16bits[idata]; //--- Eraldo
	    //	    std::cout <<" ++++ first long word where channel number is stored "<<chNo<<" idata "<<idata<<std::endl;
	  }
            
	  if ( ((idata-4)%(adcbin+1))==adcbin ) {
	    int stripNo = this->StripMapping(chNo) ;
	    //	    std::cout <<" +++ last long word creating adc time profile for strip "<<stripNo<<std::endl;
	    if (timeBinADCs.size()>0) {
	      float tval= 0;
	      float tmax= 0;
	      for (std::vector<float>::iterator i=timeBinADCs.begin();i<timeBinADCs.end();i++){
		tval+=*i;
		if (tmax < *i) tmax = *i;
	      }
	      tval = tval / float(timeBinADCs.size());
	      float adcs = tmax;
	      //		  float adcs = * (TMath::LocMax(timeBinADCs.begin(), timeBinADCs.end())) ;
	      if(adcs < 0) adcs = 0 ;
	      SRSHit * apvHit = new SRSHit() ;
	      apvHit->SetAPVID(fAPVID) ;
	      apvHit->IsHitFlag(true) ;
	      apvHit->SetDetector(fDetector) ;
	      apvHit->SetDetectorType(fDetectorType) ;
	      apvHit->SetReadoutBoard(fReadoutBoard) ;
	      apvHit->SetPadDetectorMap(fPadDetectorMap) ;
	      apvHit->SetPlane(plane) ;
	      apvHit->SetPlaneSize(fPlaneSize) ;
	      apvHit->SetTrapezoidDetRadius(fTrapezoidDetInnerRadius, fTrapezoidDetOuterRadius) ;
	      apvHit->SetAPVOrientation(fAPVOrientation) ;
	      apvHit->SetAPVstripmapping(fAPVstripmapping) ;
	      apvHit->SetAPVIndexOnPlane(fAPVIndexOnPlane) ;
	      apvHit->SetNbAPVsFromPlane(fNbOfAPVsFromPlane) ;
	      apvHit->SetTimeBinADCs(timeBinADCs) ;
	      apvHit->SetHitADCs(fZeroSupCut, adcs, fIsHitMaxOrTotalADCs) ;
	      apvHit->SetPadNo(padNo) ;
	      apvHit->SetStripNo(stripNo) ;
	      //cout << "Filling..." << stripNo << endl;
	      listOfHits.push_back(apvHit) ;
	    }
	    //	    std::cout <<" ++++ clearing the adc time profile"<<std::endl;
	    timeBinADCs.clear();
	    //cout << "timeBinADCs cleared" << endl;
	  }
	  
	  // Getting the rawdata (Taking into account signed integers) and save them
	  else{
	    float rawdata = ((float) fRawData16bits[idata]) ;
	    if ((fRawData16bits[idata]>>8)!=0) rawdata = rawdata - 65536;
	    rawdata = - rawdata ;
	    //	  std::cout <<" ++++ normal time long word "<<rawdata<<std::endl;
	    timeBinADCs.push_back(rawdata) ;
	  }
          
            // Event Dump
            //cout << "fRawData16bits.size = " << size << " channel number = " << StripNo << " rawdata = " << fRawData16bits[idata] << " idata = " << idata << " hex is : " << "0x" << std::hex << fRawData16bits[idata] << std::dec << " MSB = " << (fRawData16bits[idata]>>8) << endl;
            
            //cin.get();
            
	  startDataFlag = false ;
	  continue ;
	}
    }
    //  printf("  SRSAPVEvent::ComputeTimeBinCommonMode()==>exit \n") ;
    return listOfHits ;
}


//========================================================================================================================
void SRSAPVEvent::ComputeMeanTimeBinRawPedestalData() { 
    //  printf("SRSAPVEvent::ComputeMeanTimeBinRawPedestalData() \n") ;
    
    fPedSubFlag = false ;
    fIsCosmicRunFlag = false ;
    
    fIsPedestalRunFlag = false ;
    fIsRawPedestalRunFlag = true ;
    
    ComputeRawData16bits() ;
    ComputeTimeBinCommonMode() ;
    fRawPedestalData.clear() ;
    
    std::vector<float>  meanTimeBinRawPedestalDataVect ;
    std::pair<std::multimap<int, float>::iterator, std::multimap<int, float>::iterator> stripSetOfTimeBinRawData ;
    for (int stripNo = 0;  stripNo < NCH; stripNo++) {
        stripSetOfTimeBinRawData = fapvTimeBinDataMap.equal_range(stripNo) ;
	std::multimap<int, float>::iterator timebin_it ;
        int timebin = 0 ;
        for (timebin_it = stripSetOfTimeBinRawData.first; timebin_it != stripSetOfTimeBinRawData.second; ++timebin_it) {
            float rawdata = timebin_it->second ;
            
            if ((fReadoutBoard == "UV_ANGLE") &&  (fDetectorType == "EICPROTO1") ) {
                if(stripNo < 64) meanTimeBinRawPedestalDataVect.push_back(rawdata - fCommonModeOffsets_even[timebin]) ;
                else             meanTimeBinRawPedestalDataVect.push_back(rawdata - fCommonModeOffsets_odd[timebin]) ; 
            }
            else {
                meanTimeBinRawPedestalDataVect.push_back(rawdata - fCommonModeOffsets[timebin]) ;
            }
            timebin++ ;
        }
	float avte = 0;
	for (std::vector<float>::iterator i=meanTimeBinRawPedestalDataVect.begin();i<meanTimeBinRawPedestalDataVect.end();i++){
	  avte += *i;
	}
	avte = avte/ float(meanTimeBinRawPedestalDataVect.size());
	//        fRawPedestalData.push_back(TMath::Mean(meanTimeBinRawPedestalDataVect.begin(),meanTimeBinRawPedestalDataVect.end())) ;
        fRawPedestalData.push_back(avte) ;
        meanTimeBinRawPedestalDataVect.clear() ;
    }
}

//========================================================================================================================
void SRSAPVEvent::ComputeMeanTimeBinPedestalData() {
    //  printf("SRSAPVEvent::ComputeMeanTimeBinPedestalData() \n") ;
    
    fPedSubFlag        = false ;
    fIsCosmicRunFlag   = false ;
    
    fIsPedestalRunFlag = true ;
    fIsRawPedestalRunFlag = false ;
    
    ComputeRawData16bits() ;
    ComputeTimeBinCommonMode() ;
    fPedestalData.clear() ;
    
    std::vector<float>  meanTimeBinPedestalDataVect ;
    std::pair<std::multimap<int, float>::iterator, std::multimap<int, float>::iterator> stripSetOfTimeBinRawData ;
    for (int stripNo = 0;  stripNo < NCH; stripNo++) {
        stripSetOfTimeBinRawData = fapvTimeBinDataMap.equal_range(stripNo) ;
	std::multimap<int, float>::iterator timebin_it ;
        int timebin = 0 ;
        for (timebin_it = stripSetOfTimeBinRawData.first; timebin_it != stripSetOfTimeBinRawData.second; ++timebin_it) {
            float rawdata = timebin_it->second ;
            
            if ((fReadoutBoard == "UV_ANGLE") &&  (fDetectorType == "EICPROTO1") ) {
                if(stripNo < 64) meanTimeBinPedestalDataVect.push_back(rawdata - fCommonModeOffsets_even[timebin]) ;
                else             meanTimeBinPedestalDataVect.push_back(rawdata - fCommonModeOffsets_odd[timebin]) ;
            }
            else {
                meanTimeBinPedestalDataVect.push_back(rawdata - fCommonModeOffsets[timebin]) ;
            }
            timebin++ ;
        }
	float avte = 0;
	for (std::vector<float>::iterator i=meanTimeBinPedestalDataVect.begin();i<meanTimeBinPedestalDataVect.end();i++){
	  avte += *i;
	}
	avte = avte/ float(meanTimeBinPedestalDataVect.size());
	
	//        fPedestalData.push_back(TMath::Mean(meanTimeBinPedestalDataVect.begin(),meanTimeBinPedestalDataVect.end())) ;
        fPedestalData.push_back(avte) ;
        meanTimeBinPedestalDataVect.clear() ;
    }
}
