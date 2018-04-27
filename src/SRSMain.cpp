//Framework Includes
#include "SRSMain.h"
#include "SRSConfiguration.h"
#include "SRSMapping.h"
#include "SRSHit.h"
#include "SRSAPVEvent.h"
#include "SRSFECDecoder.h"
#include "SRSEventBuilder.h"
#include "SRSOutputROOT.h"

//SRSMain* SRSMain::_repro = 0;

SRSMain::~SRSMain(){
  if (!_isClosed)
    this->Close();
};
 
/*SRSMain* 
SRSMain::Reprocessor(const std::string& rawfile, const std::string& config){
  if (_repro == 0) _repro = new SRSMain(rawfile,config);
  return _repro;
}*/


void 
SRSMain::Reprocess(){

  std::ifstream ifs(_rawfile.c_str(),std::ios::binary);
  unsigned int nMax = 5000;
  unsigned long buffer[nMax];
  std::vector<unsigned int> previous;
  unsigned long event[nMax];
  unsigned int nw = 0;
  bool newevent=true;
  unsigned int nEvent = 0;
  unsigned int prenevt =0;
  unsigned int iFEC=0,nFEC = 2;
  SRSEventBuilder* eventBuilder = 0;
  bool debug_s=false;
  do {
    if (newevent){
      newevent = false;
      nEvent++;
      if (nEvent%1000==0) std::cout<<" +++ reprocess: Event # "<<nEvent<<std::endl;
    }
    unsigned int j  =0x0;
    ifs.read((char*)&j,4);
    //    if (debug_s) std::cout <<" j  "<<std::hex<<j<<std::endl;
    if (!ifs.eof()){
      buffer[nw]=j&0xffffffff;
      nw++;
      if (nw >= nMax) {
	std::cout<<" +++ reprocess: Number of long words on this event exceed the maximum "<<nMax<<" please change nMax"<<std::endl;
	this->Close();
        exit(-1);
      }
      // Each FEC data finish with 0xfafafafa and one events is componsed by NFEC fragments. The event builder collect the hit of the                                                                    
      // fragments and make the clusters.                                                                                                                                                                
      if ((j&0xffffffff) == 0xfafafafa){
        bool goodfragment=false;
        iFEC++;
	if (debug_s) std::cout <<std::dec<<" end of data fragment record mark "<<std::endl;
        // realign the words skiping event header, the event start 8 words before the first adc mark (0x41505a)
        bool first = true;                                                                                          
	int start_eventRecord=0;
	int end_eventRecord=0;
        for (unsigned int ir=0;ir<nw;ir++){
          if( ((buffer[ir] >> 8) & 0xffffff) == 0x41505a) {
            if (first) {
	      // take out the event header, keep the fec header and the adc payloads
              start_eventRecord = ir - 8;
              if (start_eventRecord < 0) {
		if (debug_s)std::cout <<"  +++ reprocess: Wrong Start of Event"<<std::endl;
		goodfragment = false;
              } else {
		goodfragment = ( ((buffer[ir-6])&0xff) == iFEC);
	      }
              first = false;
            }
          }
	}
        end_eventRecord=nw-start_eventRecord;
        int tFEC=0;
        if (previous.size()>0){
          if (goodfragment ) {
            //  build hits from previous events; 
	    //  decode previous
            //  and copy actual event into previous
            for (unsigned int i=0;i<previous.size();i++){
              event[i]=previous[i];
            }
	    if (debug_s) std::cout <<" +++ reprocess: event # "<<nEvent<<" FEC "<<2-iFEC+1<<std::endl;
            tFEC=2-iFEC+1;
            if (tFEC==1) {
              eventBuilder = new SRSEventBuilder(nEvent,_conf->GetMaxClusterSize(), _conf->GetMinClusterSize(), _conf->GetZeroSupCut(), _conf->GetRunType(), 
						 std::string(_conf->GetClusterPositionCorrectionFlag()).find("applyCorrections")!=std::string(_conf->GetClusterPositionCorrectionFlag()).npos);
            }
            SRSFECDecoder dec(eventBuilder);
            dec.decodeFEC(previous.size(),event);
	    if(debug_s) std::cout << " >>>> Good! Event  "<<std::setw(6)<<prenevt<<" iFEC "<<tFEC<<" # word  "<<previous.size()<<std::endl;
            previous.clear();
            for (unsigned int ir=0; ir < nw-start_eventRecord;ir++){
              previous.push_back(buffer[start_eventRecord+ir]);
            }
            prenevt=nEvent;
          }else{
            //  bad fragment    merge fragment with previous fregment
            iFEC--;
            for (unsigned int ir=0;ir<nw;ir++){
              previous.push_back(buffer[ir]);
            }
            prenevt=nEvent;
            tFEC=iFEC;
            if (tFEC==1) {
              eventBuilder = new SRSEventBuilder(nEvent,_conf->GetMaxClusterSize(), _conf->GetMinClusterSize(), _conf->GetZeroSupCut(), _conf->GetRunType(), 
						 std::string(_conf->GetClusterPositionCorrectionFlag()).find("applyCorrections")!=std::string(_conf->GetClusterPositionCorrectionFlag()).npos);
            }
	    if(debug_s)	std::cout << " >>>> BAD! Merging two fragment in  "<<nEvent<< " tFEC "<<iFEC<<" # word "<<previous.size()<<std::endl;
	    //decode previous                                                                                                                                                                            
            for(unsigned int i=0;i<previous.size();i++){
              event[i]=previous[i];
            }
	    if (debug_s) std::cout <<" +++ reprocess: event # "<<nEvent<<" FEC "<<tFEC<<std::endl;
            SRSFECDecoder dec(eventBuilder);
            dec.decodeFEC(previous.size(),event);
            previous.clear();
          }	  
        }else{
	  if(debug_s) std::cout <<" >>>> SKIP the decoding waiting next fragment"<<std::endl;
          for (unsigned int ir=0; ir < nw-start_eventRecord;ir++){
            previous.push_back(buffer[start_eventRecord+ir]);
          }
          prenevt=nEvent;
        }
        if (tFEC == nFEC){
          // The event is now complete and we can form the cluster
	  if (debug_s) std::cout <<" +++ reprocess: compute cluster"<<std::endl;
          eventBuilder->ComputeClustersInDetectorPlane();
          // Fill the trees
          _root->FillRootFile(eventBuilder);
          eventBuilder=0;
	  
        }
        if (iFEC == 2) {
          newevent = true;
          //  delete eventBuilder;
          iFEC = 0;
        }
        nw = 0;
      } // end of Mark control
    } // good file control...
  }while(!ifs.eof());
  ifs.close();
  if (eventBuilder != 0){
    //    delete eventBuilder;
  }
}



SRSMain::SRSMain(const std::string& rawfile, const std::string& config) : 
	_rawfile(rawfile), 
	_config(config),
	_conf(new SRSConfiguration),
	//_maps(new SRSMapping),
	_maps( SRSMapping::GetInstance() ),
	_root(new SRSOutputROOT),
	_isClosed(false)
  {
  this->Init();
}

void 
SRSMain::Init(){
 //_conf = new SRSConfiguration(_config.c_str());
 _conf->Init(_config.c_str());
 //_maps = SRSMapping::GetInstance();
 _maps->LoadDefaultMapping(_conf->GetMappingFile());
 _maps->LoadAPVtoPadMapping(_conf->GetPadMappingFile());
 _maps->PrintMapping();
 //_root = new SRSOutputROOT(_conf->GetZeroSupCut(), _conf->GetROOTDataType());
 _root->InitRootFile();
 _root->SetZeroSupCut(std::stoi(_conf->GetZeroSupCut()));
 _root->SetROOTDataType(_conf->GetROOTDataType());
}

void SRSMain::Close(){
  std::string strRunName = _rawfile;
  _isClosed = true;
  if( strRunName.find(".raw") != std::string::npos){
	strRunName.erase(strRunName.find(".raw"), strRunName.length() - strRunName.find(".raw") );
  }

  _root->SetRunName(strRunName);
  _root->WriteRootFile();
  //delete _root;
  _root.reset();
  //delete _conf;
  _conf.reset();
  _maps.reset();
}
