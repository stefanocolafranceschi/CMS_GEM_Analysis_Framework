//
//  DetectorSuperchamber.h
//  
//
//  Created by Brian L Dorney on 17/11/16.
//
//

#ifndef ____DetectorSuperchamber__
#define ____DetectorSuperchamber__

//C++ Includes
#include <map>
#include <memory>
#include <stdio.h>
#include <string>

//Framework Includes
#include "DetectorTiming.h"
#include "TimingUtilityTypes.h"

//ROOT Includes

namespace QualityControl{
    namespace Timing{
        
        //Declare class
        class DetectorSuperchamber;
        
        //Define detector iterator
        typedef
        std::map<int, std::shared_ptr<Timing::DetectorTiming> >::const_iterator imap_detIter_t;
        
        //Define detector matrix
        //Don't necessarily need to wrap this map in a struct,
        //But in case we need to add other information in the future I have.
        struct DetectorMatrix{
            //Key value is the (i,j) or (i,k) index on the stand
            //The superchamber holds either detectors or pmts
            std::map<std::pair<int,int>, std::shared_ptr<Timing::DetectorSuperchamber> > m_map_detectors;
            
            //Default Destructor
            ~DetectorMatrix(){
                clear();
            }
            
            void clear(){
                m_map_detectors.clear();
            }
            
            std::map<std::pair<int,int>, std::shared_ptr<Timing::DetectorSuperchamber> >::iterator begin(){
                return m_map_detectors.begin();
            }
            
            std::map<std::pair<int,int>, std::shared_ptr<Timing::DetectorSuperchamber> >::iterator end(){
                return m_map_detectors.end();
            }
            
            std::shared_ptr<Timing::DetectorSuperchamber> & operator[](std::pair<int,int> pair_Pos){
                return m_map_detectors[pair_Pos];
            } //End subscript operator
            
            /*DetectorMatrix & operator=(const DetectorMatrix & other){ //assignment operator
                if (this != &other) { //Protects against invalid self-assignment
                    for (auto iterSC = other.m_map_detectors.begin(); iterSC != other.m_map_detectors.end(); ++iterSC) {
                        m_map_detectors[ (*iterSC).first ] = (*iterSC).second;
                    }
                } //Protects against invalid self-assignment
                
                return *this;
            }*/ //End assignment operator
        }; //End DetectorMatrix
        
        //Define class
        class DetectorSuperchamber {
            
        public:
            //Constructors
            //------------------------------------------------------------------------------------------------------------------------------------------
            
            //Destructors
            //------------------------------------------------------------------------------------------------------------------------------------------
            /*~DetectorSuperchamber(){
                reset();
            };*/
            
            //Operators
            //------------------------------------------------------------------------------------------------------------------------------------------
            /*DetectorSuperchamber & operator=(const DetectorSuperchamber & other){
                if (this != &other) { //Protects against invalid self-assignment
                    m_strName = other.m_strName;  //Name of the detector
                    m_map_strBaseAddr2iDetPos = other.m_map_strBaseAddr2iDetPos;
                    detType = other.detType;
                    
                    for (auto iterDet = other.m_map_detectors.begin(); iterDet != other.m_map_detectors.end(); ++iterDet) {
                        m_map_detectors[ (*iterDet).first ] = (*iterDet).second;
                    }
                } //Protects against invalid self-assignment
                
                return *this;
            }*/ //End Assignment operator
            
            //Actions - Methods that Do Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            virtual void reset();
            virtual void resetChannels();
            
            //Getters - Methods that Get (i.e. Return) Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            //Get Detector - by index within superchamber
            virtual std::shared_ptr<Timing::DetectorTiming> getDetector(int iDetPos){
                return m_map_detectors[iDetPos];
            };
            
            //Get Pointer to the beginning of m_map_detectors
            virtual Timing::imap_detIter_t getDetectorPtrBegin(){
                return m_map_detectors.begin();
            };
            
            //Get Pointer to the ending of m_map_detectors
            virtual Timing::imap_detIter_t getDetectorPtrEnd(){
                return m_map_detectors.end();
            };
            
            //Get Detector - by VME base address
            virtual std::shared_ptr<Timing::DetectorTiming> getDetector(std::string & strBaseAddress);
            
            //Get Detector's name
            virtual std::string getName(){
                return m_strName;
            };
            
            //Get Number of detectors
            virtual int getNDetectors(){
                return m_map_detectors.size();
            };
            
            virtual Timing::DetType getType(){
                return detType;
            };
            
            virtual bool hasData();
            
            //Printers - Methods that Print Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            
            //Setters - Methods that Set Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            
            //Set a detector
            virtual void setDetector(int iDetPos, std::shared_ptr<Timing::DetectorTiming> inputDetector);
            
            
            //Set Detector's name
            virtual void setName(std::string strInput){
                m_strName = strInput;
                return;
            };
            
            //Set the detector's type
            virtual void setType(Timing::DetType inputDetType){
                detType = inputDetType;
                return;
            };
            
        protected:
            //Attributes - Methods that Set Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            std::string m_strName;  //Name of the detector
            
            //Key value is detector position within superchamber
            //Mapped value is base address the of the detector
            std::map<std::string, int> m_map_strBaseAddr2iDetPos;
            
            std::map<int, std::shared_ptr<Timing::DetectorTiming> > m_map_detectors;
            
            Timing::DetType detType;
            
        }; //End class DetectorSuperchamber
    } //End namespace Timing
} //End namespace QualityControl


#endif /* defined(____DetectorSuperchamber__) */