//
//  TimingUtilityFunctions.cpp
//  
//
//  Created by Brian L Dorney on 04/11/15.
//
//

#include "TimingUtilityFunctions.h"

using std::cout;
using std::endl;

using namespace Timing;

//Function Code Begins
//========================================================================================
//boolean
//----------------------------------------------------------------------------------------
bool QualityControl::Timing::convert2bool(std::string str, bool &bExitSuccess){
    //convert to upper case
    std::transform(str.begin(), str.end(), str.begin(), toupper);
    
    //Empty String?
    if ( str.empty() ) {
        bExitSuccess = false;
        return false;
    }
    
    //Input recognized?
    if ( !(str.compare("T")==0 || str.compare("TRUE")==0 || str.compare("1")==0
           || str.compare("F")==0 || str.compare("FALSE")==0 || str.compare("0")==0) ) {
        bExitSuccess = false;
        return false;
    }
    
    bExitSuccess = true;
    
    return (str.compare("T")==0 || str.compare("TRUE")==0 || str.compare("1")==0);
} //End convert2bool()

//file I/O
//----------------------------------------------------------------------------------------
std::istream & QualityControl::Timing::getlineNoSpaces(std::istream & stream, std::string & str){
    //get the line
    std::getline(stream, str);
    
    if (str.find("'",0) == std::string::npos ) { //Case: Header
        str.erase(remove(str.begin(),str.end(), ' ' ), str.end() ); //spaces
        str.erase(remove(str.begin(),str.end(), '\t' ), str.end() );//tabs
    } //End Case: Header
    else{ //Case: line input
        int iFirstQuote = str.find("'",0);
        int iLastQuote = str.rfind("'");
        
        //Spaces
        if (str.find(' ') != std::string::npos ) { //Case: Spaces Found
            str.erase(remove(str.begin(),str.begin()+iFirstQuote, ' ' ), str.begin()+iFirstQuote ); //Until first single-quote
            str.erase(remove(str.begin()+iLastQuote,str.end(), ' ' ), str.end() ); //Until last single-quote
        } //End Case: Spaces Found
        
        //tabs
        if(str.find('\t') != std::string::npos){ //Case: Tabs Found
            str.erase(remove(str.begin(),str.begin()+iFirstQuote, '\t' ), str.begin()+iFirstQuote ); //Until first single-quote
            str.erase(remove(str.begin()+iLastQuote,str.end(), '\t' ), str.end() ); //Until last single-quote
        } //End Case: Tabs Found
    } //End Case: line input
    
    return stream;
} //End getlineNoSpaces()

//Math
//----------------------------------------------------------------------------------------
//T -> Type; A -> Allocator
//Determines the difference in mean between two datasets
//template<typename T, typename A>
//float QualityControl::Timing::deltaMean( std::vector<T,A> const &vec1, std::vector<T,A> const &vec2){
//    float fMean1 = (std::accumulate(vec1.begin(), vec1.end(), 0.0) / (float) vec1.size() );
//    float fMean2 = (std::accumulate(vec2.begin(), vec2.end(), 0.0) / (float) vec2.size() );
//
//    return fMean1 - fMean2;
//} //End deltaMean

int QualityControl::Timing::getDeltaTForChannel(std::map<std::string, int> inputMap){
    int max = getMaxForChannelAND(inputMap);
    int min = getMinForChannelOR(inputMap);
    
    //cout<<"QualityControl::Timing::getDeltaTForChannel() - max = " << max << "; min = " << min << endl;
    
    if (max > 0 && min > 0 ) {
        return abs(max - min);
    }
    else{
        return -1;
    }
} //End getDeltaTForChannel()

//Gets the maximum value for two channels (both channels required to be nonzero)
int QualityControl::Timing::getMaxForChannelAND(std::map<std::string, int> inputMap){
    //Variable Declaration
    //int iRetVal;
    
    std::pair<std::string, int> min = *min_element(inputMap.begin(), inputMap.end(), CompareSecond_Min() );
    //std::pair<std::string, int> min = *min_element(inputMap.begin(), inputMap.end() );
    
    //Require All Elements to be nonzero (i.e. have a signal)
    if ( min.second > 0 ) {
        //std::pair<std::string, int> max = *max_element(inputMap.begin(), inputMap.end(), CompareSecond_Max() );
        std::pair<std::string, int> max = *max_element(inputMap.begin(), inputMap.end(), CompareSecond_Min() );
        //std::pair<std::string, int> max = *max_element(inputMap.begin(), inputMap.end() );
        
        //iRetVal = getMaxForChannelAND(inputMap);
        //iRetVal = max.second;
        return max.second;
    }
    else{
        //One or more channels off
        //iRetVal = -1;
        return -1;
    }
    //return iRetVal;
} //End getMaxForChannel

//Gets the minimum value for two channels
int QualityControl::Timing::getMinForChannelOR(std::map<std::string, int> inputMap){
    //Variable Declaration
    std::map<std::string, int>::iterator iterMap = inputMap.begin();
    std::map<std::string, int>::iterator iterMapEnd = inputMap.end();
    
    while( iterMap != inputMap.end() ){
        if( 0 == (*iterMap).second){
            //C++11 is magic
            iterMap = inputMap.erase(iterMap);
        }
        else{
            ++iterMap;
        }
    }
    
    if( 0 == inputMap.size() ){
        return -1;
    }
    else{
        std::pair<std::string, int> min = *std::min_element(inputMap.begin(), inputMap.end(), CompareSecond_Min());
        //std::pair<std::string, int> min = *std::min_element(inputMap.begin(), inputMap.end() );
        
        return min.second;
    }
}

//printers
//----------------------------------------------------------------------------------------
//Prints the status bits of an input TFile
void QualityControl::Timing::printROOTFileStatus(TFile * file_ROOT){
    std::cout << "Input ROOT File Status:\n";
    std::cout << "\tIsOpen() = " << file_ROOT->IsOpen() << std::endl;
    std::cout << "\tIsZombie() = " << file_ROOT->IsZombie() << std::endl;
    
    return;
} //End printROOTFileStatus

//Prints an error message to the user if a specific pattern (strPatternNotFound) is not found during file I/O
//Also informs the user which calling method and file the pattern is not found under
void QualityControl::Timing::printStringNotFoundMsg(std::string strCallingMethod, std::string strPatternNotFound, std::string strLine, std::string strFilename){
    std::cout<<"String '" << strPatternNotFound << "' Not Found in line:\n";
    std::cout<< strLine <<std::endl;
    std::cout<<"Exiting " << strCallingMethod << ", Cross-Check input File:\n";
    std::cout<<strFilename<<std::endl;
    
    return;
} //End printStringNotFoundMsg()

//Prints All Bit Flags for an input ifstream
void QualityControl::Timing::printStreamStatus(std::ifstream &inputStream){
    std::cout << "Input File Stream Bit Status:\n";
    std::cout << " good()=" << inputStream.good() << std::endl;
    std::cout << " eof()=" << inputStream.eof() << std::endl;
    std::cout << " fail()=" << inputStream.fail() << std::endl;
    std::cout << " bad()=" << inputStream.bad() << std::endl;
    
    return;
} //End treeProducerTDC::printStreamStatus()

//Prints All Bit Flags for an input fstream
void QualityControl::Timing::printStreamStatus(std::fstream &inputStream){
    std::cout << "Input File Stream Bit Status:\n";
    std::cout << " good()=" << inputStream.good() << std::endl;
    std::cout << " eof()=" << inputStream.eof() << std::endl;
    std::cout << " fail()=" << inputStream.fail() << std::endl;
    std::cout << " bad()=" << inputStream.bad() << std::endl;
    
    return;
} //End treeProducerTDC::printStreamStatus()

//string manipulation
//----------------------------------------------------------------------------------------
//Safe conversion to float
//Scientific notation not supported
float QualityControl::Timing::stofSafe(std::string strInputValue){
    bool bInputUnderstood = false;
    bool bManualEntry = false;
    
    float ret_float = -1;
    
    std::string strUserInput;
    
    if (strInputValue.find_first_not_of( "-0123456789." ) == std::string::npos) { //Case: only numeric data
        ret_float =  std::stof(strInputValue);
    } //End Case: only numeric data
    else{ //Case: non numeric data entered
        std::cout<<"QualityControl::Timing::stofSafe() - Sorry numeric conversion failed\n";
        std::cout<<"\tValue = " << strInputValue << std::endl;
        std::cout<<"\tNonumeric characters present (sorry scientific notation not supported)\n";
        std::cout<<"\tWould you like to enter a value manually [y/N]?\n";
        
        while (!bInputUnderstood) { //Loop requesting manual entry
            std::cin>>strUserInput;
            
            std::transform(strUserInput.begin(), strUserInput.end(), strUserInput.begin(), toupper);
            
            if (0 == strUserInput.compare("Y") || 0 == strUserInput.compare("YES") || 0 == strUserInput.compare("1") ) { //Case: User wants to correct via manual entry
                bInputUnderstood = true;
                bManualEntry = true;
            } //End Case: User wants to correct via manual entry
            else if ( 0 == strUserInput.compare("N") || 0 == strUserInput.compare("NO") || 0 == strUserInput.compare("0") ){ //Case: User does not wish for manual entry
                bInputUnderstood = true;
                bManualEntry = false;
            } //End Case: User does not wish for manual entry
            else{ //Case: User Input Not Understood
                std::cout<<"\tSorry your input was not understood\n";
                std::cout<<"\tPlease answer from the set {y,yes,1,N,NO,0}\n";
                std::cout<<"\tWould you like to enter a value manually {y,yes,1,N,NO,0}?\n";
            } //End Case: User Input Not Understood
        } //End Loop requesting manual entry
        
        if (bManualEntry) { //Case: Manual Entry
            //std::cout<<"\tOkay, for Field = " << strInputField << " Enter a Numeric Value:\n";
            std::cout<<"\tOkay, Enter a Numeric Value:\n";
            std::cin>>strUserInput;
            
            ret_float = stofSafe(strUserInput); //Recursion
        } //End Case: Manual Entry
    } //End Case: non numeric data entered
    
    return ret_float;
} //End QualityControl::Timing::stofSafe()

float QualityControl::Timing::stofSafe(std::string strInputField, std::string strInputValue){
    bool bInputUnderstood = false;
    bool bManualEntry = false;
    
    float ret_float = -1;
    
    std::string strUserInput;
    
    if (strInputValue.find_first_not_of( "-0123456789." ) == std::string::npos) { //Case: only numeric data
        ret_float =  std::stof(strInputValue);
    } //End Case: only numeric data
    else{ //Case: non numeric data entered
        std::cout<<"QualityControl::Timing::stofSafe() - Sorry numeric conversion failed\n";
        std::cout<<"\tField = " << strInputField << std::endl;
        std::cout<<"\tValue = " << strInputValue << std::endl;
        std::cout<<"\tNonumeric characters present (sorry scientific notation not supported)\n";
        std::cout<<"\tWould you like to enter a value manually [y/N]?\n";
        
        while (!bInputUnderstood) { //Loop requesting manual entry
            std::cin>>strUserInput;
            
            std::transform(strUserInput.begin(), strUserInput.end(), strUserInput.begin(), toupper);
            
            if (0 == strUserInput.compare("Y") || 0 == strUserInput.compare("YES") || 0 == strUserInput.compare("1") ) { //Case: User wants to correct via manual entry
                bInputUnderstood = true;
                bManualEntry = true;
            } //End Case: User wants to correct via manual entry
            else if ( 0 == strUserInput.compare("N") || 0 == strUserInput.compare("NO") || 0 == strUserInput.compare("0") ){ //Case: User does not wish for manual entry
                bInputUnderstood = true;
                bManualEntry = false;
            } //End Case: User does not wish for manual entry
            else{ //Case: User Input Not Understood
                std::cout<<"\tSorry your input was not understood\n";
                std::cout<<"\tPlease answer from the set {y,yes,1,N,NO,0}\n";
                std::cout<<"\tWould you like to enter a value manually {y,yes,1,N,NO,0}?\n";
            } //End Case: User Input Not Understood
        } //End Loop requesting manual entry
        
        if (bManualEntry) { //Case: Manual Entry
            std::cout<<"\tOkay, for Field = " << strInputField << " Enter a Numeric Value:\n";
            std::cin>>strUserInput;
            
            ret_float = stofSafe(strInputField, strUserInput); //Recursion
        } //End Case: Manual Entry
    } //End Case: non numeric data entered
    
    return ret_float;
} //End QualityControl::Timing::stofSafe()

int QualityControl::Timing::stoiSafe(std::string strInputValue){
    bool bInputUnderstood = false;
    bool bManualEntry = false;
    
    int ret_int = -1;
    
    std::string strUserInput;
    
    if (strInputValue.find_first_not_of( "-0123456789" ) == std::string::npos) { //Case: only numeric data
        ret_int =  std::stoi(strInputValue);
    } //End Case: only numeric data
    else{ //Case: non numeric data entered
        std::cout<<"QualityControl::Timing::stoiSafe() - Sorry numeric conversion failed\n";
        std::cout<<"\tValue = " << strInputValue << std::endl;
        std::cout<<"\tNonumeric characters present (or maybe you gave a floating point number instead?)\n";
        std::cout<<"\tWould you like to enter a value manually [y/N]?\n";
        
        while (!bInputUnderstood) { //Loop requesting manual entry
            std::cin>>strUserInput;
            
            std::transform(strUserInput.begin(), strUserInput.end(), strUserInput.begin(), toupper);
            
            if (0 == strUserInput.compare("Y") || 0 == strUserInput.compare("YES") || 0 == strUserInput.compare("1") ) { //Case: User wants to correct via manual entry
                bInputUnderstood = true;
                bManualEntry = true;
            } //End Case: User wants to correct via manual entry
            else if ( 0 == strUserInput.compare("N") || 0 == strUserInput.compare("NO") || 0 == strUserInput.compare("0") ){ //Case: User does not wish for manual entry
                bInputUnderstood = true;
                bManualEntry = false;
            } //End Case: User does not wish for manual entry
            else{ //Case: User Input Not Understood
                std::cout<<"\tSorry your input was not understood\n";
                std::cout<<"\tPlease answer from the set {y,yes,1,N,NO,0}\n";
                std::cout<<"\tWould you like to enter a value manually {y,yes,1,N,NO,0}?\n";
            } //End Case: User Input Not Understood
        } //End Loop requesting manual entry
        
        if (bManualEntry) { //Case: Manual Entry
            std::cout<<"\tOkay, Enter a Numeric Value:\n";
            std::cin>>strUserInput;
            
            ret_int = stoiSafe(strUserInput); //Recursion
        } //End Case: Manual Entry
    } //End Case: non numeric data entered
    
    return ret_int;
} //End QualityControl::Timing::stoiSafe()

int QualityControl::Timing::stoiSafe(std::string strInputField, std::string strInputValue){
    bool bInputUnderstood = false;
    bool bManualEntry = false;
    
    int ret_int = -1;
    
    std::string strUserInput;
    
    if (strInputValue.find_first_not_of( "-0123456789" ) == std::string::npos) { //Case: only numeric data
        ret_int =  std::stoi(strInputValue);
    } //End Case: only numeric data
    else{ //Case: non numeric data entered
        std::cout<<"QualityControl::Timing::stoiSafe() - Sorry numeric conversion failed\n";
        std::cout<<"\tField = " << strInputField << std::endl;
        std::cout<<"\tValue = " << strInputValue << std::endl;
        std::cout<<"\tNonumeric characters present (or maybe you gave a floating point number instead?)\n";
        std::cout<<"\tWould you like to enter a value manually [y/N]?\n";
        
        while (!bInputUnderstood) { //Loop requesting manual entry
            std::cin>>strUserInput;
            
            std::transform(strUserInput.begin(), strUserInput.end(), strUserInput.begin(), toupper);
            
            if (0 == strUserInput.compare("Y") || 0 == strUserInput.compare("YES") || 0 == strUserInput.compare("1") ) { //Case: User wants to correct via manual entry
                bInputUnderstood = true;
                bManualEntry = true;
            } //End Case: User wants to correct via manual entry
            else if ( 0 == strUserInput.compare("N") || 0 == strUserInput.compare("NO") || 0 == strUserInput.compare("0") ){ //Case: User does not wish for manual entry
                bInputUnderstood = true;
                bManualEntry = false;
            } //End Case: User does not wish for manual entry
            else{ //Case: User Input Not Understood
                std::cout<<"\tSorry your input was not understood\n";
                std::cout<<"\tPlease answer from the set {y,yes,1,N,NO,0}\n";
                std::cout<<"\tWould you like to enter a value manually {y,yes,1,N,NO,0}?\n";
            } //End Case: User Input Not Understood
        } //End Loop requesting manual entry
        
        if (bManualEntry) { //Case: Manual Entry
            std::cout<<"\tOkay, for Field = " << strInputField << " Enter a Numeric Value:\n";
            std::cin>>strUserInput;
            
            ret_int = stoiSafe(strInputField, strUserInput); //Recursion
        } //End Case: Manual Entry
    } //End Case: non numeric data entered
    
    return ret_int;
} //End QualityControl::Timing::stoiSafe()

std::pair<std::string,std::string> QualityControl::Timing::getParsedLine(std::string strInputLine, bool &bExitSuccess){
    //Variable Declaration
    int iPos_Equals = strInputLine.find("=",0);
    int iPos_End    = strInputLine.find(";",0);
    int iPos_Quote1 = strInputLine.find("'",0); //Position of first single-quote
    //int iPos_Quote2 = strInputLine.rfind("'");  //Position of last single-quote
    int iPos_Quote2 = strInputLine.find("'",iPos_Quote1+1); //Position of last single-quote
    
    std::string strFieldName;
    std::string strFieldValue;
    
    //Check to make sure iPos_Equals found in input string
    if (iPos_Equals == std::string::npos && strInputLine.find("[",0) == std::string::npos ) {
        printStringNotFoundMsg("TimingUtilityFunctions::getParsedLine","=",strInputLine,"See Next Error Msg for Filename" );
        
        bExitSuccess = false;
        
        return std::make_pair("","");
    }
    
    //Check to make sure iPos_End found in input string
    if (iPos_End == std::string::npos && strInputLine.find("[",0) == std::string::npos) {
        printStringNotFoundMsg("TimingUtilityFunctions::getParsedLine",";",strInputLine,"See Next Error Msg for Filename" );
        
        bExitSuccess = false;
        
        return std::make_pair("","");
    }
    
    //Check to make sure iPos_Quote1 found in input string
    if (iPos_Quote1 == std::string::npos && strInputLine.find("[",0) == std::string::npos) {
        printStringNotFoundMsg("TimingUtilityFunctions::getParsedLine","'",strInputLine,"See Next Error Msg for Filename" );
        
        bExitSuccess = false;
        
        return std::make_pair("","");
    }
    
    //Check to make sure iPos_Quote2 found in input string
    if (iPos_Quote2 == std::string::npos && strInputLine.find("[",0) == std::string::npos) {
        printStringNotFoundMsg("TimingUtilityFunctions::getParsedLine","'",strInputLine,"See Next Error Msg for Filename" );
        
        bExitSuccess = false;
        
        return std::make_pair("","");
    }
    
    //Get the Strings
    strFieldName    = strInputLine.substr(0,iPos_Equals);
    strFieldValue   = strInputLine.substr(iPos_Quote1+1, iPos_Quote2 - iPos_Quote1 - 1);
    
    //Set the Exit Flag
    bExitSuccess = true;
    
    //Return the Pair
    return std::make_pair(strFieldName,strFieldValue);
} //End getParsedLine()

//template method must be defined in the header file
//template<class TConversion>
//std::string QualityControl::Timing::getString(TConversion input){
//    std::stringstream sstream;
//    sstream<<input;
//    return sstream.str();
//} //End getString()

//Returns a vector of strings taken from a comma separated list; ignores whitespaces in between elements
std::vector<std::string> QualityControl::Timing::getCharSeparatedList(std::string strInput, char cDelimiter){
    //Variable Declaration
    std::vector<std::string> vec_retStrings;
    
    //Remove all white spaces and tabs
    if( !(strInput.find(' ') == std::string::npos) ) strInput.erase(remove(strInput.begin(), strInput.end(), ' '), strInput.end() );
    if( !(strInput.find('\t')== std::string::npos) ) strInput.erase(remove(strInput.begin(), strInput.end(), '\t'), strInput.end() );
    
    if (strInput.find(cDelimiter) == std::string::npos ) { //Case: Comma's not found
        vec_retStrings.push_back(strInput);
    } //End Case: Comma's not found
    else{ //Case: Comma separated list given!
        //Declare Variables
        std::stringstream sstream(strInput);
        
        std::string strEntry;
        
        //Give the input string to the string stream
        //sstream<<strInput;
        
        //Get Entries
        while ( std::getline(sstream, strEntry, cDelimiter) ){ //Loop through the string
            vec_retStrings.push_back(strEntry);
        } //End loop through the string
    } //End Case: Comma separated list given!
    
    return vec_retStrings;
} //End getCommaSeparatedList()

//Returns a vector of strings taken from a comma separated list; ignores whitespaces in between elements
/*std::vector<std::string> QualityControl::Timing::getCommaSeparatedList(std::string strInput){
 //Variable Declaration
 std::vector<std::string> vec_retStrings;
 
 //Remove all white spaces and tabs
 if( !(strInput.find(' ') == std::string::npos) ) strInput.erase(remove(strInput.begin(), strInput.end(), ' '), strInput.end() );
 if( !(strInput.find('\t')== std::string::npos) ) strInput.erase(remove(strInput.begin(), strInput.end(), '\t'), strInput.end() );
 
 if (strInput.find(",") == std::string::npos ) { //Case: Comma's not found
 vec_retStrings.push_back(strInput);
 } //End Case: Comma's not found
 else{ //Case: Comma separated list given!
 //Declare Variables
 std::stringstream sstream(strInput);
 
 std::string strEntry;
 
 //Give the input string to the string stream
 //sstream<<strInput;
 
 //Get Entries
 while ( std::getline(sstream, strEntry, ',') ){ //Loop through the string
 vec_retStrings.push_back(strEntry);
 } //End loop through the string
 } //End Case: Comma separated list given!
 
 return vec_retStrings;
 }*/ //End getCommaSeparatedList()

//Unsorted
//----------------------------------------------------------------------------------------

int QualityControl::Timing::getCyclicColor(int iInput){
    int ret_int;
    
    switch (iInput % 8) {
        case 0:     ret_int = kBlack;     break;
        case 1:     ret_int = kRed+1;     break;
        case 2:     ret_int = kGreen+1;   break;
        case 3:     ret_int = kBlue+1;    break;
        case 4:     ret_int = kCyan; break;
        case 5:     ret_int = kYellow+2;  break;
        case 6:     ret_int = kOrange;  break;
        case 7:     ret_int = kPink;  break;
        case 8:     ret_int = kMagenta;  break;
        default:
            break;
    }
    
    return ret_int;
} //End getCyclicColor()

