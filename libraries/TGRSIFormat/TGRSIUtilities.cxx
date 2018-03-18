#include "TGRSIUtilities.h"

#include <cstdlib>
#include <sys/stat.h>
#include <iostream>
#include <fstream>

#include "TObjArray.h"
#include "TObjString.h"
#include "TPRegexp.h"
#include "TString.h"

bool file_exists(const char* filename)
{
   // std::ifstream(filename);
   struct stat buffer;
   return (stat(filename, &buffer) == 0);
}

bool all_files_exist(const std::vector<std::string>& filenames)
{
   for(auto& filename : filenames) {
      if(!file_exists(filename.c_str())) {
         return false;
      }
   }
   return true;
}

int GetRunNumber(const std::string& fileName)
{
   if(fileName.length() == 0) {
      return 0;
   }
   std::size_t found = fileName.rfind(".root");
   if(found == std::string::npos) {
      return 0;
   }
   std::size_t found2 = fileName.rfind('-');
   // printf("found 2 = %i\n",found2);

   if(found2 == std::string::npos) {
      found2 = fileName.rfind('_');
   }
   std::string temp;
   if(found2 == std::string::npos || fileName.compare(found2 + 4, 5, ".root") != 0) {
      temp = fileName.substr(found - 5, 5);
   } else {
      temp = fileName.substr(found - 9, 5);
   }
   return atoi(temp.c_str());
}

int GetSubRunNumber(const std::string& fileName)
{
   if(fileName.length() == 0) {
      return -1;
   }

   std::size_t found = fileName.rfind('-');
   if(found != std::string::npos) {
      std::string temp = fileName.substr(found + 1, 3);
      return atoi(temp.c_str());
   }
   found = fileName.rfind('_');
   if(found != std::string::npos) {
      std::string temp = fileName.substr(found + 1, 3);
      return atoi(temp.c_str());
   }
   return -1;
}



void trim(std::string& line, const std::string & trimChars) {
   //Removes the the string "trimCars" from the start or end of 'line'
  if(line.length() == 0)
    return;

  std::size_t found = line.find_first_not_of(trimChars);
  if(found != std::string::npos)
    line = line.substr(found, line.length());

  found = line.find_last_not_of(trimChars);
  if(found != std::string::npos)
    line = line.substr(0, found + 1);
}


