#include "anbed.hpp"
#include "antimestamp.hpp"
#include <string>

BED::BED(){};

BED::BED(std::string c, uint32_t s, uint32_t e): chr(c), start(s), end(e){};

BED::BED(const std::string& line){
  if(line.empty()){
    std::cerr << '(' << antimestamp() << "): Cannot parse empty BED line\n";
    exit(1);
  } else {
      int index = 0;
      std::string value;
      std::istringstream stream(line);
      while(std::getline(stream, value, '\t')){
        switch(index){
          case 0: {
            chr = value;
            break;
          }
          case 1: {
            start = static_cast<uint32_t>(std::stoul(value));
            break;
          }
          case 2: {
            end = static_cast<uint32_t>(std::stoul(value));
            break;
          }
          default: {
            break;
          }
        }
        ++index;
      }       
  }
}

void BED::parse_multibed(const std::string& line)
{
  if(line.empty()){
    fprintf(stderr,"Cannot parse empty BED line\n");
    exit(1);
  } 
  else {
      int index = 0;
      std::string value;
      std::istringstream stream(line);
      while(std::getline(stream, value, ':')) {
        if(index == 0) chr = value;
        else if(index == 1){
          int index2 = 0;
          std::string value2;
          std::istringstream stream2(value);
          while(std::getline(stream2, value2, '-')) {
            if(index2 == 0) start = static_cast<uint32_t>(std::stoul(value2));
            else if(index2 == 1) end = static_cast<uint32_t>(std::stoul(value2));
            ++index2;
          }
        }
        ++index;
      }
  }
}

std::string BED::toString() const
{
  return chr + "\t" + std::to_string(start) + "\t" + std::to_string(end);
}

std::string BED::toBEDstring() const
{
  return chr + ":" + std::to_string(start) + "-" + std::to_string(end);
}


void parse_bed_file(const std::string bedfile, BedMap& bedmap)
{
  std::ifstream inputbed(bedfile.c_str());
  std::string line;
  // Read the next line from File untill it reaches the end.
  while (std::getline(inputbed, line)){
    BED annotation = BED(line);
    //check if chromosome already exists in index
      auto it = bedmap.find(annotation.chr);
      //exists, update key-value 
      if (it != bedmap.end()) {
        it->second.push_back(annotation);
      }
      //does not exist, create key,value entry
      else {
        //set temp vector of coords
        std::vector<BED> localbeds;
        //add current coord to temp vector
        localbeds.push_back(annotation);
        //add vector
        bedmap[annotation.chr] = localbeds;
      }
  }
  inputbed.close();

  int total = 0;
  for(const auto& chr : bedmap) total += (int)chr.second.size();
  std::cerr << '(' << antimestamp() << "): Loaded " << total << " total annotation(s) across " << bedmap.size() << " contig(s)\n";
}