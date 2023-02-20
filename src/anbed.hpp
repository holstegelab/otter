#ifndef BED_HPP
#define BED_HPP

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

class BED {
  public:
    std::string chr;
    uint32_t start;
    uint32_t end;
    BED(std::string, uint32_t, uint32_t);
    BED(const std::string&);
    std::string toString() const;
    std::string toBEDstring() const;
};

//type alias for map -> (chrm index, vector of Region objects)
typedef std::map<std::string, std::vector<BED>> BedMap;

void parse_bed_file(
  const std::string, 
  BedMap&
);

#endif