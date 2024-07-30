#ifndef BED_HPP
#define BED_HPP

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <optional>

class BED {
  public:
    std::string chr;
    uint32_t start;
    uint32_t end;
    BED();
    BED(std::string, int, int);
    std::string toString() const;
    std::string toScString() const;
};

std::optional<BED> parse_bed(const std::string&);
std::optional<BED> parse_sc_bed(const std::string&);

void parse_bed_file(
  const std::string, 
  std::vector<BED>&
);

#endif