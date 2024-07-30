#include <stdint.h>
#include "anbed.hpp"
#include "antimestamp.hpp"
#include <string>
#include <optional>

BED::BED(){};

BED::BED(std::string c, int s, int e): chr(c), start(s), end(e){};

std::string BED::toString() const
{
  return chr + "\t" + std::to_string(start) + "\t" + std::to_string(end);
}

//single-column bed line
std::string BED::toScString() const
{
  return chr + ":" + std::to_string(start) + "-" + std::to_string(end);
}

//standard multi-column bed line
std::optional<BED> parse_bed(const std::string& line){
  std::vector<std::string> columns;
  std::string value;
  std::istringstream stream(line);
  while(std::getline(stream, value, '\t')) columns.emplace_back(value);
  //single-column bed line
  if(columns.size() == 1) return parse_sc_bed(columns.front());
  else if(columns.size() < 3) {
    std::cerr << '(' << antimestamp() << "): Skipping ambiguous BED line: " << line << "\n";
    return std::nullopt;
  }
  else return BED(columns.front(), static_cast<uint32_t>(std::stoul(columns[1])), static_cast<uint32_t>(std::stoul(columns[2])));
}
//single-column bed line
std::optional<BED> parse_sc_bed(const std::string& line)
{
  int index = 0;
  std::string chr;
  int start = -1, end = -1;
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
  if(chr.empty() || start < 0 || end < 0) {
    std::cerr << '(' << antimestamp() << "): Skipping ambiguous multi-BED line: " << line << "\n";
    return std::nullopt;
  }
  else return BED(chr, start, end);
}

void parse_bed_file(const std::string bedfile, std::vector<BED>& bed_vector)
{
  std::ifstream inputbed(bedfile.c_str());
  std::string line;
  // Read the next line from File untill it reaches the end.
  while (std::getline(inputbed, line)) {
    if(line.empty()) std::cerr << '(' << antimestamp() << "): [WARNING] Skipping empty BED line\n";
    else if(line.front() != '#') {
      auto bed = parse_bed(line);
      if(bed) bed_vector.emplace_back(bed.value());
    }
  }
  inputbed.close();

  std::cerr << '(' << antimestamp() << "): Loaded " << bed_vector.size() << " total annotation(s)\n";
}
