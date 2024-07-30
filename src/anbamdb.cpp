#include "anbamdb.hpp"
#include "anbamfilehelper.hpp"
#include "antimestamp.hpp"
#include "sam.h"
#include <sstream>
#include <iostream>

/******************************
 *** SampleIndex defintions ***
 ******************************/
SampleIndex::SampleIndex():offset_l(1), offset_r(0){};

void SampleIndex::_init(const std::string& line)
{
	if(line.substr(0, 2) == "RG") {
		if(line.substr(3,2) == "ID") index2sample.emplace_back(line.substr(6));
		else std::cerr << '(' << antimestamp() << "): [WARNING] unable to parse sample-name from following BAM-header line: " << line << std::endl;
	}
	else if(line.substr(0,2) == "PG"){
		if(line.size() >= 15 && line.substr(0,15) == "PG\tID:otter\tOF:"){
			std::string input = line.substr(15);
	  		std::string value;
	  		std::istringstream stream(input);
	  		std::vector<std::string> columns;
	  		while(std::getline(stream, value, ',')) columns.emplace_back(value);
	  		if(columns.size() == 1) {
	  			offset_l = std::stoi(columns.front());
	  			offset_r = std::stoi(columns.front());
	  		}
	  		else if(columns.size() == 2){
	  			offset_l = std::stoi(columns[0]);
	  			offset_r = std::stoi(columns[1]);
	  		}
	  		else {
	  			std::cerr << '(' << antimestamp() << "): [ERROR] unable to parse offset value from the following BAM-header line: " << line << std::endl;
	  			exit(1);
	  		}
		}
	}
}

void SampleIndex::init(const std::string& bam)
{
	offset_l = 1;
	offset_r = 0;
	BamInstance bam_inst;
	bam_inst.init(bam, true);
	std::string tag;
	for(uint32_t i = 0; i < bam_inst.header->l_text; ++i){
		if(bam_inst.header->text[i] != '@' && bam_inst.header->text[i] != '\n') tag += bam_inst.header->text[i];
		else if(!tag.empty()){
				_init(tag);
				tag.clear();
			}
	}
	if(!tag.empty()) _init(tag);
	if(index2sample.empty()){
		std::cerr << '(' << antimestamp() << "): [ERROR] unable to parse sample-name (read-group) from the following BAM file" << bam << std::endl;
	  	exit(1);
	}
	bam_inst.destroy();
	for(int i = 0; i < (int)index2sample.size(); ++i) sample2index[index2sample[i]] = i;
}
