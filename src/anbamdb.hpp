#ifndef ANBAMDB_HPP
#define ANBAMDB_HPP

#include "sam.h"
#include <string>
#include <vector>
#include <map>

extern const std::string rg_tag;
extern const std::string pg_tag;
extern const std::string ta_tag;
extern const std::string ac_tag;
extern const std::string tc_tag;

class SampleIndex{
	public:
		int offset_l;
		int offset_r;
		std::vector<std::string> index2sample;
		std::map<std::string, int> sample2index;
		
		SampleIndex();
		void init(const std::string&);
	private:
		void _init(const std::string&);
};

void fetch_preset_offset(
	const std::string&, 
	int&
);

void output_preset_offset_tag(
	const int& offset
);

void parse_bam_allele(
	const std::string&, 
	const int&, 
	const int&, 
	const std::map<std::string,int>&, 
	bam1_t*&, 
	std::vector<std::string>&, 
	std::vector<int>&
);

bool is_multi_sample(
	const std::vector<int>&
);

void set_sample_intervals(
	const std::vector<int>&, 
	std::vector<std::pair<int, std::pair<int, int>>>&
);


#endif