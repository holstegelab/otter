#ifndef BAM_DB_HPP
#define BAM_DB_HPP

#include <htslib/sam.h>
#include <string>
#include <vector>
#include <map>

extern const std::string rg_tag;
extern const std::string ta_tag;
extern const std::string ac_tag;
extern const std::string tc_tag;

void index_read_groups(
	const std::string&, 
	std::map<std::string, int>&, 
	std::vector<std::string>&
);

void parse_bam_allele(
	const std::string&, 
	const int&, const std::map<std::string,int>&, 
	bam1_t*&, 
	std::vector<std::string>&, 
	std::vector<int>&
);

bool is_multi_sample(
	const std::vector<int>&
);

void sort_bam_alleles(
	std::vector<int>&, 
	std::vector<std::string>&
);

void set_sample_intervals(
	const std::vector<int>&, 
	std::vector<std::pair<int, std::pair<int, int>>>&
);


#endif