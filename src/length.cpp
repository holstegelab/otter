#include "antimestamp.hpp"
#include "bam_db.hpp"
#include "length.hpp"
#include "anbed.hpp"
#include "anbamfilehelper.hpp"
#include <htslib/sam.h>
#include <string>
#include <vector>
#include <iostream>

void length(const std::string& bam, const std::string& bed, const int& ac_mincov, const int& tc_mincov)
{
	std::vector<BED> regions;
	parse_bed_file(bed, regions);

	std::cerr << '(' << antimestamp() << "): Loaded " << regions.size() << " total regions\n";

	BamInstance bam_inst;
	bam_inst.init(bam, true);

	std::vector<std::string> sample_index;
	std::map<std::string, int> sample2index;
	index_read_groups(bam, sample2index, sample_index);

	std::cerr << '(' << antimestamp() << "): Found " << sample_index.size() << " samples (read-group tags)\n";

	for(const auto& region : regions){
		const std::string region_str = region.chr + ':' + std::to_string(region.start) + '-' + std::to_string(region.end);
		hts_itr_t *iter = sam_itr_querys(bam_inst.idx, bam_inst.header, region_str.c_str());
		if(iter == nullptr) std::cerr << "(" << antimestamp() << "): WARNING: query failed at region " << region_str << std::endl;
		else{
			std::vector<int> allele_sample_indeces;
			std::vector<std::string> alleles;
			while(sam_itr_next(bam_inst.fp, iter, bam_inst.read) > 0) parse_bam_allele(region_str, ac_mincov, tc_mincov, sample2index, bam_inst.read, alleles, allele_sample_indeces);
			for(int i = 0; i < (int)allele_sample_indeces.size(); ++i) std::cout << region_str << '\t' << sample_index[allele_sample_indeces[i]] << '\t' << alleles[i].size() << '\n';
		}
	}

	bam_inst.destroy();

}