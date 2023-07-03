#include "bam_db.hpp"
#include "anbamfilehelper.hpp"
#include "antimestamp.hpp"
#include <htslib/sam.h>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

const std::string rg_tag = "RG";
const std::string ta_tag = "ta";
const std::string ac_tag = "ac";
const std::string tc_tag = "tc";


void index_read_groups(const std::string& bam, std::map<std::string, int>& sample2index, std::vector<std::string>& sample_index)
{
	BamInstance bam_inst;
	bam_inst.init(bam, true);
	const char* result;
	int index = 0;
	result = sam_hdr_line_name(bam_inst.header, rg_tag.c_str(), index++);
	while(result != nullptr){
		std::string local_sample = result;
		sample_index.emplace_back(local_sample);
		result = sam_hdr_line_name(bam_inst.header, rg_tag.c_str(), index++);
	}
	bam_inst.destroy();
	for(int i = 0; i < (int)sample_index.size(); ++i) sample2index[sample_index[i]] = i;
}

void parse_bam_allele(const std::string& target_region, const int& mincov, const std::map<std::string,int>& sample2index, bam1_t*& read, std::vector<std::string>& alleles, std::vector<int>& sample_indeces)
{
	std::string region;
	auto aux_ptr = bam_aux_get(read, ta_tag.c_str());
	if(aux_ptr != NULL) region = (char*)(aux_ptr+1);
	std::cout << target_region << '\t' << region << '\n';
	if(target_region == region){
		aux_ptr = bam_aux_get(read, rg_tag.c_str());
		std::string sample;
		if(aux_ptr != NULL) sample = (char*)(aux_ptr+1);
		std::cout << sample << '\n';
		auto it = sample2index.find(sample);
		if(it == sample2index.end()){
			std::cerr << "(" << antimestamp() << "): ERROR unrecognized read-group tag " << sample << std::endl;
			exit(1);
		}
		aux_ptr = bam_aux_get(read, ac_tag.c_str());
		int ac = 0, tc = 0;
		if(aux_ptr != NULL) ac = bam_aux2i(aux_ptr);
		aux_ptr = bam_aux_get(read, tc_tag.c_str());
		if(aux_ptr != NULL) tc = bam_aux2i(aux_ptr);
		if(ac >= mincov){
		    uint8_t *q = bam_get_seq(read);
		    std::string seq(read->core.l_qseq, 'N');
		    for(int i = 0; i  < read->core.l_qseq; i++) seq[i] = seq_nt16_str[bam_seqi(q, i)];
		    sample_indeces.emplace_back(it->second);
		    alleles.emplace_back(seq);
		}
	}
}

bool is_multi_sample(const std::vector<int>& allele_indeces)
{
	for(int i = 1; i < (int)allele_indeces.size(); ++i) if(allele_indeces[0] != allele_indeces[i]) return true;
	return false;
}

void update_indeces_alleles(const std::vector<int>& init_indeces, const std::vector<int>& sample_indeces, const std::vector<std::string>& alleles, std::vector<int>& _sample_indeces, std::vector<std::string>& _alleles, int& last, int& current)
{
	if(current - last > 2){
		std::cerr << "(" << antimestamp() << "): ERROR more than two alleles found " << std::endl;
		exit(1);
	}
	else {
		if(current - last == 1){
			int total_alleles = 0;
			while(total_alleles++ < 2){
				_sample_indeces.emplace_back(sample_indeces[init_indeces[last]]);
				_alleles.emplace_back(alleles[init_indeces[last]]);

			}
		}
		else{
			for(int i = last; i < current; ++i){
				_sample_indeces.emplace_back(sample_indeces[init_indeces[i]]);
				_alleles.emplace_back(alleles[init_indeces[i]]);
			}
		}
		last = current;
	}
}

void sort_bam_alleles(std::vector<int>& sample_indeces, std::vector<std::string>& alleles)
{

	std::vector<int> init_indeces(sample_indeces.size());
	for(int i = 0; i < (int)sample_indeces.size(); ++i) init_indeces[i] = i;
	sort(init_indeces.begin(), init_indeces.end(), [&sample_indeces](const int& i, const int& j){
		if(sample_indeces[i] == sample_indeces[j]) return i < j;
		else return sample_indeces[i] < sample_indeces[j];
	});
	
	std::vector<int> _sample_indeces;
	std::vector<std::string> _alleles;
	int last = 0;
	int current = 1;
	for(; current < (int)init_indeces.size(); ++current){
		if(sample_indeces[init_indeces[current]] != sample_indeces[init_indeces[last]]) update_indeces_alleles(init_indeces, sample_indeces, alleles, _sample_indeces, _alleles, last, current);
	}

	update_indeces_alleles(init_indeces, sample_indeces, alleles, _sample_indeces, _alleles, last, current);

	/**
	std::cout << _sample_indeces.size() << '\t' << _alleles.size() << '\n';
	for(const auto& i : _sample_indeces) std::cout << i << '\n';
	for(const auto& s :_alleles) std::cout << s << '\n';
	*/


	sample_indeces = _sample_indeces;
	alleles = _alleles;
}

void set_sample_intervals(const std::vector<int>& sample_indeces, std::vector<std::pair<int, std::pair<int, int>>>& sample2interval)
{
	int last = 0;
	int index = last + 1;
	while(index < (int)sample_indeces.size()){
		if(sample_indeces[index] != sample_indeces[last]){
			sample2interval.emplace_back(sample_indeces[last], std::make_pair(last, index - 1));
			last = index;
		}
		++index;
	}
	sample2interval.emplace_back(sample_indeces[last], std::make_pair(last, index - 1));
}


