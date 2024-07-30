#include "anfahelper.hpp"
#include "faidx.h"
#include <string>
#include <iostream>

void FaidxInstance::init(const std::string& reference){faidx_ptr = fai_load(reference.c_str());}

int FaidxInstance::fetch(const std::string& chr, int start, int end, std::string& seq) const
{
	int ref_l;
	char* ref_seq = faidx_fetch_seq(faidx_ptr, chr.c_str(), start, end, &ref_l);
	if(ref_l > 0){
		seq.resize(ref_l, 'N');
		for(int i = 0; i < ref_l; ++i) seq[i] = std::toupper(*(ref_seq +i));
	}
	free(ref_seq);
	return ref_l;
}

void FaidxInstance::destroy(){fai_destroy(faidx_ptr); faidx_ptr = nullptr;}