#include "anbamfilehelper.hpp"
#include "antimestamp.hpp"
#include "sam.h"
#include <string>
#include <map>
#include <vector>
#include <thread>
#include <iostream>
#include <mutex>

#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2))

void BamInstance::init(const std::string& file, bool _index)
{
	index = _index;
	fp = bgzf_open(file.c_str(), "r");
	header = bam_hdr_read(fp);
	if(header == nullptr) std::cout << "error\n";
	//load index
	if(index) idx = bam_index_load((file + ".bai").c_str());
	//init read
	read = bam_init1();
};

void BamInstance::destroy()
{
	bgzf_close(fp);
	fp = nullptr;
	if(index) hts_idx_destroy(idx);
	idx = nullptr;
	bam_hdr_destroy(header);
	header = nullptr;
	bam_destroy1(read);
	read = nullptr;
};