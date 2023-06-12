#include "fa2sam.hpp"
#include "formatter.hpp"
#include "anbed.hpp"
#include "antimestamp.hpp"
#include "kseq.h"
#include <zlib.h>
#include <iostream>

KSEQ_INIT(gzFile, gzread);

void fa2sam(const std::string& ref, const std::string& rg, const std::vector<std::string>& fastas)
{
	int seq_l;
	gzFile fp = gzopen(ref.c_str(), "r");
	kseq_t *seq = kseq_init(fp);
	while ((seq_l = kseq_read(seq)) >= 0) std::cout << "@SQ\tSN:" << seq->name.s << "\tLN:" << seq_l << '\n';
	kseq_destroy(seq);
	gzclose(fp);

	std::cout << "@RG\tID:" << rg << '\n';

	for(const auto& fa : fastas){
		gzFile fp = gzopen(fa.c_str(), "r");
		kseq_t *seq = kseq_init(fp);
		while ((seq_l = kseq_read(seq)) >= 0) {
			const std::string sequence = seq->seq.s;
			const std::string bed = seq->name.s;
			BED local_bed;
			local_bed.parse_multibed(bed);
			int acov = 1;
			int tcov = 1;
			int ic = -1;
			double se = -1.0;
			if(seq->comment.s != nullptr){
				const std::string comment = seq->comment.s;
				std::istringstream stream(comment);
				std::string value;
				int index = 0;
				while(std::getline(stream, value, ' ')){
					if(index == 0) acov = std::stoi(value);
					else if(index == 1) tcov = std::stoi(value);
					else if(index == 2) se = std::stod(value);
					else if(index == 3) se = std::stod(value);
					++index;
				}
			}
			output_fa2sam(bed, local_bed.chr, local_bed.start, local_bed.end, sequence, rg, acov, tcov, -1 , -1, ic,  se);
		}
		kseq_destroy(seq);
		gzclose(fp);
	}

}