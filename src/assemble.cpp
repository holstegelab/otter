#include "assemble.hpp"
#include "antimestamp.hpp"
#include "anbamfilehelper.hpp"
#include "parse_bam_alignments.hpp"
#include "anbed.hpp"
#include "otter_opts.hpp"
#include "pairwise_alignment.hpp"
#include "fasta_helper.hpp"
#include "BS_thread_pool.hpp"
#include <htslib/faidx.h>
#include "bindings/cpp/WFAligner.hpp"
#include <iostream>
#include <vector>
#include <thread>
#include <map>
#include <mutex>


void general_process(const OtterOpts& params, const std::string& bam, const std::vector<BED>& bed_regions, const std::string& reference, BS::thread_pool& pool)
{
	
	std::cerr << '(' << antimestamp() << "): Processing " << bam << '\n';
	std::mutex std_out_mtx;
	pool.parallelize_loop(0, bed_regions.size(),
		[&pool, &std_out_mtx, &params, &bam, &bed_regions, &reference](const int a, const int b){
			BamInstance bam_inst;
			bam_inst.init(bam, true);
			wfa::WFAlignerEdit aligner(wfa::WFAligner::Score, wfa::WFAligner::MemoryMed);
			wfa::WFAlignerGapAffine aligner2(4,6,2,wfa::WFAligner::Alignment,wfa::WFAligner::MemoryMed);
			FaidxInstance faidx_inst;
			if(!reference.empty()) faidx_inst.init(reference);

			for(int i = a; i < b; ++i) {
				AlignmentBlock alignment_block;
				const BED& local_bed = bed_regions[i];
				BED mod_bed = local_bed;
				mod_bed.start -= (uint32_t)params.offset;
				mod_bed.end += (uint32_t)params.offset;
				parse_alignments(params, mod_bed, bam_inst, alignment_block);
				if(alignment_block.names.size() == 0) std::cerr << '(' << antimestamp() << "): WARNING: no reads found at region " << local_bed.toBEDstring() << '\n';
				else if ((int)alignment_block.names.size() > params.max_cov) std::cerr << '(' << antimestamp() << "): WARNING: abnormal coverage for region " << local_bed.toBEDstring() << " (" << alignment_block.names.size() << "x)\n";
				else{
					if(!reference.empty()) otter_realignment(local_bed.chr, (int)mod_bed.start, (int)mod_bed.end, params.flank, params.min_sim, faidx_inst, alignment_block, aligner2);
					std::vector<int> spannable_indeces;
					for(int j = 0; j < (int)alignment_block.statuses.size(); ++j) if(alignment_block.statuses[j].is_spanning()) spannable_indeces.emplace_back(j);
					if(!spannable_indeces.empty()) {
						//for(int j = 0; j < (int)alignment_block.names.size(); ++j) std::cout << '>' << alignment_block.names[j] << ' ' << alignment_block.statuses[j].is_spanning() << '\n' << alignment_block.seqs[j] << '\n';
						DistMatrix distmatrix(spannable_indeces.size());
						std::vector<int> labels;
						otter_hclust(params.max_alleles, params.bandwidth, params.max_error, spannable_indeces, distmatrix, aligner, alignment_block, labels);
						//for(int j = 0; j < (int)alignment_block.names.size(); ++j) std::cout << j << '\t' << labels[j] << '\n';
						std::vector<std::string> consensus_seqs;
						otter_rapid_consensus(spannable_indeces, labels, distmatrix, aligner2, alignment_block, consensus_seqs);
						otter_nonspanning_assigment(params.min_sim, params.max_error, alignment_block, aligner, labels);
						std_out_mtx.lock();
						//for(int j = 0; j < (int)alignment_block.names.size(); ++j) std::cout << j << '\t' << labels[j] << '\t' << alignment_block.statuses[j].spanning_l << '\t' << alignment_block.statuses[j].spanning_r << '\t' << alignment_block.statuses[j].alignment_coords.first << '\t' << alignment_block.statuses[j].alignment_coords.second << '\n';
						for(int j = 0; j < (int)consensus_seqs.size(); ++j){
							int cov = 0;
							for(const auto l : labels) if(l == j) ++cov;
							std::cout << '>' << bed_regions[i].toBEDstring() << ' ' << cov << ' ' << alignment_block.names.size() <<  '\n';
							std::cout << consensus_seqs[j] << '\n';
						}
						std_out_mtx.unlock();
					}
				}
			}
			bam_inst.destroy();
			if(!reference.empty()) faidx_inst.destroy();
	}).wait();
}


void assemble(const std::vector<std::string>& bams, const std::string& bed, const std::string& reference, const OtterOpts& params)
{
 	BS::thread_pool pool(params.threads);

 	//load bed file
 	BedMap map_beds;
 	parse_bed_file(bed, map_beds);

 	std::vector<BED> bed_regions;
 	for(const auto& chr : map_beds) for(const auto& bed : chr.second) bed_regions.emplace_back(bed);
 	for(const auto& bam : bams) general_process(params, bam, bed_regions, reference, pool);
}