#include "bindings/cpp/WFAligner.hpp"
#include "assemble.hpp"
#include "antimestamp.hpp"
#include "anbed.hpp"
#include "anbamfilehelper.hpp"
#include "otter_opts.hpp"
#include "BS_thread_pool.hpp"
#include "anseqs.hpp"
#include "anfahelper.hpp"
#include "analignments.hpp"
#include "andistmat.hpp"
#include "otterclust.hpp"
#include <iostream>
#include <vector>
#include <thread>
#include <map>
#include <mutex>
#include <cmath>

uint32_t count_spanning_reads(const std::vector<ANREAD>& anread_block)
{
	int count = 0;
	for(const auto& read : anread_block) if(read.is_spanning()) ++count;
	return count;
}

void partition_valid_reads(bool& ignore_haps, std::vector<ANREAD>& anread_block, std::vector<int>& valid_indeces, std::vector<int>& invalid_indeces)
{
	for(int i = 0; i < (int)anread_block.size(); ++i) {
		if(!anread_block[i].is_spanning()) invalid_indeces.emplace_back(i);
		else {
			if(ignore_haps) valid_indeces.emplace_back(i);
			else if(anread_block[i].hpt.is_defined()) valid_indeces.emplace_back(i);
			else invalid_indeces.emplace_back(i);
		}
	}
}

void assemble_process(const OtterOpts& params, const std::string& bam, const std::vector<BED>& bed_regions, const std::string& reference, bool& reads_only, BS::thread_pool& pool)
{
	std::cerr << '(' << antimestamp() << "): Processing " << bam << " (" << params.read_group << ')' << std::endl;
	std::mutex std_out_mtx;
	pool.parallelize_loop(0, bed_regions.size(),
		[&pool, &std_out_mtx, &params, &bam, &bed_regions, &reference, &reads_only](const int a, const int b){
			BamInstance bam_inst;
			bam_inst.init(bam, true);
			FaidxInstance faidx_inst;
			if(!reference.empty()) faidx_inst.init(reference);
			wfa::WFAlignerEdit aligner(wfa::WFAligner::Score, wfa::WFAligner::MemoryMed);
			wfa::WFAlignerGapAffine aligner2(4,6,2,wfa::WFAligner::Alignment,wfa::WFAligner::MemoryMed);
			for(int bed_i = a; bed_i < b; ++bed_i) {
				/** set local BED region of interest with offsets **/
				const BED& local_bed = bed_regions[bed_i];
				//std::cout << local_bed.toScString() << std::endl;
				BED mod_bed = local_bed;
				mod_bed.start -= (int)params.offset_l;
				mod_bed.end += (int)params.offset_r;
				/** parse reads, perform local realignments where needed **/
				std::vector<ANREAD> anread_block;
				if(params.is_debug){
					std_out_mtx.lock();
					std::cerr << '(' << antimestamp() << "): [DEBUG] Processing " << local_bed.toScString() << std::endl;
					std_out_mtx.unlock();
				}
				parse_anreads(params, mod_bed, bam_inst, anread_block);
				if(params.is_debug){
					std_out_mtx.lock();
					std::cerr << '(' << antimestamp() << "): [DEBUG] Loaded " << anread_block.size() << " reads" << std::endl;
					std_out_mtx.unlock();
				}
				if((int)anread_block.size() > params.max_cov) std::cerr << '(' << antimestamp() << "): [WARNING] Skipping region with abnormal coverage: " << local_bed.toScString() << " (" << anread_block.size() << ")" << std::endl;
				else{
					if(!reference.empty()) {
						local_realignment(mod_bed.chr, mod_bed.start, mod_bed.end, params.flank, params.min_sim, faidx_inst, anread_block, aligner2);
						if(params.is_debug){
							std_out_mtx.lock();
							std::cerr << '(' << antimestamp() << "): [DEBUG] Locally realigned " << anread_block.size() << " reads" << std::endl;
							std_out_mtx.unlock();
						}
					}
					if(reads_only){
						std_out_mtx.lock();
						for(const auto& read : anread_block) {
							if(params.is_fa) read.stdout_fa(local_bed.toScString());
							else read.stdout_sam(local_bed.chr,local_bed.start, local_bed.end, params.read_group);
						}
						std_out_mtx.unlock();
					}
					else{
						//for(const auto& r : anread_block) std::cout << r.name << '\t' << r.is_spanning_l << '\t' << r.is_spanning_r << '\n' << r.seq << '\n';
						uint32_t spanning_reads = count_spanning_reads(anread_block);
						if(spanning_reads == 0) {
							std_out_mtx.lock();
							std::cerr << '(' << antimestamp() << "): [WARNING] No spanning reads for " << local_bed.toScString() << std::endl;
							std_out_mtx.unlock();
						}
						else{
							/** find valid reads (spanning and haplotagged (if user required)) **/
							//note: mutable to adjust condition
							bool local_ignore_haps = params.ignore_haps;
							std::vector<int> valid_indeces;
							std::vector<int> invalid_indeces;
							partition_valid_reads(local_ignore_haps, anread_block, valid_indeces, invalid_indeces);
							//std::cout << "valid: " << valid_indeces.size() << '\t' << "invalid:" << invalid_indeces.size() << std::endl;
							//not enough haplotagged reads, adjust to 'ignore-haps' mode
							if(valid_indeces.size() < 2){
								local_ignore_haps = true;
								valid_indeces.clear();
								invalid_indeces.clear();
								partition_valid_reads(local_ignore_haps, anread_block, valid_indeces, invalid_indeces);
								if(spanning_reads != valid_indeces.size()){
									std_out_mtx.lock();
									std::cerr << '(' << antimestamp() << "): [ERROR] Unexpected number of valid reads after switching to 'ignore-haps' mode: " << spanning_reads << " vs " << valid_indeces.size() << std::endl;
									std_out_mtx.unlock();
									exit(1);
								}
							}
							if(valid_indeces.empty()) {
								std_out_mtx.lock();
								std::cerr << '(' << antimestamp() << "): [WARNING] No spanning reads for " << local_bed.toScString() << std::endl;
								std_out_mtx.unlock();
							}
							else{
								/** create matrix of pairwise edit-distances **/
								DistMatrix distmatrix(valid_indeces.size());
								if(params.max_alleles != 1) fill_dist_matrix(local_ignore_haps, aligner, anread_block, valid_indeces, distmatrix);
								/** cluster reads and fine allele seqs **/
								ClusteringStatus clustmsg;
								otter_hclust(local_ignore_haps, params.max_alleles, params.bandwidth_short, params.bandwidth_length, params.bandwidth_long, params.max_error, params.min_cov_fraction, params.min_cov_fraction2_l, params.min_cov_fraction2_f, valid_indeces, distmatrix, aligner, anread_block, clustmsg);
								std::vector<int> labels(anread_block.size(), -1);
								for(uint32_t i = 0; i < clustmsg.labels.size(); ++i) {
									labels[valid_indeces[i]] = clustmsg.labels[i];
								}
								//for(uint32_t i = 0; i < labels.size(); ++i) std::cout << i << '\t' << labels[i] << '\t' << anread_block[i].name << std::endl;
								if(!invalid_indeces.empty()){
									//std::cout << "reassigning" << std::endl;
									invalid_reassignment(local_ignore_haps, params.min_sim, params.max_error, clustmsg.fc, anread_block, aligner, labels);
									//for(uint32_t i = 0; i < labels.size(); ++i) std::cout << i << '\t' << labels[i] << '\t' << anread_block[i].name << std::endl;
								}
								std::vector<ANALLELE> alleles(clustmsg.fc);
								rapid_consensus(local_ignore_haps, anread_block, labels, valid_indeces, clustmsg.fc, distmatrix, aligner2, alleles);
								//std::cout << "cosensus" << std::endl;
								std_out_mtx.lock();
								for(int l = 0; l < clustmsg.fc; ++l) {
									alleles[l].ic = clustmsg.ic;
									if(params.is_fa) alleles[l].stdout_fa(params.read_group, local_bed.toScString() + '#' + std::to_string(l));
									else alleles[l].stdout_sam(local_bed.toScString() + "_" + std::to_string(l), local_bed.chr, local_bed.start, local_bed.end, params.read_group);
								}
								std_out_mtx.unlock();
							}
						}
					}
				}
			}
			bam_inst.destroy();
			if(!reference.empty()) faidx_inst.destroy();
	}).wait();
}

void assemble(const std::string& bam, const std::string& bed, const std::string& reference, bool& reads_only, const OtterOpts& params)
{
 	BS::thread_pool pool(params.threads);

 	std::vector<BED> bed_regions;
 	parse_bed_file(bed, bed_regions);

 	if(!params.is_fa){
 		BamInstance bam_inst;
		bam_inst.init(bam, true);
		for(int i = 0; i < bam_inst.header->n_targets; ++i){
			std::cout << "@SQ\tSN:" << bam_inst.header->target_name[i] << "\tLN:" << bam_inst.header->target_len[i] << '\n';
		}
		std::cout << "@RG\tID:" << params.read_group << '\n';
		std::cout << "@PG\tID:otter\tOF:" << params.offset_l << ',' << params.offset_r << '\n';
		//output_preset_offset_tag(params.offset);
		bam_inst.destroy();
 	}
	assemble_process(params, bam, bed_regions, reference, reads_only, pool);
}