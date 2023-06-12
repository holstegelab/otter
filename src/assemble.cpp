#include "assemble.hpp"
#include "antimestamp.hpp"
#include "anbamfilehelper.hpp"
#include "parse_bam_alignments.hpp"
#include "anbed.hpp"
#include "otter_opts.hpp"
#include "pairwise_alignment.hpp"
#include "fasta_helper.hpp"
#include "formatter.hpp"
#include "BS_thread_pool.hpp"
#include <htslib/faidx.h>
#include "bindings/cpp/WFAligner.hpp"
#include <iostream>
#include <vector>
#include <thread>
#include <map>
#include <mutex>
#include <cmath>

double compute_se(const std::vector<double> values)
{
	if(values.empty()) return -1.0;
	else{
		double u = 0.0, n = 0.0;
		for(const auto& v : values) u += v;
		u /= values.size();
		for(const auto& v : values) n += std::pow(v - u, 2.0);
		return std::sqrt(n/(values.size() - 1));
	}
	
}

void general_process(const OtterOpts& params, const std::string& bam, const std::vector<BED>& bed_regions, const std::string& reference, const bool& reads_only, BS::thread_pool& pool)
{
	
	std::cerr << '(' << antimestamp() << "): Processing " << bam << '\n';
	std::mutex std_out_mtx;
	pool.parallelize_loop(0, bed_regions.size(),
		[&pool, &std_out_mtx, &params, &bam, &bed_regions, &reference, &reads_only](const int a, const int b){
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
					std::string region_str = local_bed.toBEDstring();
					if(!reference.empty()) otter_realignment(local_bed.chr, (int)mod_bed.start, (int)mod_bed.end, params.flank, params.min_sim, faidx_inst, alignment_block, aligner2);
					if(reads_only){
						std_out_mtx.lock();
						for(int i = 0; i < (int)alignment_block.names.size(); ++i) {
							if(params.is_sam) output_fa2sam(alignment_block.names[i], local_bed.chr, local_bed.start, local_bed.end, alignment_block.seqs[i], params.read_group, -1, alignment_block.names.size(), alignment_block.statuses[i].spanning_l, alignment_block.statuses[i].spanning_r, -1, -1.0);
							else std::cout << '>' << region_str << ' ' << alignment_block.names[i] << ' ' << alignment_block.statuses[i].is_spanning() << '\n' << alignment_block.seqs[i] << '\n';
						}
						std_out_mtx.unlock();
					}
					else{
						std::vector<int> spannable_indeces;
						for(int j = 0; j < (int)alignment_block.statuses.size(); ++j) if(alignment_block.statuses[j].is_spanning()) spannable_indeces.emplace_back(j);
						if(!spannable_indeces.empty()) {
							//for(int j = 0; j < (int)alignment_block.names.size(); ++j) std::cout << '>' << alignment_block.names[j] << ' ' << alignment_block.statuses[j].is_spanning() << '\n' << alignment_block.seqs[j] << '\n';
							DistMatrix distmatrix(spannable_indeces.size());
							std::vector<int> labels;
							int initial_clusters;
							otter_hclust(params.max_alleles, params.bandwidth, params.max_error, params.min_cov_fraction, spannable_indeces, distmatrix, aligner, alignment_block, labels, initial_clusters);
							//for(int j = 0; j < (int)alignment_block.names.size(); ++j) std::cout << j << '\t' << labels[j] << '\n';
							std::vector<std::string> consensus_seqs;
							std::vector<std::vector<double>> ses;
							otter_rapid_consensus(spannable_indeces, labels, distmatrix, aligner2, alignment_block, consensus_seqs, ses);
							otter_nonspanning_assigment(params.min_sim, params.max_error, alignment_block, aligner, labels);
							std_out_mtx.lock();
							//for(int j = 0; j < (int)alignment_block.names.size(); ++j) std::cout << j << '\t' << labels[j] << '\t' << alignment_block.statuses[j].spanning_l << '\t' << alignment_block.statuses[j].spanning_r << '\t' << alignment_block.statuses[j].alignment_coords.first << '\t' << alignment_block.statuses[j].alignment_coords.second << '\n';
							for(int j = 0; j < (int)consensus_seqs.size(); ++j){
								int spanning_cov = 0, cov = 0;
								for(int k = 0; k < (int)spannable_indeces.size();++k) if(labels[spannable_indeces[k]] == j) ++spanning_cov;
								for(const auto& l : labels) if(l == j) ++cov;
								double se; 
								if(spanning_cov == 1) se = -1;
								else if(spanning_cov > 2) se = compute_se(ses[j]);
								else {
									std::vector<int> local_cluster_indeces;
									for(int k = 0; k < (int)spannable_indeces.size();++k) if(labels[spannable_indeces[k]] == j) local_cluster_indeces.emplace_back(k);
									se = distmatrix.get_dist(local_cluster_indeces[0], local_cluster_indeces[1]);
								}
								if(params.is_sam) {
									std::string assembly_name = region_str + "_" + std::to_string(j);
									output_fa2sam(assembly_name, local_bed.chr, local_bed.start, local_bed.end, consensus_seqs[j], params.read_group, cov, alignment_block.names.size(), -1, -1, initial_clusters, se);
								}
								else{
									std::cout << '>' << region_str << ' ' << cov << ' ' << alignment_block.names.size() << ' ' << initial_clusters << ' ' << se << '\n';
									std::cout << consensus_seqs[j] << '\n';
								}
							}
							std_out_mtx.unlock();
						}
					}
				}
			}
			bam_inst.destroy();
			if(!reference.empty()) faidx_inst.destroy();
	}).wait();
}


void assemble(const std::vector<std::string>& bams, const std::string& bed, const std::string& reference, const bool& reads_only, const OtterOpts& params)
{
 	BS::thread_pool pool(params.threads);

 	//load bed file
 	BedMap map_beds;
 	parse_bed_file(bed, map_beds);

 	std::vector<BED> bed_regions;
 	for(const auto& chr : map_beds) for(const auto& bed : chr.second) bed_regions.emplace_back(bed);
 	if(params.is_sam){
		BamInstance bam_inst;
		bam_inst.init(bams.front(), true);
		for(int i = 0; i < bam_inst.header->n_targets; ++i){
			std::cout << "@SQ\tSN:" << bam_inst.header->target_name[i] << "\tLN:" << bam_inst.header->target_len[i] << '\n';
		}
		if(!params.read_group.empty()) std::cout << "@RG\tID:" << params.read_group << '\n';
		bam_inst.destroy();
	}
 	for(const auto& bam : bams) general_process(params, bam, bed_regions, reference, reads_only, pool);
}