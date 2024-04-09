#include "assemble.hpp"
#include "antimestamp.hpp"
#include "anbamfilehelper.hpp"
#include "parse_bam_alignments.hpp"
#include "bam_db.hpp"
#include "anbed.hpp"
#include "otter_opts.hpp"
#include "pairwise_alignment.hpp"
#include "fasta_helper.hpp"
#include "formatter.hpp"
#include "BS_thread_pool.hpp"
#include <htslib/faidx.h>
#include "bindings/cpp/WFAligner.hpp"
#include "interval_tree.h"
#include "opinterval.hpp"
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

void general_process(const OtterOpts& params, const std::string& bam, const std::vector<BED>& bed_regions, const std::string& reference, const std::string& read_group, bool& reads_only, BS::thread_pool& pool, std::mutex& std_out_mtx, BamInstance& bam_inst, FaidxInstance& faidx_inst, wfa::WFAlignerEdit& aligner, wfa::WFAlignerGapAffine& aligner2, const int& a, const int& b)
{
	for(int bed_i = a; bed_i < b; ++bed_i) {
		AlignmentBlock alignment_block;
		const BED& local_bed = bed_regions[bed_i];
		BED mod_bed = local_bed;
		mod_bed.start -= (uint32_t)params.offset;
		mod_bed.end += (uint32_t)params.offset;
		parse_alignments(params, mod_bed, bam_inst, alignment_block);
		if(alignment_block.names.size() == 0) {
			std_out_mtx.lock();
			std::cerr << '(' << antimestamp() << "): WARNING: no reads found at region " << local_bed.toBEDstring() << '\n';
			std_out_mtx.unlock();
		}
		else if ((int)alignment_block.names.size() > params.max_cov) std::cerr << '(' << antimestamp() << "): WARNING: abnormal coverage for region " << local_bed.toBEDstring() << " (" << alignment_block.names.size() << "x)\n";
		else{
			std::string region_str = local_bed.toBEDstring();
			if(!reference.empty()) otter_realignment(local_bed.chr, (int)mod_bed.start, (int)mod_bed.end, params.flank, params.min_sim, faidx_inst, alignment_block, aligner2);
			if(reads_only){
				std_out_mtx.lock();
				for(int i = 0; i < (int)alignment_block.names.size(); ++i) {
					if(params.is_sam) output_fa2sam(alignment_block.names[i], local_bed.chr, local_bed.start, local_bed.end, alignment_block.seqs[i], read_group, -1, alignment_block.names.size(), alignment_block.statuses[i].spanning_l, alignment_block.statuses[i].spanning_r, -1, -1.0, alignment_block.hps[i].ps, alignment_block.hps[i].hp);
					else std::cout << '>' << region_str << ' ' << alignment_block.names[i] << ' ' << alignment_block.statuses[i].is_spanning() << '\n' << alignment_block.seqs[i] << '\n';
				}
				std_out_mtx.unlock();
			}
			else{
				std::vector<int> spannable_indeces;
				for(int j = 0; j < (int)alignment_block.statuses.size(); ++j) if(alignment_block.statuses[j].is_spanning()) spannable_indeces.emplace_back(j);
				if(spannable_indeces.empty()) {
					std_out_mtx.lock();
					std::cerr << '(' << antimestamp() << "): WARNING: spanning reads found " << local_bed.toBEDstring() << '\n';
					std_out_mtx.unlock();
				}
				else {
					//for(int j = 0; j < (int)alignment_block.names.size(); ++j) std::cout << '>' << alignment_block.names[j] << ' ' << alignment_block.statuses[j].is_spanning() << '\n' << alignment_block.seqs[j] << '\n';
					DistMatrix distmatrix(spannable_indeces.size());
					std::vector<int> labels;
					int initial_clusters;
					otter_hclust(params.ignore_haps, params.max_alleles, params.bandwidth, params.max_error, params.min_cov_fraction, params.min_cov_fraction2_l, params.min_cov_fraction2_f, spannable_indeces, distmatrix, aligner, alignment_block, labels, initial_clusters);
					//for(int j = 0; j < (int)alignment_block.names.size(); ++j) std::cout << j << '\t' << labels[j] << '\t' << alignment_block.hps[j] << '\n';
					if(!params.ignore_haps) {
						int max_label = 0;
						for(int l = 0; l < (int)labels.size(); ++l) if(labels[l] > max_label) max_label = labels[l];
						++max_label;
						for(int l = 0; l < (int)max_label; ++l){
							std::vector<int> cluster_indeces;
							for(int k = 0; k < (int)spannable_indeces.size(); ++k) if(labels[spannable_indeces[k]] == l) cluster_indeces.emplace_back(k);
							sort(cluster_indeces.begin(), cluster_indeces.end());
							otter_pairwise_dist(!params.ignore_haps, cluster_indeces, true, alignment_block, aligner, distmatrix);
						}
					}
					std::vector<std::string> consensus_seqs;
					std::vector<std::vector<double>> ses;
					otter_rapid_consensus(spannable_indeces, labels, distmatrix, aligner2, alignment_block, consensus_seqs, ses);
					otter_nonspanning_assigment(params.ignore_haps, params.min_sim, params.max_error, alignment_block, aligner, labels);
					//for(int j = 0; j < (int)alignment_block.names.size(); ++j) std::cout << j << '\t' << labels[j] << '\t' << alignment_block.statuses[j].spanning_l << '\t' << alignment_block.statuses[j].spanning_r << '\t' << alignment_block.statuses[j].alignment_coords.first << '\t' << alignment_block.statuses[j].alignment_coords.second << '\n';
					for(int j = 0; j < (int)consensus_seqs.size(); ++j){
						int spanning_cov = 0, cov = 0, local_hp = -1, local_ps = -1;
						bool conflict_hp = false, conflict_ps = false;
						for(int k = 0; k < (int)spannable_indeces.size();++k) {
							if(labels[spannable_indeces[k]] == j) {
								++spanning_cov;
							}
						}
						for(int k = 0; k < (int)labels.size(); ++k) {
							if(labels[k] == j && alignment_block.hps[k].hp >= 0) {
								if(local_hp >= 0 && alignment_block.hps[k].hp != local_hp) conflict_hp = true;
								if(local_ps >= 0 && alignment_block.hps[k].ps != local_ps) conflict_ps = true;
								local_hp = alignment_block.hps[k].hp;
								local_ps = alignment_block.hps[k].ps;
							}
						}
						if(conflict_ps){
							std_out_mtx.lock();
							std::cerr << '(' << antimestamp() << "): WARNING: conflicting PS-tag for reads in " << local_bed.toBEDstring() << '\n';
							std_out_mtx.unlock();
							local_ps = -1;
						}
						if(conflict_hp){
							std_out_mtx.lock();
							std::cerr << '(' << antimestamp() << "): WARNING: conflicting HP-tag for reads in " << local_bed.toBEDstring() << '\n';
							std_out_mtx.unlock();
							local_hp = -1;
						}
						if(params.max_alleles == 1) local_hp = -1;
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
							std_out_mtx.lock();
							output_fa2sam(assembly_name, local_bed.chr, local_bed.start, local_bed.end, consensus_seqs[j], read_group, cov, alignment_block.names.size(), -1, -1, initial_clusters, se, local_ps, local_hp);
							std_out_mtx.unlock();
						}
						else{
							std_out_mtx.lock();
							std::cout << '>' << region_str << ' ' << cov << ' ' << alignment_block.names.size() << ' ' << initial_clusters << ' ' << se << '\n';
							std::cout << consensus_seqs[j] << '\n';
							std_out_mtx.unlock();
						}
					}
				}
			}
		}
	}

}

void sing_sample_process(const OtterOpts& params, const std::string& bam, const std::vector<BED>& bed_regions, const std::string& reference, const std::string& read_group, bool& reads_only, BS::thread_pool& pool)
{
	std::cerr << '(' << antimestamp() << "): Processing " << bam << " (" << read_group << ')' << std::endl;
	std::mutex std_out_mtx;
	pool.parallelize_loop(0, bed_regions.size(),
		[&pool, &std_out_mtx, &params, &bam, &bed_regions, &reference, &reads_only, &read_group](const int a, const int b){
			BamInstance bam_inst;
			bam_inst.init(bam, true);
			wfa::WFAlignerEdit aligner(wfa::WFAligner::Score, wfa::WFAligner::MemoryMed);
			wfa::WFAlignerGapAffine aligner2(4,6,2,wfa::WFAligner::Alignment,wfa::WFAligner::MemoryMed);
			FaidxInstance faidx_inst;
			if(!reference.empty()) faidx_inst.init(reference);
			general_process(params, bam, bed_regions, reference, read_group, reads_only, pool, std_out_mtx, bam_inst, faidx_inst, aligner, aligner2, a, b);
			bam_inst.destroy();
			if(!reference.empty()) faidx_inst.destroy();
	}).wait();
}

void multi_sample_process(const OtterOpts& params, const std::vector<std::string>& bams, const std::vector<BED>& bed_regions, const std::string& reference, const std::vector<std::string>& read_groups, bool& reads_only, BS::thread_pool& pool)
{
	std::mutex std_out_mtx;
	pool.parallelize_loop(0, bams.size(), [&pool, &std_out_mtx, &params, &bams, &bed_regions, &reference, &reads_only, &read_groups](const int a, const int b){
		for(int bam_i = a; bam_i < b; ++bam_i) {
			std_out_mtx.lock();
			std::cerr << '(' << antimestamp() << "): Processing " << bams[bam_i] << " (" << read_groups[bam_i] << ')' << std::endl;
			std_out_mtx.unlock();
			BamInstance bam_inst;
			bam_inst.init(bams[bam_i], true);
			wfa::WFAlignerEdit aligner(wfa::WFAligner::Score, wfa::WFAligner::MemoryMed);
			wfa::WFAlignerGapAffine aligner2(4,6,2,wfa::WFAligner::Alignment,wfa::WFAligner::MemoryMed);
			FaidxInstance faidx_inst;
			if(!reference.empty()) faidx_inst.init(reference);
			general_process(params, bams[bam_i], bed_regions, reference, read_groups[bam_i], reads_only, pool, std_out_mtx, bam_inst, faidx_inst, aligner, aligner2, 0, bed_regions.size());
			bam_inst.destroy();
			if(!reference.empty()) faidx_inst.destroy();
		}
	}).wait();
}


void assemble(const std::vector<std::string>& bams, const std::string& bed, const std::string& reference, bool& reads_only, const OtterOpts& params)
{
 	BS::thread_pool pool(params.threads);

 	std::vector<BED> bed_regions;
 	parse_bed_file(bed, bed_regions);

 	std::vector<std::string> read_groups(bams.size());
 	if(params.is_sam){
		BamInstance bam_inst;
		bam_inst.init(bams.front(), true);
		for(int i = 0; i < bam_inst.header->n_targets; ++i){
			std::cout << "@SQ\tSN:" << bam_inst.header->target_name[i] << "\tLN:" << bam_inst.header->target_len[i] << '\n';
		}
		for(int i = 0; i < (int)read_groups.size(); ++i) {
			if(params.read_group.empty()) read_groups[i] = bams[i].substr(bams[i].find_last_of("/\\") + 1);
			else read_groups[i] = params.read_group;
		}
		for(const auto& rg : read_groups) std::cout << "@RG\tID:" << rg << '\n';
		output_preset_offset_tag(params.offset);
		bam_inst.destroy();
	}

	{
		std::vector<std::string> read_groups_tmp = read_groups;
		sort(read_groups_tmp.begin(), read_groups_tmp.end());
		for(int i = 1; i < (int)read_groups_tmp.size(); ++i){
			if(read_groups_tmp[i] == read_groups_tmp[i-1]) std::cerr << '(' << antimestamp() << "): WARNING: non-unique sample identifier across multiple BAM-files: " << read_groups_tmp[i] << std::endl;
		}
	}

	if(bams.size() == 1) sing_sample_process(params, bams[0], bed_regions, reference, read_groups[0], reads_only, pool);
	else multi_sample_process(params, bams, bed_regions, reference, read_groups, reads_only, pool);
}
