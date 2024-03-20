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

void general_process(const OtterOpts& params, const std::string& bam, const std::vector<BED>& bed_regions, const std::string& reference, const std::string& read_group, bool& reads_only, BS::thread_pool& pool)
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
							if(params.is_sam) output_fa2sam(alignment_block.names[i], local_bed.chr, local_bed.start, local_bed.end, alignment_block.seqs[i], read_group, -1, alignment_block.names.size(), alignment_block.statuses[i].spanning_l, alignment_block.statuses[i].spanning_r, -1, -1.0);
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
							otter_hclust(params.max_alleles, params.bandwidth, params.max_error, params.min_cov_fraction, params.min_cov_fraction2_l, params.min_cov_fraction2_f, spannable_indeces, distmatrix, aligner, alignment_block, labels, initial_clusters);
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
									output_fa2sam(assembly_name, local_bed.chr, local_bed.start, local_bed.end, consensus_seqs[j], read_group, cov, alignment_block.names.size(), -1, -1, initial_clusters, se);
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

void construct_bed_interval_tree(const uint32_t& offset, const std::vector<BED>& bed_regions, IntervalTree<int, int>& bed_tree)
{
    std::vector<Interval<int, int>> bed_intervals;
	for(int i = 0; i < (int)bed_regions.size(); ++i) {
		int start = bed_regions[i].start - offset;
		int end = bed_regions[i].end + offset;
		bed_intervals.emplace_back(start, end, i);
	}
	bed_tree = IntervalTree<int, int>(std::move(bed_intervals));
	std::cerr << '(' << antimestamp() << "): Constructed interval tree for " << bed_regions.size() << " target regions\n";
}

void wga_genotyper(const OtterOpts& params, const std::string& bam, const std::vector<BED>& bed_regions, BS::thread_pool& pool)
{
	IntervalTree<int, int> bed_tree;
	construct_bed_interval_tree(params.offset, bed_regions, bed_tree);

	std::vector<std::string> ref_chrms;
	{
		BamInstance bam_inst;
		bam_inst.init(bam, false);
		for(int i = 0; i < bam_inst.header->n_targets; ++i){
			std::string chr = bam_inst.header->target_name[i];
			uint32_t chr_l = bam_inst.header->target_len[i];
			ref_chrms.emplace_back(chr + ":1-" + std::to_string(chr_l));
		}
		bam_inst.destroy();
	}

	std::cerr << "(" << antimestamp() << "): Parallelising across " << ref_chrms.size() << " contigs" << std::endl;
	std::mutex std_out_mtx;

	pool.parallelize_loop(0, ref_chrms.size(), [&pool, &std_out_mtx, &bam, &bed_regions, &bed_tree, &ref_chrms, &params](const int a, const int b){
		BamInstance bam_inst;
		bam_inst.init(bam, true);
		if(!bam_inst.index) std::cerr << "(" << antimestamp() << "): WARNING: index not found. Skipping " << bam << std::endl;
		else{
			for(int chr_i = a; chr_i < b; ++chr_i){
				hts_itr_t *iter = sam_itr_querys(bam_inst.idx, bam_inst.header, ref_chrms[chr_i].c_str());
				if(iter == nullptr) std::cerr << "(" << antimestamp() << "): WARNING: query failed at contig " <<  ref_chrms[chr_i] << std::endl;
				else{
					while(sam_itr_next(bam_inst.fp, iter, bam_inst.read) > 0){
						uint32_t ref_end_pos = bam_endpos(bam_inst.read);
						auto bed_overlaps = bed_tree.findOverlapping(bam_inst.read->core.pos, ref_end_pos);
						auto it = bed_overlaps.begin();
						while(it != bed_overlaps.end()){
							if(bed_regions[it->value].chr == bam_inst.header->target_name[chr_i]) ++it;
							else it = bed_overlaps.erase(it);
						}
						if(!bed_overlaps.empty()){
							std::string name = (char*)bam_inst.read->data;
							//std::cout <<"#contig:" << name << '\t' << bam_inst.read->core.pos << '-' << ref_end_pos << '\n';
							std::vector<std::pair<int,int>> ref_positions;
							std::vector<OpInterval> query_positions;
							get_op_intervals(bam_inst.read, ref_positions, query_positions);
							if(ref_positions.size() != query_positions.size()){
								std::cerr << antimestamp() << "): Unexpected number of ref and query OP-intervals: " << ref_positions.size() << " vs " << query_positions.size() << '\n';
								exit(1);
							}
							std::vector<Interval<int, int>> op_intervals;
							for(int i = 0; i < (int)query_positions.size();++i) op_intervals.emplace_back(ref_positions[i].first, ref_positions[i].second, i);
							IntervalTree<int, int> op_tree;
							op_tree = IntervalTree<int, int>(std::move(op_intervals));
							for(const auto& overlap : bed_overlaps) {
								const auto& local_bed = bed_regions[overlap.value];
								//std::cout << "#" << local_bed.toBEDstring() << std::endl;
								auto bed_op_overlaps = op_tree.findOverlapping(overlap.start, overlap.stop);
								std::sort(bed_op_overlaps.begin(), bed_op_overlaps.end(), [](auto const& x, auto const& y){
									if(x.start == y.start) return x.stop < y.stop; 
									else return x.start < y.start;
								});
								bool clipped_l = false, clipped_r = false;
								int query_start = 0, query_end = 0;
								for(int i = 0; i < (int)bed_op_overlaps.size(); ++i){
									const auto& op_ref = bed_op_overlaps[i];
									const auto& op_query = query_positions[op_ref.value];
									//std::cout << local_bed.chr << ':' << overlap.start << '-' << overlap.stop << '\t' << op_ref.start << '-' << op_ref.stop << ':' << op_ref.value << '\t' << op_query.start << '-' << op_query.end << ':' << op_query.op << '\t' << query_start << '-' << query_end << '\n';
									if(op_query.op == BAM_CSOFT_CLIP || op_query.op == BAM_CHARD_CLIP) {
										if(i == 0){
											clipped_l = true;
											query_start = op_query.end;
										}
										else {
											clipped_r = true;
											query_end = op_query.start;
										}
									}
									else{
										if(i == 0) {
											if(op_query.op == BAM_CDEL){
												if(op_ref.start <= overlap.start && op_ref.stop >= overlap.stop) break;
												else query_start = op_query.start;
											}
											else query_start = op_query.start + (overlap.start - op_ref.start);
										}
										if(i + 1 == (int)bed_op_overlaps.size()) {
											if(op_query.op == BAM_CDEL) query_end = op_query.end;
											else query_end = op_query.end - (op_ref.stop - overlap.stop);
										}
									}
								}

								//std::cout << bam_inst.read->core.l_qseq << '\t' << (query_end - query_start) << '\t' << query_start << '-' << query_end << std::endl;
								uint8_t *q = bam_get_seq(bam_inst.read);
								std::string seq(query_end - query_start == 0 ? 1 : query_end - query_start, 'N');
								for(int i = query_start; i < query_end; i++) seq[i - query_start] = seq_nt16_str[bam_seqi(q, i)];

								//std::cout << overlap.start << '-' << overlap.stop << '\n';
								//std::cout << local_bed.toBEDstring() << '\t' << name << ':' << query_start << '-' << query_end << '\n';
								//std::cout << seq << '\n';
								std_out_mtx.lock();
								if(params.is_sam) output_fa2sam(name, local_bed.chr, local_bed.start, local_bed.end, seq, params.read_group, -1, -1, int(!clipped_l), int(!clipped_r), -1, -1);
								else{
									std::cout << '>' << local_bed.toBEDstring() << '\t' << name << '\t' << int(!(clipped_l || clipped_r)) << std::endl;
									std::cout << seq << std::endl;
								}
								std_out_mtx.unlock();
							
							}
						}
					}
				}
				hts_itr_destroy(iter);
			}
		}
		bam_inst.destroy();
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
		bam_inst.destroy();
	}

	{
		std::vector<std::string> read_groups_tmp = read_groups;
		sort(read_groups_tmp.begin(), read_groups_tmp.end());
		for(int i = 1; i < (int)read_groups_tmp.size(); ++i){
			if(read_groups_tmp[i] == read_groups_tmp[i-1]) std::cerr << '(' << antimestamp() << "): WARNING: non-unique sample identifier across multiple BAM-files: " << read_groups_tmp[i] << std::endl;
		}
	}

 	for(int i = 0; i < (int)read_groups.size(); ++i){
 		if(params.is_wga) wga_genotyper(params, bams[i], bed_regions, pool);
 		else general_process(params, bams[i], bed_regions, reference, read_groups[i], reads_only, pool);
 	}
}