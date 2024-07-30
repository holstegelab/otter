#include "wgat.hpp"
#include "antimestamp.hpp"
#include "anbamfilehelper.hpp"
#include "anbed.hpp"
#include "otter_opts.hpp"
#include "anfahelper.hpp"
#include "anseqs.hpp"
#include "BS_thread_pool.hpp"
#include "faidx.h"
#include "interval_tree.h"
#include "opinterval.hpp"
#include <iostream>
#include <vector>
#include <thread>
#include <map>
#include <mutex>
#include <cmath>

void construct_bed_interval_tree(const uint32_t& offset_l, const uint32_t& offset_r, const std::vector<BED>& bed_regions, IntervalTree<int, int>& bed_tree)
{
    std::vector<Interval<int, int>> bed_intervals;
	for(int i = 0; i < (int)bed_regions.size(); ++i) {
		int start = bed_regions[i].start - offset_l;
		int end = bed_regions[i].end + offset_r;
		bed_intervals.emplace_back(start, end, i);
	}
	bed_tree = IntervalTree<int, int>(std::move(bed_intervals));
	std::cerr << '(' << antimestamp() << "): Constructed interval tree for " << bed_regions.size() << " target regions\n";
}

void wga_bam_genotyper_process(const OtterOpts& params, const std::vector<BED>& bed_regions, const IntervalTree<int, int>& bed_tree, const std::vector<std::string>& ref_chrms, const int& a, const int& b, std::mutex& std_out_mtx, BS::thread_pool& pool, BamInstance& bam_inst)
{
	for(int chr_i = a; chr_i < b; ++chr_i){
		hts_itr_t *iter = bam_itr_querys(bam_inst.idx, bam_inst.header, ref_chrms[chr_i].c_str());
		if(iter == nullptr) std::cerr << "(" << antimestamp() << "): WARNING: query failed at region " <<  ref_chrms[chr_i] << std::endl;
		else{
			int alignment_index = 0;
			while(bam_itr_next(bam_inst.fp, iter, bam_inst.read) > 0){
				if(bam_inst.read->core.l_qseq > 0){
					uint32_t ref_end_pos = bam_inst.read->core.pos + bam_cigar2rlen(bam_inst.read->core.n_cigar, bam_get_cigar(bam_inst.read));
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
							//std::cout << "overlaps: " << bed_op_overlaps.size() << '\n';
							for(int i = 0; i < (int)bed_op_overlaps.size(); ++i){
								const auto& op_ref = bed_op_overlaps[i];
								const auto& op_query = query_positions[op_ref.value];
								//std::cout << i << ' ' << local_bed.chr << ':' << overlap.start << '-' << overlap.stop << '\t' << op_ref.start << '-' << op_ref.stop << ':' << op_ref.value << '\t' << op_query.start << '-' << op_query.end << ':' << op_query.op << '\t' << bam_inst.read->core.l_qseq << '\t' << query_start << '-' << query_end << '\n';
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
							//std::cout << "final " << local_bed.chr << ':' << overlap.start << '-' << overlap.stop << ' ' <<  query_start << '-' << query_end << ' ' << clipped_l << ' ' << clipped_r << '\n';
							if(clipped_l || clipped_r) std::cerr << antimestamp() << "): [WARNING] skipping non-spanning whole-genome alignment at " << local_bed.toScString() << " for: " << name << '\n';
							else{
								//std::cout << bam_inst.read->core.l_qseq << '\t' << (query_end - query_start) << '\t' << query_start << '-' << query_end << std::endl;
								uint8_t *q = bam_get_seq(bam_inst.read);
								std::string seq(query_end - query_start == 0 ? 1 : query_end - query_start, 'N');
								for(int i = query_start; i < query_end; i++) seq[i - query_start] = seq_nt16_str[bam_seqi(q, i)];
								ANALLELE allele(seq);
								//std::cout << overlap.start << '-' << overlap.stop << '\n';
								//std::cout << local_bed.toBEDstring() << '\t' << name << ':' << query_start << '-' << query_end << '\n';
								//std::cout << seq << '\n';
								std_out_mtx.lock();
								if(params.is_fa) allele.stdout_fa(params.read_group, name + '#' + local_bed.toScString() + '#' + std::to_string(alignment_index), true, int(!clipped_l), int(!clipped_r));
								else allele.stdout_sam(name + '#' + local_bed.toScString() + "_" + std::to_string(alignment_index), local_bed.chr, local_bed.start, local_bed.end, params.read_group, true, int(!clipped_l), int(!clipped_r));
								std_out_mtx.unlock();
							}
						}
					}
					++alignment_index;
				}
			}
		}
		bam_itr_destroy(iter);
	}
}


void wga_bam_genotyper(const OtterOpts& params, const std::string& bam, const std::vector<BED>& bed_regions)
{
	BS::thread_pool pool(params.threads);
	IntervalTree<int, int> bed_tree;
	construct_bed_interval_tree(params.offset_l, params.offset_r, bed_regions, bed_tree);

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

	std::mutex std_out_mtx;

	std::cerr << "(" << antimestamp() << "): Parallelising across " << ref_chrms.size() << " contigs" << std::endl;
	pool.parallelize_loop(0, ref_chrms.size(), [&pool, &std_out_mtx, &bam, &bed_regions, &bed_tree, &ref_chrms, &params](const int a, const int b){
		BamInstance bam_inst;
		bam_inst.init(bam, true);
		if(!bam_inst.index) {
			std::cerr << "(" << antimestamp() << "): [ERROR] index not found " << bam << std::endl;
			exit(1);
		}
		else wga_bam_genotyper_process(params, bed_regions, bed_tree, ref_chrms, a, b, std_out_mtx, pool, bam_inst);
		bam_inst.destroy();
	}).wait();
}


void wgat(const OtterOpts& params, const std::string& input, const std::string& bed_file){

 	std::vector<BED> bed_regions;
 	parse_bed_file(bed_file, bed_regions);
 	if(input.substr(input.find_last_of(".") + 1) == "bam") {
 		if(!params.is_fa){
 			BamInstance bam_inst;
			bam_inst.init(input, true);
			for(int i = 0; i < bam_inst.header->n_targets; ++i){
				std::cout << "@SQ\tSN:" << bam_inst.header->target_name[i] << "\tLN:" << bam_inst.header->target_len[i] << '\n';
			}
			std::cout << "@RG\tID:" << params.read_group << '\n';
			std::cout << "@PG\tID:otter\tOF:" << params.offset_l << ',' << params.offset_r << '\n';
			bam_inst.destroy();
 		}
 		wga_bam_genotyper(params, input, bed_regions);
 	} 

}