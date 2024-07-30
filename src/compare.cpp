#include "bindings/cpp/WFAligner.hpp"
#include "antimestamp.hpp"
#include "anbamdb.hpp"
#include "anbamfilehelper.hpp"
#include "anseqs.hpp"
#include "anbed.hpp"
#include "BS_thread_pool.hpp"
#include "otter_opts.hpp"
#include "compare.hpp"
#include "sam.h"
#include <map>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

class DistCompare{
	public:
		int i;
		int j;
		double edit;
		double ops;
		DistCompare(int _i, int _j, double _e, double _o): i(_i), j(_j), edit(_e), ops(_o){};
};

void local_parse_analleles(const BamInstance& bam_inst, const BED& bed, const std::map<std::string, int>& sample2index, std::vector<ANALLELE>& anallele_block, std::vector<int>& allele_sample_indeces, std::vector<int>& spannings)
{
	hts_itr_t *iter = bam_itr_querys(bam_inst.idx, bam_inst.header, bed.toScString().c_str());
	if(iter == nullptr) std::cerr << "(" << antimestamp() << "): WARNING: query failed at region " <<  bed.toScString() << std::endl;
	else{
		while(bam_itr_next(bam_inst.fp, iter, bam_inst.read) > 0){
			std::string name = (char*)bam_inst.read->data;
			if(name.substr(0, bed.chr.size()) == bed.chr){
				char spanning = 'u';
				auto aux_ptr = bam_aux_get(bam_inst.read, sp_tag.c_str());
				if(aux_ptr != NULL) spanning = bam_aux2A(aux_ptr);
				ANALLELE anallele;
				parse_anallele(bed.toScString(), sample2index, bam_inst.read, anallele_block, allele_sample_indeces);
				if(spanning == 'u') spannings.emplace_back(-1);
				else if(spanning == 'b') spannings.emplace_back(0);
				else if(spanning == 'l') spannings.emplace_back(1);
				else if(spanning == 'r') spannings.emplace_back(2);
				else if(spanning == 'n') spannings.emplace_back(3);
			}
		}
	}
	bam_itr_destroy(iter);
}

void get_distances(wfa::WFAligner& aligner, std::vector<ANALLELE>& subjs, std::vector<ANALLELE>& querys, std::vector<DistCompare>& distances)
{
	for(uint32_t i = 0; i < subjs.size(); ++i){
		auto& subj = subjs[i].seq;
		for(uint32_t j = 0; j < querys.size(); ++j){
			auto& query = querys[j].seq;
			if(subj == query || (subj == "N" && query == "NDNNN") || (query == "N" && subj == "NDNNN")) distances.emplace_back(i, j, 0, query.size());
			else if(subj == "N" || query == "N" || subj == "NDNNN" || query == "NDNNN") distances.emplace_back(i, j, query.size() - 1, query.size());
			else {
				if(subj.size() > query.size()) aligner.alignEnd2End(subj, query); else aligner.alignEnd2End(query, subj);
				int edit = aligner.getAlignmentScore();
				int ops = aligner.getAlignmentCigar().size();
				distances.emplace_back(i, j, edit, ops);
			}
		}
	}
}

void compare(const OtterOpts& params, const std::string& bed_file, const std::string& reference, const std::string& target)
{
	std::vector<BED> regions;
	parse_bed_file(bed_file, regions);

	BamInstance bam_ref, bam_target;
	bam_ref.init(reference, true);
	bam_target.init(target, true);

	std::map<std::string,int> sample2index;
	{
		SampleIndex si;
		si.init(reference);
		sample2index[si.index2sample.front()] = 0;

	}
	{
		SampleIndex sit;
		sit.init(target);
		sample2index[sit.index2sample.front()] = 1;
	}

	BS::thread_pool pool(params.threads);
	std::mutex stdout_mtx;

	pool.parallelize_loop(0, regions.size(),
		[&pool, &stdout_mtx, &params, &regions, &sample2index, &bam_ref, &bam_target](const int a, const int b){
			wfa::WFAlignerEdit aligner(wfa::WFAligner::Alignment, wfa::WFAligner::MemoryMed);
			for(int region_i = a; region_i < b; ++region_i) {
				const BED& region = regions[region_i];
				const std::string region_str = region.toScString();
				//std::cout << region_str << '\n';
				std::vector<ANALLELE> reference_alleles;
				std::vector<int> reference_spannings;
				std::vector<ANALLELE> query_alleles;
				std::vector<int> allele_sample_indeces;
				local_parse_analleles(bam_ref, region, sample2index, reference_alleles, allele_sample_indeces, reference_spannings);
				parse_analleles(params, bam_target, region, sample2index, query_alleles, allele_sample_indeces);
				if(query_alleles.size() == 1) query_alleles.emplace_back(query_alleles.front());
				if(reference_alleles.size() > 2) std::cerr << "(" << antimestamp() << "): WARNING: skipping region due to multiple expected alignments (>2) for region: " << region.toScString() <<std::endl;
				else if(reference_alleles.size() == 1) std::cerr << "(" << antimestamp() << "): WARNING: skipping region due to single expected alignment for region: " << region.toScString() <<std::endl;
				else if(reference_alleles.size() == 0) std::cerr << "(" << antimestamp() << "): WARNING: skipping region due no expected alignments for region: " << region.toScString() <<std::endl;
				else if(query_alleles.size() == 0) std::cerr << "(" << antimestamp() << "): WARNING: skipping region due no query alleles for region: " << region.toScString() <<std::endl;
				else {
					std::vector<DistCompare> dist_edges;
					get_distances(aligner, reference_alleles, query_alleles, dist_edges);
					std::sort(dist_edges.begin(), dist_edges.end(), [](const auto& x, const auto& y){
						if(x.edit == y.edit) return x.ops < y.ops; 
						else return x.edit < y.edit;
					});
					/**
					std::cout << region_str << '\t' << reference_alleles.size() << '\t' << query_alleles.size() << '\t' << allele_sample_indeces.size();
					for(const auto& e : dist_edges) std::cout << '\t' << e.i << ',' << e.j << ',' << e.edit << ',' << e.ops;
					std::cout << '\n';
					*/
					
					uint32_t edge_1_j = 1;
					const auto& edge_0 = dist_edges.front();
					for(; edge_1_j < dist_edges.size(); ++edge_1_j){
						const auto& edge_1 = dist_edges[edge_1_j];
						//if(reference_alleles.size() == 1 && edge_1.j != edge_0.j) break;
						if(edge_1.i != edge_0.i && edge_1.j != edge_0.j) break;
					}
					const auto& edge_1 = dist_edges[edge_1_j];
					//std::cout << "min:\t" << edge_0.i << ',' << edge_0.j << ',' << edge_0.edit << ',' << edge_0.ops << '\t' << edge_1.i << ',' << edge_1.j << ',' << edge_1.edit << ',' << edge_1.ops << '\n';
					std::vector<int> min_edges_indeces{0, (int)edge_1_j};
					for(const auto& i : min_edges_indeces){
						auto& min_edge = dist_edges[i];
						/**
						if(reference_alleles.size() == 1 && i > 0) {
							min_edge.edit = reference_alleles.front().seq.size() - 1;
							min_edge.ops = reference_alleles.front().seq.size();
						}
						*/
						std::cout << region_str << '\t' << reference_alleles[min_edge.i].seq.size() << '\t' << query_alleles[min_edge.j].seq.size() << '\t' << reference_spannings[min_edge.i] << '\t' << min_edge.edit << '\t' << min_edge.ops << '\n';
					}
				}
				
			}
	}).wait();
	bam_ref.destroy();
	bam_target.destroy();
}