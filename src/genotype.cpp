#include "antimestamp.hpp"
#include "bam_db.hpp"
#include "anbamfilehelper.hpp"
#include "anbed.hpp"
#include "otter_opts.hpp"
#include "pairwise_alignment.hpp"
#include "bindings/cpp/WFAligner.hpp"
#include "fastcluster.h"
#include <iostream>
#include <vector>


void pairwise_process(const OtterOpts& params, const int& ac_mincov, const int& tc_mincov, const bool& is_summary, const std::vector<BED>& regions, const std::string& bam, const std::vector<std::string>& index2sample, const std::map<std::string,int>& sample2index)
{
	BS::thread_pool pool(params.threads);
	std::mutex stdout_mtx;
	pool.parallelize_loop(0, regions.size(),
		[&pool, &stdout_mtx, &params, &ac_mincov, &tc_mincov,&is_summary, &bam, &regions, &index2sample, &sample2index](const int a, const int b){
			BamInstance bam_inst;
			bam_inst.init(bam, true);
			wfa::WFAlignerEdit aligner(wfa::WFAligner::Score, wfa::WFAligner::MemoryMed);
			for(int region_i = a; region_i < b; ++region_i) {
				const BED& region = regions[region_i];
				const std::string region_str = region.chr + ':' + std::to_string(region.start) + '-' + std::to_string(region.end);
				hts_itr_t *iter = sam_itr_querys(bam_inst.idx, bam_inst.header, region_str.c_str());
				if(iter == nullptr) std::cerr << "(" << antimestamp() << "): WARNING: query failed at region " << region_str << std::endl;
				else{
					std::vector<int> allele_sample_indeces;
					std::vector<std::string> alleles;
					while(sam_itr_next(bam_inst.fp, iter, bam_inst.read) > 0) parse_bam_allele(region_str, ac_mincov, tc_mincov, sample2index, bam_inst.read, alleles, allele_sample_indeces);
					if(alleles.size() > 1 && is_multi_sample(allele_sample_indeces)){
						sort_bam_alleles(allele_sample_indeces, alleles);;
						std::vector<std::pair<int, std::pair<int,int>>> sample2intervals;
						set_sample_intervals(allele_sample_indeces, sample2intervals);
						DistMatrix matrix(alleles.size());
						for(int i = 0; i < (int)alleles.size(); ++i){
							int i_l = alleles[i].size();
							for(int j = i + 1; j < (int)alleles.size(); ++j){
								int j_l = alleles[j].size();
								double max_l;
								if(i_l > j_l){
									max_l = i_l;
									aligner.alignEnd2End(alleles[j], alleles[i]);
								}
								else{
									max_l = j_l;
									aligner.alignEnd2End(alleles[i], alleles[j]);
								}
								matrix.set_dist(i, j, aligner.getAlignmentScore() / max_l);
							}
						}
						int* merge = new int[2*(alleles.size()-1)];
	    				double* height = new double[alleles.size()-1];
	    				hclust_fast(alleles.size(), matrix.values.data(), HCLUST_METHOD_AVERAGE, merge, height);
	    				int* labels = new int[alleles.size()];
		    			cutree_cdist(alleles.size(), merge, height, params.max_error, labels);
		    			if(is_summary){
		    				int max_label = 0;
		    				for(int i = 0; i < (int)alleles.size(); ++i) if(labels[i] > max_label) max_label = labels[i];
		    				++max_label;
		    				for(int label = 0; label < max_label; ++label){
		    					int n = 0;
		    					double size = 0;
		    					for(int i = 0; i < (int)alleles.size(); ++i) {
		    						if(labels[i] == label){
		    							++n;
		    							size += alleles[i].size();
		    						}
		    					}
		    					std::cout << region_str << '\t' << label << '\t' << n << '\t' << (size/n) << '\n';
		    				}
		    			}
		    			else{
		    				for(int i = 0; i < (int)sample2intervals.size(); ++i){
			    				int a1 = *(labels+sample2intervals[i].second.first);
			    				int a2 = *(labels+sample2intervals[i].second.second);
			    				bool is_a1_min = a1 < a2;
			    				std::cout << index2sample[sample2intervals[i].first] << '\t' << region_str << '\t';
			    				if(is_a1_min) std::cout << a1 << '/' << a2; else std::cout << a2 << '/' << a1;
			    				std::cout << '\n';
			    			}
		    			}
					}
				}
			}
			bam_inst.destroy();
	}).wait();

}

void genotype(const std::string& bam, const std::string& bed, const OtterOpts& params, const int& ac_mincov, const int& tc_mincov, const bool& is_summary)
{
	std::vector<BED> regions;
	parse_bed_file(bed, regions);

	std::vector<std::string> sample_index;
	std::map<std::string, int> sample2index;
	index_read_groups(bam, sample2index, sample_index);

	std::cerr << '(' << antimestamp() << "): Found " << sample_index.size() << " samples (read-group tags)\n";

	pairwise_process(params, ac_mincov, tc_mincov, is_summary, regions, bam, sample_index, sample2index);
}