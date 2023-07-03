#include "antimestamp.hpp"
#include "bam_db.hpp"
#include "anbamfilehelper.hpp"
#include "anbed.hpp"
#include "otter_opts.hpp"
#include "pairwise_alignment.hpp"
#include "bindings/cpp/WFAligner.hpp"
#include <iostream>
#include <vector>


void pairwise_process(const OtterOpts& params, const std::vector<BED>& regions, const std::string& bam, const std::vector<std::string>& index2sample, const std::map<std::string,int>& sample2index)
{
	BS::thread_pool pool(params.threads);
	std::mutex stdout_mtx;
	pool.parallelize_loop(0, regions.size(),
		[&pool, &stdout_mtx, &bam, &regions, &index2sample, &sample2index](const int a, const int b){
			BamInstance bam_inst;
			bam_inst.init(bam, true);
			wfa::WFAlignerEdit aligner(wfa::WFAligner::Score, wfa::WFAligner::MemoryMed);
			for(int i = a; i < b; ++i) {
				const BED& region = regions[i];
				const std::string region_str = region.chr + ':' + std::to_string(region.start) + '-' + std::to_string(region.end);
				hts_itr_t *iter = sam_itr_querys(bam_inst.idx, bam_inst.header, region_str.c_str());
				if(iter == nullptr) std::cerr << "(" << antimestamp() << "): WARNING: query failed at region " << region_str << std::endl;
				else{
					std::vector<int> allele_sample_indeces;
					std::vector<std::string> alleles;
					while(sam_itr_next(bam_inst.fp, iter, bam_inst.read) > 0) parse_bam_allele(region_str, 1, sample2index, bam_inst.read, alleles, allele_sample_indeces);
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

					}
				}
			}
			bam_inst.destroy();
	}).wait();

}

void genotype(const std::string& bam, const std::string& bed, const OtterOpts& params)
{
	std::vector<BED> regions;
	parse_bed_file(bed, regions);

	std::vector<std::string> sample_index;
	std::map<std::string, int> sample2index;
	index_read_groups(bam, sample2index, sample_index);

	std::cerr << '(' << antimestamp() << "): Found " << sample_index.size() << " samples (read-group tags)\n";

	pairwise_process(params, regions, bam, sample_index, sample2index);
}