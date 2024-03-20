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

void gt_pairwise_alignment(const double& max_error, wfa::WFAlignerEdit& aligner, DistMatrix& matrix, std::vector<std::string>& alleles, std::vector<int>& labels)
{
	for(int i = 0; i < (int)alleles.size(); ++i){
		int i_l = alleles[i].size();
		for(int j = i + 1; j < (int)alleles.size(); ++j){
			int j_l = alleles[j].size();
			if(i_l == j_l && alleles[i] == alleles[j]) matrix.set_dist(i, j, 0.0);
			else{
				int dist = i_l > j_l ? i_l - j_l : j_l - i_l;
				int max_l = i_l > j_l ? i_l : j_l;
				if(((double)dist / max_l) > max_error) matrix.set_dist(i, j, max_error + 0.01);
				else{
					if(i_l > j_l)aligner.alignEnd2End(alleles[j], alleles[i]);
					else aligner.alignEnd2End(alleles[i], alleles[j]);
					matrix.set_dist(i, j, (double)aligner.getAlignmentScore() / max_l);
				}
			}
		}
	}
	int* merge = new int[2*(alleles.size()-1)];
	double* height = new double[alleles.size()-1];
	hclust_fast(alleles.size(), matrix.values.data(), HCLUST_METHOD_AVERAGE, merge, height);
	cutree_cdist(alleles.size(), merge, height, max_error, labels.data());
	delete[] height;
	delete[] merge;
}

void recompute_pairwise_alignment(wfa::WFAlignerEdit& aligner, DistMatrix& matrix, std::vector<std::string>& alleles, std::vector<int>& reps, bool& length_dist)
{
	for(int i = 0; i < (int)reps.size(); ++i){
		int i_l = alleles[reps[i]].size();
		for(int j = i + 1; j < (int)reps.size(); ++j){
			int j_l = alleles[reps[j]].size();
			bool is_i_max = i_l > j_l ;
			int max_l = is_i_max ? i_l : j_l;
			int dist_l = (is_i_max ? i_l - j_l : j_l - i_l);
			if(length_dist) matrix.set_dist(reps[i], reps[j], dist_l/double(max_l));
			else{
				if(i_l > j_l) aligner.alignEnd2End(alleles[reps[j]], alleles[reps[i]]);
				else aligner.alignEnd2End(alleles[reps[i]], alleles[reps[j]]);
				matrix.set_dist(reps[i], reps[j], (double)aligner.getAlignmentScore() / max_l);
			}
		}
	}
}

void get_cluster_size(std::vector<int>& labels, std::vector<int>& label_counts)
{
		for(int i = 0; i < (int)labels.size(); ++i) ++label_counts[labels[i]];
}

void get_summary(const std::string& region_str, std::mutex& stdout_mtx, std::vector<int>& labels, std::vector<std::string>& alleles, std::vector<double>& avg_sizes)
{
	for(int label = 0; label < (int)avg_sizes.size(); ++label){
		int n = 0;
		double size = 0;
		for(int i = 0; i < (int)alleles.size(); ++i) {
			if(labels[i] == label){
				++n;
				size += alleles[i].size();
			}
		}
		avg_sizes[label] = (size/n);
	}
}

void find_reps(DistMatrix& matrix, std::vector<std::string>& alleles, std::vector<int>& labels, std::vector<int>& reps)
{
	for(int label = 0; label < (int)reps.size(); ++label){
		std::vector<uint32_t> cluster_indeces;
		for(uint32_t i = 0; i < labels.size(); ++ i) if(labels[i] == label) cluster_indeces.emplace_back(i);
		uint32_t rep_i = matrix.get_min_dist_i(cluster_indeces);
		reps[label] = (int)rep_i;
	}
}


int get_next_node(const std::vector<bool>& visited){
	for(int i = 0; i < (int)visited.size(); ++i) if(!visited[i]) return i;
	return -1;
}

void scaled_genotype(wfa::WFAlignerEdit& aligner, DistMatrix& matrix, std::vector<int>& reps, std::vector<int>& counts, std::vector<std::string>& alleles, std::vector<double>& scaled_labels, bool& length_dist_function)
{
	recompute_pairwise_alignment(aligner, matrix, alleles, reps, length_dist_function);
	int start = 0;
	for(int i = 1; i < (int)reps.size(); ++i) if(alleles[reps[i]].size() > alleles[reps[start]].size()) start = i;

	std::vector<bool> visited(reps.size(), false);
	std::vector<int> path;
	int u = start;
	path.emplace_back(u);
	std::vector<double> dists;
	double total_dist = 0;
	dists.emplace_back(0);
	while(u >= 0){
		//std::cout << '#' << u << '\n';
		visited[u] = true;
		int next_node = -1;
		double next_dist = 1000000.0;
		for(int v = 0; v < (int)reps.size(); ++v){
			if(!visited[v]){
				double alt_dist = matrix.get_dist(reps[u],reps[v]);
				if(alt_dist < next_dist) {
					next_dist = alt_dist;
					next_node = v;
				}
			}
		}
		if(next_node < 0) break;
		else {
			//std::cout << "next node " << next_node << " " << next_dist << '\n';
			path.emplace_back(next_node);
			visited[next_node] = true;
			u = next_node;
			total_dist += next_dist;
			dists.emplace_back(total_dist);
		}
	}
	for(int i = 0; i < (int)path.size(); ++i) scaled_labels[path[i]] = dists[i]/total_dist;
}



void pairwise_process(const OtterOpts& params, const int& ac_mincov, const int& tc_mincov, const bool& is_summary, const bool& is_length, const std::vector<BED>& regions, const std::string& bam, const std::vector<std::string>& index2sample, const std::map<std::string,int>& sample2index)
{
	BS::thread_pool pool(params.threads);
	std::mutex stdout_mtx;
	pool.parallelize_loop(0, regions.size(),
		[&pool, &stdout_mtx, &params, &ac_mincov, &tc_mincov, &is_summary, &is_length, &bam, &regions, &index2sample, &sample2index](const int a, const int b){
			BamInstance bam_inst;
			bam_inst.init(bam, true);
			wfa::WFAlignerEdit aligner(wfa::WFAligner::Score, wfa::WFAligner::MemoryMed);
			//aligner.setHeuristicBandedAdaptive(50, 50, 1);
			for(int region_i = a; region_i < b; ++region_i) {
				const BED& region = regions[region_i];
				const std::string region_str = region.chr + ':' + std::to_string(region.start) + '-' + std::to_string(region.end);
				hts_itr_t *iter = sam_itr_querys(bam_inst.idx, bam_inst.header, region_str.c_str());
				if(iter == nullptr) std::cerr << "(" << antimestamp() << "): WARNING: query failed at region " << region_str << std::endl;
				else{
					std::vector<int> allele_sample_indeces;
					std::vector<std::string> alleles;
					while(sam_itr_next(bam_inst.fp, iter, bam_inst.read) > 0) parse_bam_allele(region_str, ac_mincov, tc_mincov, sample2index, bam_inst.read, alleles, allele_sample_indeces);
					if(alleles.size() > 1){
						sort_bam_alleles(allele_sample_indeces, alleles);
						std::vector<std::pair<int, std::pair<int,int>>> sample2intervals;
						set_sample_intervals(allele_sample_indeces, sample2intervals);
						if(is_length){
							for(const auto& si_pair : sample2intervals){
								int a1_l = alleles[si_pair.second.first].size();
								int a2_l = alleles[si_pair.second.second].size();
								stdout_mtx.lock();
								std::cout << region_str << '\t' << index2sample[si_pair.first] << '\t' << (a1_l < a2_l ? a1_l : a2_l) << '\t' << (a1_l > a2_l ? a1_l : a2_l) << '\t' << (a1_l + a2_l) << '\n';
								stdout_mtx.unlock();
							}
						}
						else{
							if(is_multi_sample(allele_sample_indeces)){
								DistMatrix matrix(alleles.size());
								std::vector<int> labels(alleles.size());
								gt_pairwise_alignment(params.max_error,aligner, matrix, alleles, labels);
								int max_label = 0;
								for(int i = 0; i < (int)alleles.size(); ++i) if(labels[i] > max_label) max_label = labels[i];
								++max_label;
								std::vector<int> label_counts(max_label, 0);
				    			get_cluster_size(labels, label_counts);
								std::vector<double> avg_sizes(max_label);
								get_summary(region_str, stdout_mtx, labels, alleles, avg_sizes);
			    				if(is_summary) {
			    					for(int label = 0; label < max_label; ++label){
			    						stdout_mtx.lock();
										std::cout << region_str << '\t' << label << '\t' << label_counts[label] << '\t' << avg_sizes[label] << '\n';
										stdout_mtx.unlock();
			    					}
			    				}
			    				else{
				    				std::vector<int> reps(max_label);
				    				find_reps(matrix, alleles, labels, reps);
				    				std::vector<double> scaled_labels(max_label);
				    				double u = 0.0, n = 0.0;
									for(const auto& v : avg_sizes) u += v;
									u /= avg_sizes.size();
									for(const auto& v : avg_sizes) n += ((v - u)*(v-u));
									double se = std::sqrt(n/(avg_sizes.size() - 1));
									double large;
									for(const auto& size : avg_sizes) if(size > 1000) ++large;
									large /= index2sample.size();
									bool length_dist_function = false;
									if(large >= 0.5 && se > 500) {
										stdout_mtx.lock();
										std::cerr << "(" << antimestamp() << "): WARNING: " << region_str << " triggered length distance function (F1Kbp=" << large << ",SE=" << se << ')' << std::endl;
										stdout_mtx.unlock();
										length_dist_function = true;
									}
									scaled_genotype(aligner, matrix, reps, label_counts, alleles, scaled_labels, length_dist_function);
					    			for(int i = 0; i < (int)sample2intervals.size(); ++i){
					    				double a1 = scaled_labels[labels[sample2intervals[i].second.first]];
					    				double a2 = scaled_labels[labels[sample2intervals[i].second.second]];
					    				bool is_a1_min = a1 < a2;
					    				stdout_mtx.lock();
					    				std::cout << index2sample[sample2intervals[i].first] << '\t' << region_str << '\t';
					    				if(is_a1_min) std::cout << a1 << '/' << a2; else std::cout << a2 << '/' << a1;
					    				std::cout << '\n';
					    				stdout_mtx.unlock();
					    			}
				    			}
							}
						}
					}
				}
			}
			bam_inst.destroy();
	}).wait();

}

void genotype(const std::string& bam, const std::string& bed, const OtterOpts& params, const int& ac_mincov, const int& tc_mincov, const bool& is_summary, const bool& is_length)
{
	std::vector<BED> regions;
	parse_bed_file(bed, regions);

	std::vector<std::string> sample_index;
	std::map<std::string, int> sample2index;
	index_read_groups(bam, sample2index, sample_index);

	std::cerr << '(' << antimestamp() << "): Found " << sample_index.size() << " samples (read-group tags)\n";

	pairwise_process(params, ac_mincov, tc_mincov, is_summary, is_length, regions, bam, sample_index, sample2index);
}