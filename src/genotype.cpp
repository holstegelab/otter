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


void append_ref_sample(const BED& region, const int& offset, const int& ref_allele_i, std::mutex& stdout_mtx, FaidxInstance& faidx_inst, std::vector<std::string>& alleles, std::vector<int>& sample_indeces)
{
	std::string ref_allele;
	faidx_inst.fetch(region.chr, region.start - offset, region.end + offset - 1, ref_allele);
	if(ref_allele.empty()){
		stdout_mtx.lock();
		std::cerr << "(" << antimestamp() << "): WARNING: could not find reference sequence allele for " << region.toBEDstring() << std::endl;
		stdout_mtx.unlock();
	}
	else{
		alleles.emplace_back(ref_allele);
		sample_indeces.emplace_back(ref_allele_i);
	}
}

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
	auto matrix_values_cpy = matrix.values;
	hclust_fast(alleles.size(), matrix_values_cpy.data(), HCLUST_METHOD_AVERAGE, merge, height);
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

void output_vcf_header(const std::string& bam, const std::vector<std::string>& sample_index, const std::string& ref_name)
{
	BamInstance bam_inst;
	bam_inst.init(bam, true);
	std::cout << "##fileformat=VCFv4.2\n";
	for(int i = 0; i < bam_inst.header->n_targets; ++i) std::cout << "##contig=<ID=" << bam_inst.header->target_name[i] << ",length=" << bam_inst.header->target_len[i] << ">\n";
	bam_inst.destroy();
	std::cout << 
	"##INFO=<ID=GTS,Number=1,Type=String,Description=\"Multi-Allelic Scaled Genotype Array\">\n" <<
	"##ALT=<ID=DEL,Description=\"Deletion\">\n" <<
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	for(const auto& sample : sample_index) if(sample != ref_name) std::cout << '\t' << sample;
	std::cout << '\n';
}

void output_vcf(const BED& region, const int& offset, const std::vector<std::string>& index2sample, const std::map<std::string,int>& sample2index, const std::vector<std::pair<int, std::pair<int,int>>>& sample2intervals, const std::vector<int>& labels, const std::vector<double>& scaled_labels, const std::vector<std::string>& alleles, const std::vector<int>& reps, const int& ref_allele_i)
{
	std::vector<int> index2sampleintervals(index2sample.size(), -1);
	for(int i = 0; i < (int)sample2intervals.size(); ++i) index2sampleintervals[sample2intervals[i].first] = i;
	std::cout << region.chr << '\t' << (1 + region.start - offset) << '\t' << region.toBEDstring() << '\t';
	if(ref_allele_i < 0) std::cout << "N\t"; else std::cout << alleles[ref_allele_i] << '\t';
	if(reps.size() == 1) {if(ref_allele_i < 0) std::cout << alleles[reps[0]]; else std::cout << '.';}
	else{
		int ref_label = -1;
		if(ref_allele_i >= 0) ref_label = labels[ref_allele_i];
		bool past_first = false;
		for(int i = 0; i < (int)reps.size(); ++i){
			if(i != ref_label) {
				if(past_first) std::cout << ',';
				else past_first = true;
				const std::string& local_allele = alleles[reps[i]];
				std::cout << (local_allele == "N" ? "<DEL>" : local_allele);
			}
		}
	}
	std::cout << "\t.\t.\tGTS=";
	for(int i = 0; i < (int)reps.size(); ++i){
		if(i > 0) std::cout << ',';
		std::cout << scaled_labels[i];
	}
	std::cout << "\tGT\t";
	int samples_processed = 0;
	for(int i = 0; i < (int)index2sample.size(); ++i){
		if(i != ref_allele_i){
			const auto& interval_i = index2sampleintervals[i];
			if(samples_processed > 0) std::cout << '\t';
			if(interval_i < 0) std::cout << "./.";
			else {
				const auto& interval = sample2intervals[interval_i];
				const int& a1 = labels[interval.second.first];
				const int& a2 = labels[interval.second.second];
				std::cout << (a1 < a2 ? a1 : a2) << '/' << (a1 > a2 ? a1 : a2);
			}
			++samples_processed;
		}
	}
	std::cout << '\n';
}

void center_ref_labels(const int& max_label, const int& ref_label, std::vector<int>& labels)
{
	std::vector<int> new_labels(max_label);
	for(int i = 0; i < max_label; ++i) new_labels[i] = i;
	std::sort(new_labels.begin(), new_labels.end(),[&ref_label](const int& x, const int& y){ return (x == ref_label || (y == ref_label ? false : x < y)); });
	std::vector<int> old2new(max_label);
	for(int i = 0; i < max_label; ++i) old2new[new_labels[i]] = i;
	for(int i = 0; i < (int)labels.size(); ++i) labels[i] = old2new[labels[i]];
}


void pairwise_process(const OtterOpts& params, const int& ac_mincov, const int& tc_mincov, const bool& is_summary, const bool& is_length, const std::vector<BED>& regions, const std::string& bam, const std::vector<std::string>& index2sample, const std::map<std::string,int>& sample2index, const std::string& reference, const int& ref_allele_i)
{
	BS::thread_pool pool(params.threads);
	std::mutex stdout_mtx;
	pool.parallelize_loop(0, regions.size(),
		[&pool, &stdout_mtx, &params, &ac_mincov, &tc_mincov, &is_summary, &is_length, &bam, &regions, &index2sample, &sample2index, &reference, &ref_allele_i](const int a, const int b){
			BamInstance bam_inst;
			bam_inst.init(bam, true);
			wfa::WFAlignerEdit aligner(wfa::WFAligner::Score, wfa::WFAligner::MemoryMed);
			FaidxInstance faidx_inst;
			bool ref_is_fa = ref_allele_i >= 0 && !reference.empty();
			if(ref_is_fa) faidx_inst.init(reference);
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
						if(ref_is_fa) append_ref_sample(region, params.offset, ref_allele_i, stdout_mtx, faidx_inst, alleles, allele_sample_indeces);
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
								int local_ref_allele_i = -1;
								for(const auto& refinterval : sample2intervals) if(refinterval.first == ref_allele_i) local_ref_allele_i = refinterval.second.first;	
								if(local_ref_allele_i >= 0) center_ref_labels(max_label, ref_allele_i, labels);
								else if(local_ref_allele_i >= 0) std::cerr << '(' << antimestamp() << "): WARNING: reference sequence for " << region_str << " not found\n";
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
				    				if(max_label <= 2) for(int i = 0; i < max_label; ++i) scaled_labels[i] = i;
				    				else{
				    					double u = 0.0, n = 0.0;
										for(const auto& v : avg_sizes) u += v;
										u /= avg_sizes.size();
										for(const auto& v : avg_sizes) n += ((v - u)*(v-u));
										double se = std::sqrt(n/(avg_sizes.size() - 1));
										double large = 0.0;
										for(const auto& size : avg_sizes) if(size > params.dist_length_large) ++large;
										large /= index2sample.size();
										bool length_dist_function = false;
										if(large >= params.dist_length_frac && params.dist_length_se > 500) {
											stdout_mtx.lock();
											std::cerr << "(" << antimestamp() << "): WARNING: " << region_str << " triggered length distance function (F1Kbp=" << large << ",SE=" << se << ')' << std::endl;
											stdout_mtx.unlock();
											length_dist_function = true;
										}
										scaled_genotype(aligner, matrix, reps, label_counts, alleles, scaled_labels, length_dist_function);
				    				}
				    				if(ref_allele_i >= 0){
				    					stdout_mtx.lock();
				    					output_vcf(region, params.offset, index2sample, sample2index, sample2intervals, labels, scaled_labels, alleles, reps, local_ref_allele_i);
				    					stdout_mtx.unlock();
				    				}
				    				else{
				    					for(int i = 0; i < (int)sample2intervals.size(); ++i){
						    				double a1 = labels[sample2intervals[i].second.first];
						    				double a2 =labels[sample2intervals[i].second.second];
						    				double sa1 = scaled_labels[a1];
						    				double sa2 = scaled_labels[a2];
						    				//bool is_a1_min = a1 < a2;
						    				stdout_mtx.lock();
						    				std::cout << index2sample[sample2intervals[i].first] << '\t' << region_str << '\t' << a1 << '/' << a2 << '\t' << sa1 << '/' << sa2 << '\n';
						    				//if(is_a1_min) std::cout << a1 << '/' << a2; else std::cout << a2 << '/' << a1;
						    				//std::cout << '\n';
						    				stdout_mtx.unlock();
						    			}
				    				}
				    			}
							}
						}
					}
				}
			}
			faidx_inst.destroy();
			bam_inst.destroy();
	}).wait();

}

void genotype(const std::string& bam, const std::string& bed, OtterOpts& params, const int& ac_mincov, const int& tc_mincov, const bool& is_summary, const bool& is_length, const std::string& reference, const std::string& vcf)
{
	std::vector<BED> regions;
	parse_bed_file(bed, regions);

	std::vector<std::string> sample_index;
	std::map<std::string, int> sample2index;
	index_read_groups(bam, sample2index, sample_index);
	std::cerr << '(' << antimestamp() << "): Found " << sample_index.size() << " samples (read-group tags)\n";

	int offset = -1;
	fetch_preset_offset(bam, offset);
	if(offset < 0) params.offset = 0;
	else{
		std::cerr << '(' << antimestamp() << "): Found preset offset size of " << offset <<'\n';
		params.offset = offset;
	}

	std::string ref_name;
	if(!reference.empty()) ref_name = antimestamp();
	else if(!vcf.empty()) ref_name = vcf;

	int reference_i = -1;
	if(!reference.empty()){
		auto it = sample2index.find(ref_name);
		if(it != sample2index.end()) {
			std::cerr << '(' << antimestamp() << "): ERROR: duplicate reference name: " << ref_name << '\n';
			exit(1);
		}
		else {
			reference_i = (int)sample_index.size();
			sample2index[ref_name] = reference_i;
			sample_index.emplace_back(ref_name);
		};
	}
	else if(!vcf.empty()){
		auto it = sample2index.find(ref_name);
		if(it == sample2index.end()) std::cerr << '(' << antimestamp() << "): WARNING: could not find reference genome (" << ref_name << ") in input BAM file\n";
		else reference_i = it->second;
	}

	if(reference_i >= 0) {
		std::cerr << '(' << antimestamp() << "): Genotyping with reference genome " << (reference.empty() ? ref_name : reference) << '\n';
		output_vcf_header(bam, sample_index, ref_name);
	}

	pairwise_process(params, ac_mincov, tc_mincov, is_summary, is_length, regions, bam, sample_index, sample2index, reference, reference_i);
}