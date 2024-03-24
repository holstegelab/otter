#include "pairwise_alignment.hpp"
#include "parse_bam_alignments.hpp"
#include "kde.hpp"
#include "dist_matrix.hpp"
#include "ansparc.hpp"
#include "fasta_helper.hpp"
#include <string>
#include <iostream>
#include "bindings/cpp/WFAligner.hpp"
#include "fastcluster.h"


DecisionBound::DecisionBound(double x, double y, double z): dist0(x),dist1(y),cut0(z){};

double otter_seq_dist(const bool& ignore_haps, wfa::WFAligner& aligner, std::string& x, ParsingStatus& x_status, int& x_hp, std::string& y, ParsingStatus& y_status, int& y_hp)
{
	if(!x_status.is_spanning() && !y_status.is_spanning()) return -1;
	else{
		if(!ignore_haps && x_hp >= 0 && y_hp >= 0){
			if(x_hp == y_hp) return 0; 
			else return 1.0;
		} 
		else if(x == y) return 0.0;
		else{
			bool x_smallest = x.size() < y.size();
			double largest = x_smallest ? (double)y.size() : (double)x.size();
			
			if(x_smallest) aligner.alignEnd2End(x, y); else aligner.alignEnd2End(y, x);

			int dist = aligner.getAlignmentScore();
			if(x_status.is_spanning() && y_status.is_spanning()) return dist/largest;
			else {
				int size_diff = x_smallest ? (int)y.size() - (int)x.size() : (int)x.size() - (int)y.size();
				return (dist - size_diff)/largest;
			}
		}
	}	
}

void otter_pairwise_dist(const bool& ignore_haps, const std::vector<int>& indeces, const bool& is_indeces_subset, AlignmentBlock& sequences, wfa::WFAligner& aligner, DistMatrix& distmatrix)
{
	for(int i = 0; i < (int)indeces.size(); ++i){
		for(int j = i + 1; j < (int)indeces.size(); ++j){
			double dist = otter_seq_dist(ignore_haps, aligner, sequences.seqs[indeces[i]], sequences.statuses[indeces[i]], sequences.hps[indeces[i]].hp, sequences.seqs[indeces[j]], sequences.statuses[indeces[j]], sequences.hps[indeces[j]].hp);
			if(is_indeces_subset) distmatrix.set_dist((uint32_t)indeces[i], (uint32_t)indeces[j], dist); 
			else distmatrix.set_dist((uint32_t)i, (uint32_t)j, dist);
		}
	}

}

DecisionBound otter_find_clustering_dist(const double& bandwidth, const std::vector<int>& indeces, AlignmentBlock& sequences, wfa::WFAligner& aligner, DistMatrix& distmatrix)
{
	KDE kde(bandwidth);
	for(const auto& v : distmatrix.values) if(v >= 0.0) kde.values.emplace_back(v);
 	std::vector<double> densities;
	for(double x = 0.0; x <= 1.0; x += bandwidth) densities.emplace_back(kde.f(x));
	//for(const auto& d : densities) std::cout << d << '\n';
	std::vector<std::pair<int,double>> maximas;
	std::vector<std::pair<int,double>> minimas;
	kde.maximas(densities, maximas, minimas);
 	if(maximas.size() == 1) return DecisionBound(maximas[0].first*0.01,maximas[0].first*0.01,-1.0);
	else if(maximas.size() == 2) return DecisionBound(maximas[0].first*0.01, maximas[1].first*0.01, minimas[0].first*0.01); 
	else {
		int m_i = 1;
		for(int i = 2; i < (int)maximas.size(); ++i){
			if(maximas[i].second - maximas[m_i].second > 1.5 && maximas[i].second > maximas[m_i].second) {
				m_i = i;
			}
		}
		return DecisionBound(maximas[0].first*0.01, maximas[m_i].first*0.01, minimas[m_i - 1].first*0.01); 
	}
}

void otter_hclust(const bool& ignore_haps, const int& max_alleles, const double& bandwidth, const double& max_tolerable_diff, const double& min_cov_fraction, const int& min_cov_fraction2_l, const double& min_cov_fraction2_f,const std::vector<int>& spannable_indeces, DistMatrix& distmatrix, wfa::WFAligner& aligner, AlignmentBlock& sequences, std::vector<int>& cluster_labels, int& initial_clusters)
{
	cluster_labels.resize(sequences.names.size(), -1);
	if(spannable_indeces.size() == 1) cluster_labels[spannable_indeces[0]] = 0;
	else if(spannable_indeces.size() == 2){
		double dist = otter_seq_dist(ignore_haps, aligner, sequences.seqs[0], sequences.statuses[0], sequences.hps[0].hp, sequences.seqs[1], sequences.statuses[1], sequences.hps[1].hp);
		cluster_labels[spannable_indeces[0]] = 0;
		cluster_labels[spannable_indeces[1]] = dist > max_tolerable_diff ? 1 : 0;
	}
	else{
		//std::cout << "clustering dist of " << dists.dist0 << ' ' <<  dists.dist1 << ' ' << dists.cut0 << '\n';
		/**
		double* distmat = distmatrix.values.data();
		uint32_t k,i,j;
		for(i=k=0; i < spannable_indeces.size(); ++i) {
			for(j=i+1; j < spannable_indeces.size(); ++j) {
				distmat[k] = distmatrix.get_dist(i, j);
				++k;
			}
		}
		*/
		otter_pairwise_dist(ignore_haps, spannable_indeces, false, sequences, aligner, distmatrix);
		if(max_alleles == 1) {
			for(int i = 0; i < (int)spannable_indeces.size(); ++i) cluster_labels[spannable_indeces[i]] = 0;
			initial_clusters = 1;
		}
		else{
			DecisionBound dists = otter_find_clustering_dist(bandwidth, spannable_indeces, sequences, aligner, distmatrix);
			int* merge = new int[2*(spannable_indeces.size()-1)];
		    double* height = new double[spannable_indeces.size()-1];
		    auto distmatrix_cpy = distmatrix.values;
		    hclust_fast(spannable_indeces.size(), distmatrix_cpy.data(), HCLUST_METHOD_AVERAGE, merge, height);
		    if(dists.dist1 - dists.dist0 <= max_tolerable_diff) {
		    	for(const auto& i : spannable_indeces) cluster_labels[i] = 0;
		    	initial_clusters = 1;
		    }
		    else {
		    	if(dists.cut0 < 0){
		    		std::cerr << "ERROR: unexpected clustering boundary: " << dists.dist0 << ',' << dists.dist1 << ',' << dists.cut0 << '\n';
		    		exit(1);
		    	}
		    	int* labels = new int[spannable_indeces.size()];
		    	double dist_final = dists.dist1 == 0.01 ? dists.dist1 : dists.cut0 + 0.01;
			    cutree_cdist(spannable_indeces.size(), merge, height, dist_final, labels);
			    int total_alleles = 0;
			    for(int i = 0; i < (int)spannable_indeces.size(); ++i) if(labels[i] > total_alleles) total_alleles = labels[i];
			    ++total_alleles;
				initial_clusters = total_alleles;
				int min_cov1 = int(spannable_indeces.size()*min_cov_fraction + 0.5);
				int min_cov2 = int(spannable_indeces.size()*min_cov_fraction2_f + 0.5);
				if(max_alleles != 0) {
					std::vector<int> label_counts(total_alleles);
					std::vector<int> label_max_sizes(total_alleles);

					for(int i = 0; i < (int)spannable_indeces.size(); ++i) {
						++label_counts[labels[i]];
						if((int)sequences.seqs[spannable_indeces[i]].size() > label_max_sizes[labels[i]]) label_max_sizes[labels[i]] = sequences.seqs[spannable_indeces[i]].size();
					}
					int label_max_cov = 0;
					std::vector<int> label_min_covs(total_alleles);
					for(int l = 0; l < (int)label_max_sizes.size(); ++l){
						if(label_counts[l] > label_max_cov) label_max_cov = label_counts[l];
						if(label_max_sizes[l] < min_cov_fraction2_l) label_min_covs[l] = min_cov1;
						else label_min_covs[l] = min_cov2;
					}

					bool is_only_singletons = true;
					for(int l = 0; l < (int)label_min_covs.size(); ++l){
						if(label_counts[l] >= label_min_covs[l]) {
							is_only_singletons = false;
							break;
						}
					}

					if(is_only_singletons) cutree_k(spannable_indeces.size(), merge, max_alleles, labels);
					else{
						int outlier_clusters_n = 0;
						int seed_clusters_n = 0;
						for(int l = 0; l < (int)label_counts.size(); ++l){
							if(label_counts[l] < label_min_covs[l]) ++outlier_clusters_n;
							else ++seed_clusters_n;
						}
						if(seed_clusters_n == 0 || seed_clusters_n > max_alleles) cutree_k(spannable_indeces.size(), merge, max_alleles, labels);
						else{
							std::vector<int> outlier_clusters;
							std::vector<int> seed_clusters;
							for(int l = 0; l < total_alleles; ++l) {
								if(label_counts[l] < label_min_covs[l]) outlier_clusters.emplace_back(l);
								else seed_clusters.emplace_back(l);
							}
							//std::cout << seed_clusters.size() << '\t' << outlier_clusters.size() << '\n';
							/**
							for(const auto& s : seed_clusters) {
								for(int i = 0; i < (int)spannable_indeces.size(); ++i) if(labels[i] == s) std::cout << "seed\t" << label_counts[s] << '\t' << sequences.seqs[spannable_indeces[i]] << '\n';
							}
							for(const auto& s : outlier_clusters) {
								for(int i = 0; i < (int)spannable_indeces.size(); ++i) if(labels[i] == s) std::cout << "outlier\t" << label_counts[s] << '\t' << sequences.seqs[spannable_indeces[i]] << '\n';
							}
						   */

							for(int i = 0; i < (int)spannable_indeces.size(); ++i){
								for(const auto& o : outlier_clusters) {
									if(labels[i] == o){
										labels[i] = -1;
										break;
									}
								}
							}
							
							std::vector<int> readjusted_seed_cluster_labels(seed_clusters.size());
							for(int i = 0; i < (int)seed_clusters.size(); ++i) readjusted_seed_cluster_labels[i] = i;
							for(int i = 0; i < (int)spannable_indeces.size(); ++i){
								for(int j = 0; j < (int)seed_clusters.size(); ++j) {
									if(labels[i] == seed_clusters[j]){
										labels[i] = readjusted_seed_cluster_labels[j];
										break;
									}
								}
							}

							for(int i = 0; i < (int)spannable_indeces.size(); ++i){
								if(labels[i] == -1){
									int closest_j;
									double min_dist = 100000.0;
									for(int j = 0; j < (int)spannable_indeces.size(); ++j){
										if(i != j && labels[j] != -1){
											double j_dist = distmatrix.get_dist(i, j);
											if(j_dist < min_dist){
												closest_j = j;
												min_dist = j_dist;
											}
										} 
									}
									labels[i] = labels[closest_j];
								}
							}
						}
					}
				}
			    for(int i = 0; i < (int)spannable_indeces.size(); ++i) {
			    	cluster_labels[spannable_indeces[i]] = labels[i];
			    }
			    delete[] labels;
		    }

		    delete[] merge;
		    delete[] height;
		}
	}
}

void otter_realignment(const std::string& chr, const int& start, const int& end, const int& flank, const double& min_sim, const FaidxInstance& faidx_inst, AlignmentBlock& sequences, wfa::WFAligner& aligner)
{
	std::string ref_left = "";
	std::string ref_right = "";
	for(int i = 0; i < (int)sequences.names.size(); ++i){
		auto status = sequences.statuses[i];
		if(!status.is_spanning() && (status.spanning_l || status.spanning_r)){
			bool left_realignment = status.spanning_r && status.alignment_coords.first >= flank;
			bool right_realignment = status.spanning_l && (int)sequences.seqs[i].size() - status.alignment_coords.second >= flank;
			std::string subseq = "";
			if(left_realignment) {
				if(ref_left.empty()) faidx_inst.fetch(chr, start - flank, start, ref_left);
				subseq = sequences.seqs[i].substr(0, status.alignment_coords.first);
				aligner.alignEnd2End(ref_left, subseq);
			}
			else if(right_realignment){
				if(ref_right.empty()) faidx_inst.fetch(chr, end, end + flank, ref_right);
				subseq = sequences.seqs[i].substr(status.alignment_coords.second);
				aligner.alignEnd2End(ref_right, subseq);
			}

			if(!subseq.empty()){
				std::vector<int> scores(subseq.size(), 0);
				int j = 0;
				for(const auto& op : aligner.getAlignmentCigar()){
					if(op != 'D'){
						int penalty = op == 'M' ? 1 : -1;
						if(penalty > 0){
							if(j == 0) scores[j] = penalty;
							else scores[j] = scores[j-1] + penalty;
						}
						else if(j > 0 && scores[j-1] > 0) scores[j] = scores[j-1] + penalty;
						++j;
					}
				}
				int max_sum_i = 0;
				for(j = 0; j < (int)scores.size(); ++j) if(scores[j] > scores[max_sum_i]) max_sum_i = j;
				int start_i = max_sum_i;
				while(start_i > 0 && scores[start_i] > 0) --start_i;
				if((scores[max_sum_i] / (double)flank) >= min_sim){
					if(left_realignment) sequences.seqs[i] = sequences.seqs[i].substr(max_sum_i);
					else if(right_realignment) sequences.seqs[i] = sequences.seqs[i].substr(0, status.alignment_coords.second + start_i);
					sequences.statuses[i].spanning_l = true;
					sequences.statuses[i].spanning_r = true;
				}
			}
		}
	}
}

void _set_unique_labels(std::vector<int>& unique_labels)
{
	sort(unique_labels.begin(), unique_labels.end());
	unique_labels.erase(std::unique( unique_labels.begin(), unique_labels.end()), unique_labels.end());
	auto it = unique_labels.begin();
	while(it != unique_labels.end()){
		if(*it == -1) it = unique_labels.erase(it);
		else ++it;
	}
}

void otter_nonspanning_assigment(const double& ignore_haps, const double& min_sim, const double& max_error, AlignmentBlock& sequences, wfa::WFAligner& aligner, std::vector<int>& labels)
{
	std::vector<int> unique_labels = labels;
	_set_unique_labels(unique_labels);
	for(int i = 0; i < (int)labels.size(); ++i){
		if(labels[i] < 0){
			std::vector<double> max_sim(unique_labels.size(), 0.0);
			int i_l = sequences.seqs[i].size();
			const auto& i_s = sequences.statuses[i];
			for(int j = 0; j < (int)labels.size(); ++j){
				if(i != j && labels[j] >= 0) {
					int j_l = sequences.seqs[j].size();
					if(i_l >= j_l) aligner.alignEnd2End(sequences.seqs[j], sequences.seqs[i]);
					else {
						if(i_s.spanning_l) aligner.alignEndsFree(sequences.seqs[i], 0, 0, sequences.seqs[j], 0, j_l - i_l);
						else if(i_s.spanning_r) aligner.alignEndsFree(sequences.seqs[i], 0, 0, sequences.seqs[j], j_l - i_l, 0);
						else aligner.alignEndsFree(sequences.seqs[i], 0, 0, sequences.seqs[j], (j_l - i_l)/2, (j_l - i_l)/2);
					}
					double sim = 1 - (aligner.getAlignmentScore() / (double)i_l);
					if(!ignore_haps && sequences.hps[i].hp >= 0 && sequences.hps[j].hp >= 0) {
						if(sequences.hps[i].hp == sequences.hps[j].hp) sim = 1.1;
					}
					if(sim > max_sim[labels[j]]) max_sim[labels[j]] = sim;
				}
			}
			//for(const auto& s : max_sim) std::cout << '\t' << s; std::cout << '\n';
			int max_sim_label = 0;
			for(int j = 1; j < (int)max_sim.size(); ++j) if(max_sim[j] > max_sim[max_sim_label]) max_sim_label = j;
			if(max_sim[max_sim_label] >= min_sim){
				double min_diff = 1.0;
				for(int j = 0; j < (int)max_sim.size(); ++j) {
					if(max_sim_label != j) {
						double diff = max_sim[max_sim_label] - max_sim[j];
						if(diff < min_diff) min_diff =  diff;
					}
				}
				if(min_diff >= max_error) labels[i] = max_sim_label;
			}
		}
	}
}

void otter_rapid_consensus(const std::vector<int>& spannable_indeces, const std::vector<int>& labels, const DistMatrix& distmatrix, wfa::WFAligner& aligner, AlignmentBlock& sequences, std::vector<std::string>& consensus_seqs, std::vector<std::vector<double>>& ses)
{
	std::vector<int> unique_labels = labels;
	_set_unique_labels(unique_labels);
	consensus_seqs.resize(unique_labels.size());
	ses.resize(unique_labels.size());
	for(int label_i = 0; label_i < (int)unique_labels.size(); ++label_i){
		int label = unique_labels[label_i];
		if(label >= 0){
			std::vector<int> spannable_indeces_cluster;
			std::vector<int> spannable_indeces_cluster_l;
			for(int i = 0; i < (int)spannable_indeces.size(); ++i) {
				if(labels[spannable_indeces[i]] == label) {
					spannable_indeces_cluster.emplace_back(i);
					spannable_indeces_cluster_l.emplace_back(sequences.seqs[i].size());
				}
			}
			std::vector<double> local_distances;
			if(spannable_indeces_cluster.size() < 3) {
				consensus_seqs[label] = sequences.seqs[spannable_indeces[spannable_indeces_cluster.front()]];
			}
			else {
				for(int i = 0; i < (int)spannable_indeces_cluster.size(); ++i){
					for(int j = i + 1; j < (int)spannable_indeces_cluster.size(); ++j){
						local_distances.emplace_back(distmatrix.get_dist(spannable_indeces_cluster[i], spannable_indeces_cluster[j]));
					}
				}
				
				ses[label] = local_distances;

				bool is_k_short = false;
				for(const auto& l : spannable_indeces_cluster_l) {
					if(l < 5){
						is_k_short = true;
						break;
					}
				}
				std::vector<uint32_t> local_indeces;
				for(const auto& i : spannable_indeces_cluster) local_indeces.emplace_back((uint32_t)i);
				uint32_t representative_i = distmatrix.get_min_dist_i(local_indeces);
				uint32_t k = 3;
				uint32_t g = 2;
				uint32_t l = 2;
				if(is_k_short){
					k = 1;
					g = 1;
					l = 1;
				}
				ansparc_graph graph;
				graph.init(sequences.seqs[spannable_indeces[representative_i]], k, g, l);
				for(const auto& i : local_indeces){
					if(i != representative_i){
						aligner.alignEnd2End(sequences.seqs[spannable_indeces[representative_i]], sequences.seqs[spannable_indeces[i]]);
						std::string cigar = aligner.getAlignmentCigar();
						graph.insert_alignment(cigar, sequences.seqs[spannable_indeces[i]]);
					}
				}
				//adjust weights
				float c = 2.0;
				if(spannable_indeces_cluster.size() < 4) c = 1.0f;
				float t = 0.2;
				graph.adjust_weights(c, t);
				graph.consensus(consensus_seqs[label]);
			}
		}
	}
}