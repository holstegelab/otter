#include "otterclust.hpp"
#include "bindings/cpp/WFAligner.hpp"
#include "fastcluster.h"
#include "ankde.hpp"
#include "anseqs.hpp"
#include "andistmat.hpp"
#include <vector>
#include <string>
#include <list>
#include <cmath>

ClusteringStatus::ClusteringStatus(): ic(0), fc(0){};

Genotype::Genotype(): gt(-1),gt_l(-1),gt_k(-1), hsd(-1){};

void ClusteringStatus::set_global_label(int l){for(uint32_t i = 0; i < labels.size(); ++i) labels[i] = l;};

DecisionBound::DecisionBound(double x, double y, double z): dist0(x),dist1(y),cut0(z){};

DecisionBound otter_find_clustering_dist(const int& radius, const double& dinterval, const double& bandwidth, DistMatrix& distmatrix)
{
	KDE kde(bandwidth);
	for(const auto& v : distmatrix.values) kde.values.emplace_back(v);

 	std::vector<double> densities;
	for(double x = 0.0; x <= 1.0; x += dinterval) {
		densities.emplace_back(kde.f(x));
	}
	double total = 0.0;
	for(const auto& d : densities) total += d;
	for(int i = 0; i < (int)densities.size(); ++i) {
		densities[i] = densities[i]/total;
		//std::cout << (i*dinterval) << '\t' << densities[i] << '\n';
	}

	std::vector<std::pair<int,double>> maximas;
	std::vector<std::pair<int,double>> minimas;
	kde.maximas(radius, densities, maximas, minimas);
	if(maximas.empty()){
 		std::cerr << "ERROR: failed to obtain maximas" << std::endl;
		exit(1);
 	}

	/**
	for(const auto& m : maximas) std::cerr << '[' << (m.first*dinterval) << ',' << m.second << "] ";
	std::cerr << std::endl;
	for(const auto& m : minimas) std::cerr << '[' << (m.first*dinterval) << ',' << m.second << "] ";
	std::cerr << std::endl;
	*/
	
 	if(maximas.size() == 1) return DecisionBound(maximas[0].first*dinterval,maximas[0].first*dinterval,-1.0);
 	else {
 		if(minimas.empty()){
 			std::cerr << "ERROR: failed to obtain minimas" << std::endl;
		    exit(1);
 		}
 		if(maximas.size() == 2) return DecisionBound(maximas[0].first*dinterval,maximas[1].first*dinterval,minimas[0].first*dinterval);
		else {
			std::vector<int> sorted_maximas(maximas.size());
			for(int i = 0; i < (int)sorted_maximas.size(); ++i) sorted_maximas[i] = i;
			sort(sorted_maximas.begin(), sorted_maximas.end(), [&maximas](const auto& a, const auto& b){ 
				double diff = maximas[a].second - maximas[b].second;
				diff = diff > 0 ? diff : -diff;
				if(diff <= 0.01) return maximas[a].first < maximas[b].first;
				else return maximas[a].second > maximas[b].second;
			});

			/**
			for(const auto& s : sorted_maximas) std::cerr << '[' << (maximas[s].first*dinterval) << ',' << maximas[s].second << "] ";
			std::cerr << "\n";
			*/
			
			int last_i = 0;
			int acc_i = 1;
			while(acc_i < (int)sorted_maximas.size()){
				//std::cerr << "in\t" << maximas[sorted_maximas[acc_i]].first*dinterval << '\n';
				int index_diff = acc_i > last_i ? acc_i - last_i : last_i - acc_i;
				double f_diff = maximas[sorted_maximas[acc_i]].second - maximas[sorted_maximas[last_i]].second;
				f_diff = f_diff < 0 ? -f_diff : f_diff;
				//std::cerr << index_diff << '\t' << f_diff << '\n';
				if(index_diff == 1 && f_diff <= 0.01){
					//std::cerr << "deleting\t" << maximas[sorted_maximas[acc_i]].first*dinterval << '\n';
					sorted_maximas.erase(sorted_maximas.begin()+acc_i);
					last_i = acc_i;
				}
				++acc_i;
			}

			if(sorted_maximas.size() < 2) return DecisionBound(maximas[0].first*dinterval,maximas[1].first*dinterval,minimas[0].first*dinterval);
			else{
				int m_first_i = sorted_maximas[0];
				int m_second_i = sorted_maximas[1];
				if(m_first_i > m_second_i){
					int tmp = m_first_i;
					m_first_i = m_second_i;
					m_second_i = tmp;
				}
				int boundary_i = m_second_i - 1;
				if(boundary_i < 0 || boundary_i >= (int)minimas.size()){
					std::cerr << "ERROR: unexpected index for minimas: " << boundary_i << std::endl;
		    		exit(1);
				}
				if(m_second_i - m_first_i > 1 && m_second_i - 2 >= 0 && (maximas[m_second_i].first*dinterval - minimas[boundary_i].first*dinterval <= 0.01)) {
					boundary_i = m_second_i - 2;
					if(boundary_i < 0 || boundary_i >= (int)minimas.size()){
						std::cerr << "ERROR: unexpected index for minimas after correction: " << boundary_i << std::endl;
		    			exit(1);
					}
				}
				//std::cout << m_first_i << '\t' << m_second_i << '\t' << (m_second_i - m_first_i)<< '\n';
				//std::cout << maximas[m_first_i].first*dinterval << '\t' << maximas[m_second_i].first*dinterval << '\t' << minimas[m_first_i + (m_second_i - m_first_i)/2].first*dinterval << '\n';
				return DecisionBound(maximas[m_first_i].first*dinterval, maximas[m_second_i].first*dinterval, minimas[m_first_i + (m_second_i - m_first_i)/2].first*dinterval); 
			}
		}
 	}
}

void otter_hclust(const bool& ignore_haps, const int& max_alleles, const double& bandwidth_short, const int& bandwidth_length, const double& bandwidth_long, const double& max_tolerable_diff, const double& min_cov_fraction, const int& min_cov_fraction2_l, const double& min_cov_fraction2_f,const std::vector<int>& indeces, DistMatrix& distmatrix, wfa::WFAligner& aligner, std::vector<ANREAD>& reads, ClusteringStatus& clustering)
{
	//note: labels is indexed against indeces
	clustering.labels.resize(indeces.size(), -1);
	/** single read situation **/
	if(indeces.size() == 1) {
		clustering.labels[0] = 0;
		clustering.ic = 1;
		clustering.fc = 1;
	}
	/** only two read situation **/
	else if(indeces.size() == 2){
		clustering.labels[0] = 0;
		clustering.labels[1] = 0;
		if(max_alleles == 1) {
			clustering.ic = 1;
			clustering.fc = 1;
		}
		else{
			double dist = distmatrix.get_dist(0, 1);
			if(dist <= max_tolerable_diff) {
				clustering.ic = 1;
				clustering.fc = 1;
			}
			else {
				clustering.labels[1] = 1;
				clustering.ic = 2;
				clustering.fc = 2;
			}
		}
	}
	/** intitiate clustering (three or more reads) **/
	else{
		//no need to cluster
		if(max_alleles == 1) {
			clustering.set_global_label(0);
			clustering.ic = 1;
			clustering.fc = 1;
		}
		else{
			const double error_intervals = 0.0025;
			//radius to calculate KDE surrounding an error rate
			int radius = int(max_tolerable_diff/error_intervals);
			radius = radius < 1 ? 1 : radius;
			double bandwidth = bandwidth_short;
			for(const auto& i : indeces) {
				if((int)reads[i].seq.size() >= bandwidth_length){
					bandwidth = bandwidth_long;
					break;
				}
			}
			DecisionBound dists = otter_find_clustering_dist(radius, error_intervals, bandwidth, distmatrix);
			//std::cout << dists.dist0 << ',' << dists.dist1 << ',' << dists.cut0 << '\n';
			//too small to separate based on user input
			if(dists.dist1 - dists.dist0 <= max_tolerable_diff){
				clustering.set_global_label(0);
				clustering.ic = 1;
				clustering.fc = 1;
			}
			else{
				int* labels = new int[indeces.size()];
				int* merge = new int[2*(indeces.size()-1)];
			    double* height = new double[indeces.size()-1];
			    auto distmatrix_cpy = distmatrix.values;
			    hclust_fast(indeces.size(), distmatrix_cpy.data(), HCLUST_METHOD_AVERAGE, merge, height);
			    //set distatnce to cut tree
		    	double dist_final = dists.dist1 == bandwidth ? dists.dist1 : dists.cut0 + 0.0025;
				cutree_cdist(indeces.size(), merge, height, dist_final, labels);
				/** set initial clusters (alleles) **/
				int total_alleles = 0;
			    for(int i = 0; i < (int)indeces.size(); ++i) if(labels[i] > total_alleles) total_alleles = labels[i];
			    ++total_alleles;
				clustering.ic = total_alleles;
				int min_cov1 = int(indeces.size()*min_cov_fraction + 0.5);
				int min_cov2 = int(indeces.size()*min_cov_fraction2_f + 0.5);

				if(max_alleles != 0) {
					//total reads per cluter label
					std::vector<int> label_counts(total_alleles);
					//max seq size observed for each cluster label
					std::vector<int> label_max_sizes(total_alleles);
					//assigned coverage threshold per cluster label 
					std::vector<int> label_required_covs(total_alleles);
					for(int i = 0; i < (int)indeces.size(); ++i) {
						++label_counts[labels[i]];
						if((int)reads[indeces[i]].seq.size() > label_max_sizes[labels[i]]) label_max_sizes[labels[i]] = reads[indeces[i]].seq.size();
					}

					/** set the required minimum coverage threshold per cluster based on the max seq size observed */
					int label_max_cov = 0;
					for(int l = 0; l < (int)label_max_sizes.size(); ++l){
						if(label_counts[l] > label_max_cov) label_max_cov = label_counts[l];
						if(label_max_sizes[l] < min_cov_fraction2_l) label_required_covs[l] = min_cov1;
						else label_required_covs[l] = min_cov2;
					}

					//for(int i = 0; i < total_alleles; ++i) std::cout << i << '\t' << label_counts[i] << '\t' << label_max_sizes[i] << '\t' << label_required_covs[i] << '\n';

					/** check if there are no clusters meeting coverage requirement */
					bool is_only_singletons = true;
					for(int l = 0; l < (int)label_required_covs.size(); ++l){
						if(label_counts[l] >= label_required_covs[l]) {
							is_only_singletons = false;
							break;
						}
					}
					//if above is true, cluster and report as-is
					if(is_only_singletons) {
						int* labels2 = new int[indeces.size()];
						cutree_k(indeces.size(), merge, max_alleles, labels2);
						clustering.fc = max_alleles;
						for(uint32_t i = 0; i < indeces.size(); ++i) labels[i] = labels2[i];
						delete[] labels2;
					}
					else{
						/** count total seed_clusters (valid clusters) and outlier clusters (invalid cluster that need to be adjusted) */
						int outlier_clusters_n = 0;
						int seed_clusters_n = 0;
						for(int l = 0; l < (int)label_counts.size(); ++l){
							if(label_counts[l] < label_required_covs[l]) ++outlier_clusters_n;
							else ++seed_clusters_n;
						}
						//cluster and report as-is (no valid clusters or total valid clusters > user input)
						if(seed_clusters_n == 0 || seed_clusters_n > max_alleles) {
							cutree_k(indeces.size(), merge, max_alleles, labels);
							clustering.fc = max_alleles;
						}
						else{
							/** identify cluster labels that are valid (seed_clusters) or invalid (outlier_clusters) **/
							std::vector<int> outlier_clusters;
							std::vector<int> seed_clusters;
							for(int l = 0; l < total_alleles; ++l) {
								if(label_counts[l] < label_required_covs[l]) outlier_clusters.emplace_back(l);
								else seed_clusters.emplace_back(l);
							}

							for(int i = 0 ; i < (int)indeces.size(); ++i){
								for(const auto& outlier_label : outlier_clusters){
									if(labels[i] == outlier_label){
										labels[i] = -1;
										break;
									}
								}
							}
							
							/**
							std::cout << seed_clusters.size() << '\t' << outlier_clusters.size() << '\n';
							
							for(const auto& s : seed_clusters) {
								for(int i = 0; i < (int)indeces.size(); ++i) if(labels[i] == s) std::cout << "seed\t" << label_counts[s] << '\t' << reads[indeces[i]].seq << '\n';
							}
							for(const auto& s : outlier_clusters) {
								for(int i = 0; i < (int)indeces.size(); ++i) if(labels[i] == s) std::cout << "outlier\t" << label_counts[s] << '\t' << reads[indeces[i]].seq << '\n';
							}
							*/

							/** reassign membership of reads from outlier clusters to closest seed cluster **/

							/** first, readjust seed cluster labels to assure linearly increasing label ints **/
							std::vector<int> readjusted_seed_cluster_labels(seed_clusters.size());
							for(int i = 0; i < (int)seed_clusters.size(); ++i) readjusted_seed_cluster_labels[i] = i;
							for(int i = 0; i < (int)indeces.size(); ++i){
								for(int j = 0; j < (int)seed_clusters.size(); ++j) {
									if(labels[i] == seed_clusters[j]){
										labels[i] = readjusted_seed_cluster_labels[j];
										break;
									}
								}
							}

							/** reassign reads to closest cluster **/
							for(int i = 0; i < (int)indeces.size(); ++i){
								if(labels[i] == -1){
									int closest_j;
									double min_dist = 100000.0;
									for(int j = 0; j < (int)indeces.size(); ++j){
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

							clustering.fc = seed_clusters_n;
							
						}
					}
				}

			    for(int i = 0; i < (int)indeces.size(); ++i) clustering.labels[i] = labels[i];

			    delete[] labels;
			    delete[] merge;
		    	delete[] height;
			}
		}
	}
}

double length_dist(const uint32_t& x, const uint32_t& y)
{
	bool is_x_smallest = x < y;
	double dist = is_x_smallest ? y - x : x - y;
	return is_x_smallest ? dist / y : dist / x;
}

void cluter_to_e(const double& max_error, const uint32_t& total_alleles, DistMatrix& distmatrix, std::vector<std::vector<int>>& clusters)
{
	/** cluster **/
	int* labels = new int[total_alleles];
	int* merge = new int[2*(total_alleles-1)];
    double* height = new double[total_alleles-1];
    auto distmatrix_cpy = distmatrix.values;
    hclust_fast(total_alleles, distmatrix_cpy.data(), HCLUST_METHOD_AVERAGE, merge, height);
    cutree_cdist(total_alleles, merge, height, max_error, labels);
    /** seperate clusters into vectors **/
    uint32_t total_clusters = 0;
    for(uint32_t i = 0; i < total_alleles; ++i) if(labels[i] > (int)total_clusters) total_clusters = labels[i];
    ++total_clusters;
    clusters.resize(total_clusters);
    for(uint32_t l = 0; l < total_clusters; ++l){
    	for(uint32_t i = 0; i < total_alleles; ++i) if(labels[i] == (int)l) clusters[l].emplace_back(i);
    }
	delete[] labels;
    delete[] merge;
	delete[] height;
}

void remap_cluster_indeces(const DistMatrix& distmatrix, const std::vector<int>& indeces, const std::vector<std::vector<int>>& input_clusters, std::vector<std::vector<int>>& output_clusters, std::vector<int>& medoids)
{
	//reindexed to alleles
	for(const auto& cluster : input_clusters){
		output_clusters.emplace_back();
		auto& mapped_cluster = output_clusters.back();
		for(const auto& i : cluster) mapped_cluster.emplace_back(indeces[i]);
		if(mapped_cluster.size() <= 2) medoids.emplace_back(mapped_cluster[0]); 
		else {
			std::vector<uint32_t> tmpvec(mapped_cluster.size(), -1);
			for(uint32_t i = 0; i < mapped_cluster.size(); ++i) tmpvec[i] = mapped_cluster[i];
			medoids.emplace_back(distmatrix.get_medoid(tmpvec));
		}
	}
}

void anallele_cluster_length(const double& max_error, const std::vector<ANALLELE>& alleles, const std::vector<int>& indeces, DistMatrix& distmatrix, std::vector<std::vector<int>>& clusters, std::vector<int>& reps)
{
	//indexed to indeces
	for(uint32_t i = 0; i < indeces.size(); ++i){
		const auto& i_l = alleles[indeces[i]].seq.size();
		for(uint32_t j = i+1; j < indeces.size(); ++j){
			const auto& j_l = alleles[indeces[j]].seq.size();
			distmatrix.set_dist(i, j, length_dist(i_l, j_l));
		}
	}
	//cluster using indeces
	std::vector<std::vector<int>> _clusters;
	cluter_to_e(max_error, indeces.size(), distmatrix, _clusters);
	//reindexed to alleles and get representative
	remap_cluster_indeces(distmatrix, indeces, _clusters, clusters, reps);
}

void generate_kusage(const uint32_t& k, const KmerEncoding& encoding, const std::vector<ANALLELE>& alleles, const std::vector<int>& indeces, std::vector<KUSAGE>& kusages)
{
	for(uint32_t i = 0; i < indeces.size(); ++i){
		const auto& allele = alleles[indeces[i]];
		std::vector<double> tmpvec;
		seq2kcounts(k, encoding, allele.seq, tmpvec);
		/**
		if(allele.seq.size() >= k){
			for(uint32_t j = 0; j < allele.seq.size() - k + 1; ++j){
				uint64_t k_j = encoding.kmer2index(allele.seq.c_str() + j, k);
				++tmpvec[k_j];
			}
		}
		*/
		kusages.emplace_back(tmpvec);
	}
}

void anallele_cluster_kusage(const double& max_error, const uint32_t& k, const std::vector<ANALLELE>& alleles, const std::vector<int>& indeces, DistMatrix& distmatrix, std::vector<KUSAGE>& kusages, std::vector<std::vector<int>>& clusters, std::vector<int>& reps)
{
	KmerEncoding encoding;
	generate_kusage(k, encoding, alleles, indeces, kusages);

	for(int i = 0; i < (int)kusages.size(); ++i){
		const auto& i_k = kusages[i];
		for(int j = i+1; j < (int)kusages.size(); ++j){
			const auto& j_k = kusages[j];
			double dist = 1.0 - ((std::isnan(i_k.vnorm) || std::isnan(j_k.vnorm)) ? 0 : (std::round(i_k.cosine_sim(j_k)*1000.0)/1000.0));
			distmatrix.set_dist(i, j, dist);
		}
	}
	//cluster using indeces
	std::vector<std::vector<int>> _clusters;
	cluter_to_e(max_error, kusages.size(), distmatrix, _clusters);
	//reindexed to alleles and get representative
	remap_cluster_indeces(distmatrix, indeces, _clusters, clusters, reps);
}

/**
int _anallele_cluster(const double& max_error_l, const double& max_error_c, const std::vector<ANALLELE>& alleles, std::vector<Genotype>& genotypes, std::vector<int>& gt_reps)
{
	std::vector<int> allele_indeces(alleles.size());
	for(uint32_t i = 0; i < allele_indeces.size(); ++i) allele_indeces[i] = i;
	std::vector<std::vector<int>> length_clusters;
	anallele_cluster_length(max_error_l, alleles, allele_indeces, length_clusters);
	int acc_gt = 0;
	for(uint32_t i = 0; i < length_clusters.size(); ++i){
		const auto& length_cluster = length_clusters[i];
		for(const auto& m : length_cluster) genotypes[m].gt_l = i;
		if(length_cluster.size() == 1) {
			auto& gt = genotypes[length_cluster[0]];
			gt.gt_k = 0;
			gt.gt = acc_gt;
			gt_reps.emplace_back(length_cluster[0]);
			++acc_gt;
		}
		else {
			std::vector<std::vector<int>> kusage_clusters;
			std::vector<int> kusage_reps;
			anallele_cluster_kusage(max_error_c, 3, alleles, length_cluster, kusage_clusters, kusage_reps);
			if(kusage_reps.size() != kusage_clusters.size()){
				std::cerr << "[ERROR] unexpected representative alleles (" << gt_reps.size() << ") for " << kusage_clusters.size() << " kusage clusters" << std::endl;
				exit(1);
			}
			for(uint32_t m = 0; m < kusage_clusters.size(); ++m){
				for(const auto& n : kusage_clusters[m]) {
					auto& gt = genotypes[n];
					gt.gt_k = m;
					gt.gt = acc_gt;
				}
				gt_reps.emplace_back(kusage_reps[m]);
				++acc_gt;
			}
		}
	}
	return acc_gt;
}
*/

int anallele_cluster(const double& max_error_l, const double& max_error_c, const std::vector<ANALLELE>& alleles, std::vector<Genotype>& genotypes, std::vector<int>& gt_reps)
{
	std::vector<int> allele_indeces(alleles.size());
	for(uint32_t i = 0; i < allele_indeces.size(); ++i) allele_indeces[i] = i;

	/** cluster by length **/
	std::vector<std::vector<int>> length_clusters;
	std::vector<int> length_reps;
	DistMatrix distmatrix_length(allele_indeces.size());
	anallele_cluster_length(max_error_l, alleles, allele_indeces, distmatrix_length, length_clusters, length_reps);
	if(length_reps.size() != length_clusters.size()){
		std::cerr << "[ERROR] unexpected number of representative alleles (" << gt_reps.size() << ") for " << length_clusters.size() << " length clusters" << std::endl;
		exit(1);
	}
	for(uint32_t i = 0; i < length_clusters.size(); ++i){
		const auto& length_cluster = length_clusters[i];
		for(const auto& j : length_cluster) genotypes[j].gt_l = i;
	}
	
	/** cluster by kusage **/
	std::vector<std::vector<int>> kusage_clusters;
	std::vector<int> kusage_reps;
	std::vector<KUSAGE> kusages;
	DistMatrix distmatrix_kusage(allele_indeces.size());
	anallele_cluster_kusage(max_error_c, 3, alleles, allele_indeces, distmatrix_kusage, kusages, kusage_clusters, kusage_reps);
	if(kusage_reps.size() != kusage_clusters.size()){
		std::cerr << "[ERROR] unexpected representative alleles (" << gt_reps.size() << ") for " << kusage_clusters.size() << " kusage clusters" << std::endl;
		exit(1);
	}
	for(uint32_t i = 0; i < kusage_clusters.size(); ++i){
		const auto& kusage = kusage_clusters[i];
		for(const auto& j : kusage) {
			genotypes[j].gt_k = i;
			genotypes[j].hsd = kusages[j].hsdiv();
		}
	}

	/** define final clusters **/
	std::list<int> remaining_alleles;
	for(int i = 0; i < (int)alleles.size(); ++i) remaining_alleles.push_back(i);
	std::vector<std::vector<int>> final_clusters;
	while(!remaining_alleles.empty()){
		int i = remaining_alleles.front();
		final_clusters.emplace_back();
		auto& local_cluster = final_clusters.back();
		auto it_j = remaining_alleles.begin();
		while(it_j != remaining_alleles.end()){
			if(genotypes[i].gt_l == genotypes[*it_j].gt_l && genotypes[i].gt_k == genotypes[*it_j].gt_k){
				local_cluster.emplace_back(*it_j);
				it_j = remaining_alleles.erase(it_j);
			}
			else ++it_j;
		}
	}
	for(int i = 0; i < (int)final_clusters.size(); ++i) {
		std::vector<uint32_t> tmpvec;
		for(const auto& j : final_clusters[i]) {
			genotypes[j].gt = i;
			tmpvec.emplace_back(j);
		}
		gt_reps.emplace_back((int)distmatrix_length.get_medoid(tmpvec));
	}
	
	return (int)final_clusters.size();
}


