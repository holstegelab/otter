#include "bindings/cpp/WFAligner.hpp"
#include "anseqs.hpp"
#include "anfahelper.hpp"
#include "andistmat.hpp"
#include <vector>
#include <string>

/**
 * @brief  Attempt to realign flanks around region of interest (ROI) in non-spanning reads due to clipping
 * 
 * Fetches corresponding subseq of size *flank* in the ROI in the reference, and aligns it 
 * to one or both ends of non-spanning reads. A read becomes spanning if the flanking alignments 
 * have a minimum seq similarity of *min_sim*, and the read seq is adjusted accordingly.
 * 
 * @param  chr:  Chromosome name of ROI
 * @param  start:  Start coord of ROI
 * @param  end:  End coord of ROI
 * @param  flank:  Size of flanking seq to align
 * @param  min_sim: Minimum seq similarity
 * @param  faidx_inst: FaidxInstance object of the ref genome
 * @param  reads: vector of ANREAD objects
 * @param  aligner: WFAligner object
 * 
 * @return void
 */
void local_realignment(
	const std::string& chr, 
	const int& start, 
	const int& end, 
	const int& flank, 
	const double& min_sim, 
	const FaidxInstance& faidx_inst, 
	std::vector<ANREAD>& reads, 
	wfa::WFAligner& aligner
);

double get_dist_anreads(
    const bool& ignore_haps, 
    wfa::WFAligner& aligner, 
    ANREAD& read_x, 
    ANREAD& read_y
);

void fill_dist_matrix(
    const bool& ignore_haps, 
    wfa::WFAligner& aligner, 
    std::vector<ANREAD>& reads,
    std::vector<int>& indeces, 
    DistMatrix& distmatrix
);

void invalid_reassignment(
    const bool& ignore_haps, 
    const double& min_sim, 
    const double& max_error, 
    const int& total_alleles, 
    std::vector<ANREAD>& reads, 
    wfa::WFAligner& aligner, 
    std::vector<int>& labels
);

void rapid_consensus(
    const bool& ignore_haps,
    std::vector<ANREAD>& reads, 
    std::vector<int>& labels, 
    std::vector<int>& valid_indeces, 
    int& total_alleles, 
    DistMatrix& valid_distmatrix, 
    wfa::WFAligner& aligner,
    std::vector<ANALLELE>& alleles
);
