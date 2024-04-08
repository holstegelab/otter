#include "formatter.hpp"
#include <string>
#include <iostream>
#include "otter_opts.hpp"
#include "parse_bam_alignments.hpp"
#include "antimestamp.hpp"

void output_fa2sam(const std::string& read, const std::string& chr, const int& start, const int& end, const std::string& seq, const std::string& rg, const int& acov, const int& tc, const int& spanning_l, const int& spanning_r, const int& initial_clusters, const double& se, std::vector<Haplotag>& hps)
{
	std::string pseudo_qual(seq.size(), '!');
	if(!rg.empty()) pseudo_qual += "\tRG:Z:" + rg;
	std::cout << read << "\t0\t" << chr << '\t' << start << "\t0\t" << seq.size() << "M\t*\t0\t0\t" << seq << '\t' << pseudo_qual << '\t' << "ta:Z:" << chr << ':' << start << '-' << end;
	if(acov > 0) std::cout << "\tac:i:" << acov;
	std::cout << "\ttc:i:" << tc;
	if(spanning_l >= 0) {
		char sp;
		if(spanning_l && spanning_r) sp = 'b';
		else if(spanning_l) sp = 'l';
		else if(spanning_r) sp = 'r';
		else sp = 'n';
		std::cout << "\tsp:A:" << sp;
	}
	if(initial_clusters > -1) std::cout << "\tic:i:" << initial_clusters;
	std::cout << "\tse:f:" << se;
    if(hps.size() > 0 && hps[0].is_defined())
    {
        std::cout << "\tPS:i:" << hps[0].ps;
        std::cout << "\tHP:i:" << hps[0].hp;
    }
    //in some situations, an assembled read connects two phase sets
    if(hps.size() > 1 && hps[1].is_defined())
    {
        std::cout << "\txs:i:" << hps[1].ps;
        std::cout << "\txp:i:" << hps[1].hp;
    }
   
	std::cout << '\n';
}
