#ifndef ANBAMTOOLS_HPP
#define ANBAMTOOLS_HPP

#include "parse_bam_alignments.hpp"
#include "antimestamp.hpp"
#include "anbamfilehelper.hpp"
#include "otter_opts.hpp"
#include "anbed.hpp"
#include <htslib/sam.h>
#include <string>
#include <memory>
#include <vector>

/**
 * Struct that accompanies the 'abg_generate' function (see below), for determining both whether a 
 * an AB-graph generation was successful, and the spanning-type.
 * 
 * @param bool    Whether an AB-graph was successfully generated
 * @param bool    Whether the corresponding read spans left-breakpoint
 * @param bool    Whether the corresponding read spans right-breakpoint
 */
class ParsingStatus {
    public:
    	bool successful;
    	bool spanning_l;
    	bool spanning_r;
        std::pair<int, int> alignment_coords;
    	ParsingStatus();
    	bool is_spanning() const;
};

class Haplotag{
    public:
        int hp;
        int ps;
        Haplotag();
        Haplotag(int, int);
        bool operator==(const Haplotag &x) const;
        bool is_defined() const;
};

class AlignmentBlock {
    public:
        std::vector<std::string> names;
        std::vector<std::string> seqs;
        std::vector<Haplotag> hps;
        std::vector<ParsingStatus> statuses;

        bool n_spannable(int) const;
};

/**
 * Function to project reference-coordinates on a given read. The projection is stored in the given 
 * vector such that each nucleotide in the read has a corresponding reference coordinate. A value of
 * '-1' indicates that the nucleotide does not exist in the reference. Additionally, the given pair 
 * of bool variables are changed to reflect if the read-alignment is clipped.
 * 
 * @param bam1_t       Read-alignment
 * @param bool         Read is left-clipped
 * @param bool         Read is right-clipped
 * @param vector<int>  Projection of reference-coordinatees on the read  
 */
void project_positions(
    bam1_t* alignment, 
    bool& clipped_l,
    bool& clipped_r,
    std::vector<int>& refcoords
);

    /**
 * Function to find corresponding start/end coordinates (breakpoints) of a given read given clipping
 * information and the projected reference coordinates. Breakpoints as stored in the given pointer 
 * of int-pairs. Failure to find the breakpoints as well as the spanning-type are reflected in the
 * given 'abg_generate_msg' variable.
 * 
 * @param int                         Start coordinate of breakpoint
 * @param int                         End coordinate of breakpoint
 * @param ParsingStatus               Status message
 * @param vector<int>                 Projected positions
 * @param unique_ptr<pair<int,int>>   Variable to store corresponding breakpoints
 * 
 */
void get_breakpoints(
    const int& start, 
    const int& end,
    ParsingStatus& msg,
    std::vector<int>& projectedcoords,
    std::unique_ptr<std::pair<int, int>>& ptr
);

void parse_alignment(
    const int& rstart, 
    const int& rend, 
    bam1_t* read, 
    ParsingStatus& msg, 
    std::string& seq
);

void parse_alignments(
    const OtterOpts& params,
    const BED& bed_region, 
    const BamInstance& bam_inst,
    AlignmentBlock& alignment_block
);

#endif