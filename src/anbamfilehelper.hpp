#ifndef ANBAMFILEHELPER_HPP
#define ANBAMFILEHELPER_HPP

#include "sam.h"
#include <vector>
#include <string>
#include <map>
#include <memory>

const char init_header[] = "@HD\tVN:1.4\tSO:unknown\n";
const uint8_t base_qual = (uint8_t)0;

/**
 * Self-contained struct for (independently) loading a BAM file (i.e. multi-threading). It contains an
 * 'init' and 'destroy' methods for cleaning.
 * 
 * @param samFile     SAM/BAM-file file pointer
 * @param bam_hdr_t   SAM/BAM header
 * @param bool        Whether SAM/BAM is indexed
 * @param hts_idx_t   SAM/BAM header pointer
 * @param bam1_t      SAM/BAM read-alignment pointer
 */
class BamInstance {
	public:
		BGZF *fp;
		bam_hdr_t *header;
		bool index;
		hts_idx_t *idx;
		bam1_t* read;

		//ANSBAM();
		void init(const std::string&,bool);
		void destroy();
};


#endif