/**
 * An implementation of the sparc algorithm to rapidly compute a consensus sequence given 
 * a sequence 'backbone' and one or more sampled sequences.
 * 
 * See https://peerj.com/articles/2016/ for details.
 * 
 * Last updated: 2023/07/03
 */

#ifndef ANSPARC_HPP
#define ANSPARC_HPP

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <set>

/** ##### SPARC graph  class definitions #####*/

struct ansparc_node {
	uint32_t id;
	std::string seq;
	ansparc_node(uint32_t& _id, std::string& _seq): id(_id), seq(_seq){};
};

struct ansparc_edge {
	uint32_t source;
	uint32_t sink;
	float weight;
	ansparc_edge(uint32_t& _source, uint32_t& _sink, float _weight = 1.0f): source(_source), sink(_sink), weight(_weight){};
};

struct ansparc_path{
	float weight;
	std::vector<uint32_t> path;
	ansparc_path(){};
	ansparc_path(float& w, std::vector<uint32_t>& p): weight(w), path(p){};
};

class ansparc_graph {

	public:
		std::string backbone;
		uint32_t k;
		uint32_t g;
		uint32_t l;
		//vector of indeces each representing a possible node in the backbone seq
		std::map<uint32_t, std::string> nodes;
		//vector of edges; each index is corresponds to outgoing edges
		std::vector<ansparc_edge> edges;
		std::map<uint32_t, std::vector<uint32_t>> alt_nodes;

		ansparc_graph(){};
		ansparc_graph(std::string&, uint32_t&, uint32_t&, uint32_t&);
		void init(std::string&, uint32_t&, uint32_t&, uint32_t&);

		void insert_node(uint32_t&, std::string&);
		void insert_edge(uint32_t&, uint32_t&);
		void insert_alignment(std::string&, std::string&);

		void adjust_weights(float&, float&);
		void consensus(std::string&);

		void printDOT();
		void clear();


	private:
		uint32_t last_id = 0;
		std::vector<uint32_t> backbone_pos;
		std::vector<std::pair<uint32_t, uint32_t>> alts;
		std::set<uint32_t> ending_nodes;
		void set_heaviest_path(const uint32_t&, std::map<uint32_t, ansparc_path>&);
		uint32_t true_seq_size(std::string&) const;

};

inline ansparc_graph::ansparc_graph(std::string& _b, uint32_t& _k, uint32_t& _g, uint32_t& _l){init(_b, _k, _g, _l);};

inline void ansparc_graph::init(std::string& _b, uint32_t& _k, uint32_t& _g, uint32_t& _l)
{
	backbone = _b;
	k = _k;
	g = _g;
	l = _l;
	uint32_t pos = 0;
	uint32_t last_pos = 0;
	std::string kmer;

	while(pos < backbone.size()){
		int d = (pos + k) - backbone.size();
		if(d <= 0) kmer = backbone.substr(pos, k); //pos + k <= backbone.size()
		else {
			kmer = '-';
			for(int i = 1; i < d; ++i) kmer += '-';
			kmer += backbone.substr(pos);
		}
		insert_node(pos, kmer);
		backbone_pos.emplace_back(pos);
		if(pos > 0) {
			insert_edge(last_pos, pos);
			last_pos = pos;
		}
		pos += g;
	}
};

inline void ansparc_graph::insert_node(uint32_t& id, std::string& seq)
{
	if(id < last_id) {
		std::cerr << "ERROR: unexpected node ID when inserting: input=" << id << "expeced=" << last_id << std::endl;
		exit(1);
	}
	nodes.insert({id, seq});
	last_id = id + 1;
};

inline void ansparc_graph::insert_edge(uint32_t& source, uint32_t& sink) {
	auto it = edges.begin();
	while(it != edges.end()){
		if(it->source != source || it->sink != sink) ++it;
		else {
			it->weight += 1.0f;
			return;
		}
	}
	edges.emplace_back(source, sink);
};

inline void ansparc_graph::insert_alignment(std::string& cigar, std::string& seq)
{
	//project corresponding coordinates of cigar to bakebone (b) and seq (s)
	std::vector<std::unique_ptr<uint32_t>> b_matched(cigar.size());
	std::vector<std::unique_ptr<uint32_t>> s_matched(cigar.size());

	//iterate and set corresponding index per sequence
	uint32_t b_i = 0;
	uint32_t s_i = 0;
	for(uint32_t i = 0; i < cigar.size(); ++i){
		auto& op = cigar[i];
		if(op == 'M' || op == 'X') {
			b_matched[i] = std::make_unique<uint32_t>(b_i);
			++b_i;
			s_matched[i] = std::make_unique<uint32_t>(s_i);
			++s_i;
		}
		else if(op == 'I'){
			s_matched[i] = std::make_unique<uint32_t>(s_i);
			++s_i;
		}
		else if(op == 'D'){
			b_matched[i] = std::make_unique<uint32_t>(b_i);
			++b_i;
		}
	}


	//iterate through each kmer in backbone and obtain corresponding kmer in seq
	uint32_t last_cigar_i = 0;
	std::vector<uint32_t> captured_nodes;
	for(uint32_t b_start_i = 0; b_start_i < backbone_pos.size(); ++b_start_i){
		uint32_t b_start = backbone_pos[b_start_i];
		std::string s_seq;
		//set end pos for current kmer
		uint32_t b_end = b_start + k;
		//identify corresponding cigar indeces for kmer
		uint32_t c_start = 0;
		uint32_t c_end = 0;
		for(uint32_t cigar_i = last_cigar_i; cigar_i < cigar.size(); ++cigar_i){
			if(b_matched[cigar_i] != nullptr && *b_matched[cigar_i] <= b_start) {
				c_start = cigar_i;
				++last_cigar_i;
			}
			if(b_matched[cigar_i] != nullptr && *b_matched[cigar_i] < b_end) c_end = cigar_i;
		}		
		if(b_start + 1 >= backbone.size() && cigar.back() == 'I') ++c_end;
		uint32_t inserted_node;
		//identify corresponding indeces in seq
		if(s_matched[c_start] != nullptr || s_matched[c_end] != nullptr) {
			//adjust starting index for seq
			if(s_matched[c_start] == nullptr) {
				uint32_t z = c_start;
				while(z <= c_end) {
					if(s_matched[z] != nullptr) s_seq += seq[*s_matched[z]];
					else s_seq += '-';
					++z;
				}
			}
			//set seq to end
			else if(s_matched[c_end] == nullptr) {
				for(uint32_t i = c_start; i < c_end; ++i) {
					if(cigar[i] != 'D') s_seq += '-'; else s_seq += seq[*s_matched[i]];
				}
				/**
				 * 				uint32_t s_pos = *s_matched[c_start];
				std::cout << "here " << s_pos << ' ' << seq.substr(c_start, c_end - c_start + 1) << '\n';
				if(s_pos + 1 < seq.size()) {
					s_seq += seq[s_pos];
					for(uint32_t i = 1; i < k; ++i) s_seq += '-';
				}
				else {
					std::cerr << "ERROR: unexpected start/end coordinates: " << s_pos << ',' << seq.size() << '\n';
					exit(1);
				}
				&*/
			}
			//set seq accordingly
			else {		
				uint32_t length = 1 + *s_matched[c_end]  - *s_matched[c_start];
				if(length < k) for(uint32_t i = 0; i < (k - length); ++i) s_seq += '-';
				s_seq += seq.substr(*s_matched[c_start], length);
			}
			auto it = nodes.find(b_start);
			if(it->second == s_seq) {
				captured_nodes.emplace_back(b_start);
				inserted_node = b_start;
			}
			//alt node required
			else {
				auto it_a = alt_nodes.find(b_start);
				uint32_t new_node_id = last_id;
				//no available alt nodes based on starting pos, add node
				if(it_a == alt_nodes.end()) {
					alt_nodes.insert({b_start, {new_node_id}});
					insert_node(new_node_id, s_seq);
				}
				//potential alt nodes exist
				else {
					bool missing = true;
					for(const auto& a : it_a->second){
						if(nodes.find(a)->second == s_seq) {
							//set new node id
							new_node_id = a;
							missing = false;
							break;
						}
					}
					if(missing) insert_node(new_node_id, s_seq);
				}
				captured_nodes.emplace_back(new_node_id);
				inserted_node = new_node_id;
			}
			s_seq.clear();
		}

		if(backbone_pos.size() - b_start_i <= l) ending_nodes.insert(inserted_node);
	}

	for(uint32_t j = 1; j < captured_nodes.size(); ++j) insert_edge(captured_nodes[j-1], captured_nodes[j]);
}

inline void ansparc_graph::adjust_weights(float& c, float& t)
{
	for(auto& e : edges){
		float t_applied = t * e.weight;
		float final_weight = c > t_applied ? c : t_applied;
		e.weight = e.weight - final_weight;
	}
}

/** ### PSEUDOCODE FOR CONSENSUS SEQUENCE ###

path: 2-tuple {weight, path} //define path struct
H <- {} //map of node -> path

def heaviest_path(node): //function to recursively backtrack and find heaviest path for a given node
	I <- incoming_nodes(node) // incoming edges
	h_weight <- 0 //init heaviest weight
	h_path <- {} //init heaviest path
	for all i in I: //itereate through every incoming node
		if H[i] = null: //traverse incoming node first
			h = heaviest_path(i) // recursive call
			H.insert{i, h.weight, h.path} //insert heaviest info of incoming node
		if H[i].weight > heaviest_weight: //incoming node is so far the heaviest
			h_weight = H[i].weight
			h_path = H[i].path
	return {h_weight, h_path} //return heaviest 

for each node: //iterate and set heavieast path for all nodes
	if(H[n] == {}): //heaviest path not yet init
		h = heaviest_path(node)
		H.insert{n, h.weight, h.path}

return node/path with highest weight in H
*/

inline void ansparc_graph::consensus(std::string& consensus_seq)
{
	std::map<uint32_t, ansparc_path> map_heaviest;
	for(const auto& n : nodes) set_heaviest_path(n.first, map_heaviest);

	uint32_t h_node;
	ansparc_path h_path;
	bool not_init = true;

	for(const auto& h : map_heaviest) {
		if(ending_nodes.find(h.first) != ending_nodes.end()){
			if(not_init || h.second.weight > h_path.weight){
				if(not_init) not_init = false;
				h_node = h.first;
				h_path = h.second;
			}
		}
	}

	h_path.path.emplace_back(h_node);

	for(uint32_t i = 0; i < h_path.path.size() - 1; ++i){
		//std::cout << consensus_seq << '\n';
		auto& seq = nodes.find(h_path.path[i])->second;
		uint32_t seq_l = true_seq_size(seq);
		if(seq_l > k) consensus_seq += seq.substr(0, seq_l - k + g);
		else for(uint32_t j = 0; j < g; ++j) if(seq[j] != '-') consensus_seq += seq[j];
	}
	/**
	for(int i = 1; i < h_path.path.size(); ++i){
		auto& seq_last = nodes.find(h_path.path[i-1])->second;
		uint32_t seq_last_l = true_seq_size(seq_last);
		auto& seq = nodes.find(h_path.path[i])->second;
		uint32_t seq_l = true_seq_size(seq);
		std::cout << seq_last << "->" << seq << '\n';
		if(seq_last_l < k) {
			if(seq.size() < k) consensus_seq += seq_last;
			else consensus_seq += seq_last.substr(1);
		}
		else{
			if(seq.size() < k) consensus_seq += seq_last;
			else if(seq_last.size() == k) consensus_seq += seq_last.substr(0, g);
			else consensus_seq += seq_last.substr(0, seq_last.size() - k + g);	
		}
	}
	*/

	auto& last =  nodes.find(h_path.path.back())->second;
	uint32_t last_l = true_seq_size(last);
	if(last_l >= k) consensus_seq += last.substr(k - g);
	else for(uint32_t j = 0; j < last.size(); ++j) if(last[j] != '-') consensus_seq += last[j];
}

inline void ansparc_graph::set_heaviest_path(const uint32_t& node, std::map<uint32_t, ansparc_path>& map_heaviest)
{
	if(map_heaviest.find(node) == map_heaviest.end()){
		//get incoming edges
		std::vector<ansparc_edge> incoming;
		for(const auto& e : edges){
			if(e.sink == node) incoming.push_back(e);
		}
		//set default weight and path
		bool not_h_defined = true;
		float h_weight = 0.0f;
		std::vector<uint32_t> h_path;
		//current node is global source node, no weight/path
		if(incoming.empty()) map_heaviest.insert({node, ansparc_path(h_weight, h_path)});
		else{
			for(const auto& e : incoming){
				auto it = map_heaviest.find(e.source);
				//traverse incoming node first
				if(it == map_heaviest.end()){
					set_heaviest_path(e.source, map_heaviest);
					it = map_heaviest.find(e.source);
				}
				if(not_h_defined || it->second.weight + e.weight > h_weight){
					if(not_h_defined) not_h_defined = false;
					h_weight = it->second.weight + e.weight;
					h_path = it->second.path;
					h_path.emplace_back(e.source);
				}
			}
			map_heaviest.insert({node, ansparc_path(h_weight, h_path)});
		}
	}
}

inline uint32_t ansparc_graph::true_seq_size(std::string& seq) const
{
	uint32_t l = 0;
	for(const auto& c : seq) if(c != '-') ++l;
	return l;
}


inline void ansparc_graph::printDOT()
{
	std::cout << "digraph ansparc {\n  graph [rankdir = LR]\n";
	for(const auto& n : nodes) std::cout << "  " << n.first << "[label = \"" << n.first << "-" << n.second << "\"]\n";
	for(const auto& e : edges) std::cout << "  " << e.source << " -> " << e.sink << " [label = \"" << e.weight << "\"]\n";
	std::cout << "}\n";
}

inline void ansparc_graph::clear()
{
	backbone.clear();
	nodes.clear();
	edges.clear();
	alt_nodes.clear();
	backbone_pos.clear();
	alts.clear();
	last_id = 0;
}


#endif