#ifndef ANPPOA_HPP
#define ANPPOA_HPP

#include <string>
#include <vector>
#include <list>
#include <memory>
#include <map>
#include <set>
#include <iostream>

struct ppoa_edge {
	uint32_t source;
	uint32_t sink;
	float weight;
	ppoa_edge(uint32_t _source, uint32_t _sink, float _weight = 1.0f): source(_source), sink(_sink), weight(_weight){};
};

struct ppoa_path{
	float weight;
	std::vector<uint32_t> path;
	ppoa_path(){};
	ppoa_path(float w): weight(w){};
	ppoa_path(float w, std::vector<uint32_t>& p): weight(w), path(p){};
};

class PPOA {

	public:
		std::string backbone;
		std::vector<bool> hps;
		//vector of indeces each representing a node in the backbone seq
		std::vector<std::string> nodes;
		//vector of edges; each index is corresponds to outgoing edges
		std::vector<std::vector<ppoa_edge>> edges;
		std::map<uint32_t, std::vector<uint32_t>> alt_nodes;

		PPOA(){};
		PPOA(std::string&);
		void init(std::string&);

		void insert_node(uint32_t, const std::string&);
		void insert_edge(uint32_t, uint32_t);
		void insert_alignment(std::string&, std::string&, bool&, bool&);

		void adjust_weights(float, float);
		void consensus(std::string&);

		void printDOT();
		void clear();


	private:
		uint32_t last_id = 0;
		std::vector<uint32_t> starting_nodes;
		std::set<uint32_t> ending_nodes;
		std::set<uint32_t> consensus_starting_nodes;
		void set_heaviest_path(const uint32_t&, std::map<uint32_t, ppoa_path>&);

};

inline PPOA::PPOA(std::string& _b){init(_b);};

inline void PPOA::init(std::string& _b)
{
	backbone = _b;
	hps.resize(backbone.size(), false);
	nodes.resize(backbone.size());
	edges.resize(backbone.size());
	last_id = backbone.size();
	for(uint32_t i = 1; i < backbone.size(); ++i) {
		if(i == 1) {
			insert_node(i-1, std::string{backbone[i-1]});
			starting_nodes.emplace_back(i-1);
		}
		insert_node(i, std::string{backbone[i]});
		insert_edge(i-1, i);
		if(backbone[i] == backbone[i-1]){
			hps[i] = true;
			if(!hps[i-1]) hps[i-1] = true;
		}
	}
};

inline void PPOA::insert_node(uint32_t id, const std::string& seq)
{
	if(id < last_id) nodes[id] = seq;
	else{
		nodes.emplace_back(seq);
		edges.emplace_back();
		last_id = id + 1;
	}
};

inline void PPOA::insert_edge(uint32_t source, uint32_t sink) {
	auto& local_edges = edges[source];
	if(local_edges.empty()) local_edges.emplace_back(source, sink);
	else{
		bool found = false;
		for(auto& edge : local_edges){
			if(edge.sink == sink){
				edge.weight += 1.0;
				found = true;
				break;
			}
		}
		if(!found) local_edges.emplace_back(source, sink);
	}
};

inline void PPOA::insert_alignment(std::string& sequence, std::string& cigar, bool& is_spanning_l, bool& is_spanning_r) {
	int previous_node = 0, ref_i = 0, target_i = 0, cigar_i = 0;
	bool is_first_node = true;

	if(!is_spanning_l){
		is_first_node = false;
		while (cigar_i < cigar.size()){
			const auto& c = cigar[cigar_i];
			if(c != 'D' && c != 'I') break;
			if(c == 'D') {
				++ref_i;
				previous_node = ref_i;
			}
			else if(c == 'I') ++target_i;
			++cigar_i;
		}

		//std::cout << previous_node << ' ' << ref_i << ' ' << target_i << ' ' << cigar_i << '\n';
	}

	while(cigar_i < cigar.size()){//for(const auto& c : cigar){
		//std::cout << previous_node << ' ' << ref_i << ' ' << target_i << ' ' << cigar_i << '\n';
		const auto& c = cigar[cigar_i];
		//path traversal remains in the backbone
		std::string target_seq = std::string{sequence[target_i]};
		//std::cout << target_seq << '\n';
		if(c == 'M' || c == 'X'){
			if(c == 'M'){
				if(is_first_node || previous_node == ref_i) is_first_node = false;
				else insert_edge(previous_node, ref_i);
				previous_node = ref_i;
			}
			//new node to add
			else if(c == 'X'){
				//new start to graph
				if(is_first_node) {
					bool need_new_starting_node = true;
					for(const auto& node : starting_nodes){
						if(nodes[node] == target_seq){
							need_new_starting_node = false;
							break;
						}
					}
					if(need_new_starting_node){
						insert_node(last_id, target_seq);
						previous_node = last_id - 1;
						starting_nodes.emplace_back(previous_node);
					}
					is_first_node = false;
				}
				//check if alternate edges in backbone already exist
				else{
					//check for existing matching alternate edge
					auto& outgoing_edges = edges[previous_node];
					int matching_outgoing_edges_i = -1;
					for(int i = 0; i < (int)outgoing_edges.size(); ++i){
						if(nodes[outgoing_edges[i].sink] == target_seq && outgoing_edges[i].sink >= backbone.size()){
							matching_outgoing_edges_i = i;
							break;
						}
					}
					//found existing alternate edge
					if(matching_outgoing_edges_i >= 0) {
						auto& matching_edge = outgoing_edges[matching_outgoing_edges_i];
						++matching_edge.weight;
						previous_node = matching_edge.sink;
					}
					//adding new alternate edge
					else{
						uint32_t new_node = last_id;
						insert_node(new_node, std::string{sequence[target_i]});
						insert_edge(previous_node, new_node);
						previous_node = new_node;
					}
				}
			}
			++ref_i;
			++target_i;
		}
		if(c == 'D'){
			if(!is_first_node) ++ref_i;
			else {
				//is_first_node = false;
				++ref_i;
				previous_node = ref_i;
			}
		}
		else if(c == 'I'){
			if(is_first_node){
				insert_node(last_id, target_seq);
				previous_node = last_id - 1;
				starting_nodes.emplace_back(previous_node);
				is_first_node = false;
			}
			else{
				//check for existing matching alternate edge
				auto& outgoing_edges = edges[previous_node];
				int matching_outgoing_edges_i = -1;
				for(int i = 0; i < (int)outgoing_edges.size(); ++i){
					const auto& outgoing_edge = outgoing_edges[i];
					//insertion, hence it cannot be an outgoing edge to backbone
					if(outgoing_edge.sink >= backbone.size() && nodes[outgoing_edge.sink] == target_seq){
						matching_outgoing_edges_i = i;
						break;
					}
				}
				//found existing alternate edge
				if(matching_outgoing_edges_i >= 0) {
					auto& matching_edge = outgoing_edges[matching_outgoing_edges_i];
					++matching_edge.weight;
					//std::cout << previous_node << '\t' << matching_edge.sink << '\n';
					previous_node = matching_edge.sink;
				}
				//adding new alternate edge
				else{
					uint32_t new_node = last_id;
					insert_node(new_node, std::string{sequence[target_i]});
					insert_edge(previous_node, new_node);
					//std::cout << previous_node << '\t' << new_node << '\n';
					previous_node = new_node;
				}
			}
			++target_i;
		}
		//if(ref_i <= 10 && is_spanning_l) consensus_starting_nodes.insert(previous_node);
		if((int)backbone.size() - ref_i <= 10 && is_spanning_r) ending_nodes.insert(previous_node);
		++cigar_i;
	}
	//std::cout << "last node: " << previous_node << '\n';
}

inline void PPOA::adjust_weights(float c, float t)
{
	for(auto& local_edges : edges){
		for(auto& e : local_edges){
			float t_applied = t * e.weight;
			float final_weight = c > t_applied ? c : t_applied;
			e.weight = e.weight - final_weight;
		}
	}
}

inline void PPOA::set_heaviest_path(const uint32_t& node, std::map<uint32_t, ppoa_path>& map_heaviest)
{
	if(map_heaviest.find(node) == map_heaviest.end()){
		//get incoming edges
		std::vector<ppoa_edge> incoming;
		for(const auto& local_edges : edges){
			for(const auto& e : local_edges){
				if(e.sink == node) incoming.emplace_back(e);
			}
		}
		//set default weight and path
		bool not_h_defined = true;
		float h_weight = 0.0f;
		//current node is global source node, no weight/path
		if(incoming.empty()) map_heaviest.insert({node, ppoa_path(h_weight)});
		else{
			std::vector<uint32_t> h_path;
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
			map_heaviest.insert({node, ppoa_path(h_weight, h_path)});
		}
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

inline void PPOA::consensus(std::string& consensus_seq)
{
	std::map<uint32_t, ppoa_path> map_heaviest;
	std::list<uint32_t> remaining_nodes;
	for(uint32_t n = 0; n < nodes.size(); ++n) remaining_nodes.emplace_back(n);

	std::vector<std::vector<ppoa_edge>> incoming_edges(nodes.size());
	for(const auto& local_edges : edges){
		for(const auto& e : local_edges) incoming_edges[e.sink].emplace_back(e);
	}

	while(!remaining_nodes.empty()){
		uint32_t next = remaining_nodes.front();
		bool incoming_defined = true;
		for(const auto& e : incoming_edges[next]) {
			if(map_heaviest.find(e.source) == map_heaviest.end()) {
				incoming_defined = false;
				break;
			}
		}
		if(incoming_defined) {
			set_heaviest_path(next, map_heaviest);
			remaining_nodes.pop_front();
		}
		else {
			remaining_nodes.pop_front();
			remaining_nodes.emplace_back(next);
		}
	}

	/**
	std::map<uint32_t, ppoa_path> map_heaviest;
	for(uint32_t n = 0; n < nodes.size(); ++n) {
		std::cout << n << std::endl;
		set_heaviest_path(n, map_heaviest);
	}
	*/
	uint32_t h_node = 0;
	ppoa_path h_path;
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


	//for(uint32_t i = 0; i < nodes.size(); ++i) std::cout << i << '\t' << nodes[i] << '\n';

	h_path.path.emplace_back(h_node);

	for(uint32_t i = 0; i < h_path.path.size(); ++i){
		auto& seq = nodes[h_path.path[i]];
		consensus_seq += seq;
	}

}

inline void PPOA::printDOT()
{
	std::cout << "digraph ansparc {\n  graph [rankdir = LR]\n";
	for(uint32_t node_id = 0; node_id < last_id; ++node_id) std::cout << "  " << node_id << "[label = \"" << node_id << "-" << nodes[node_id] << "\"]\n";
	for(int source = 0; source < (int)nodes.size(); ++source) {
		const auto& local_edges = edges[source];
		for(const auto edge: local_edges) std::cout << "  " << source << " -> " << edge.sink << " [label = \"" << edge.weight << "\"]\n";
	}
	std::cout << "}\n";
}

#endif