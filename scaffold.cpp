#include <vector>
#include <algorithm>

#include "scaffold.h"

using std::vector;
using std::ofstream;
using std::endl;

Scaffold::Scaffold(Contig *contig) {
	contigs.emplace_back(contig);
	id = contig->get_id();
}

vector<Contig*>& Scaffold::get_contigs() {
	return contigs;
}

void Scaffold::merge(Scaffold *scaffold) {
	for (auto c : scaffold->contigs) {
		c->set_scaffold_id(id);
	}
	contigs.insert(contigs.end(), scaffold->get_contigs().begin(), scaffold->get_contigs().end());
}

void Scaffold::reverse() {
	std::reverse(contigs.begin(), contigs.end());
	for (auto c : contigs) {
		c->reverse();
	}
}

string Scaffold::get_merged() {
	string output = "";
	for(auto c : contigs) {
		if(!c->is_reversed())
			output += c->get_data();
		else
			output += c->get_reverse_complement();
	}

	return output;
}