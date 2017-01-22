#include <string>
#include <algorithm>
#include <iostream>

#include "contig.h"

using std::string;

Contig::Contig(string& _id, string& _data) {
	id = _id;
	data = _data;
	scaffold_id = _id;
	reversed = false;
	start = end = true;
}

void Contig::set_scaffold_id(string _scaffold_id) {
	scaffold_id = _scaffold_id;
}

void Contig::reverse() {
	reversed = !reversed;
}

void Contig::set_start(bool _start) {
	start = _start;
}

void Contig::set_end(bool _end) {
	end = _end;
}

string Contig::get_reverse_complement() {
	string c = "";
	for (auto base : data) {
		if (base == 'A') c.push_back('T');
		if (base == 'T') c.push_back('A');
		if (base == 'C') c.push_back('G');
		if (base == 'G') c.push_back('C');
	}
	std::reverse(c.begin(), c.end());
	
	return c;
}