#ifndef SCAFFOLD_H
#define SCAFFOLD_H

#include <vector>
#include <string>
#include <iostream>

#include "contig.h"

using std::vector;
using std::string;
using std::ofstream;

class Scaffold {

public:

	Scaffold(Contig *contig);

	string get_id() {
		return id;
	}

	void merge(Scaffold *scaffold);

	vector<Contig*>& get_contigs();

	void reverse();

	string get_merged();

private:

	string id;
	vector<Contig*> contigs;
};


#endif