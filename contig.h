#ifndef CONTIG_H
#define CONTIG_H

#include <string>

using std::string;


class Contig {

public:
	
	Contig(string& _id, string& _data);

	string& get_id() {
		return id;
	}

	string& get_data() {
		return data;		
	}

	string get_reverse_complement();

	int length() {
		return len;
	}

	string& get_scaffold_id() {
		return scaffold_id;
	}

	bool get_start() {
		return start;
	}

	bool get_end() {
		return end;
	}

	bool is_reversed() {
		return reversed;
	}

	void reverse();

	void set_scaffold_id(string _scaffold_id);

	void set_end(bool _end);

	void set_start(bool _start);

private:

	string id, data, scaffold_id;
	int len;
	bool reversed, start, end;
};


#endif