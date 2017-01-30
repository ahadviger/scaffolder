#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <unordered_map>
#include <map>

#include "contig.h"
#include "scaffold.h"

#define READ_THRESHOLD 0.1
#define CONTIG_THRESHOLD 0.4
#define MATCH_THRESHOLD 0.7

using std::string;
using std::ifstream;
using std::ofstream;
using std::getline;
using std::istringstream;
using std::sort;
using std::unordered_map;
using std::set;
using std::pair;
using std::cout;
using std::endl;

struct overlap {
	string read_id;
	int read_length, read_start, read_end;
	string strand;
	string contig_id;
	int contig_length, contig_start, contig_end;
	int num_matches, match_length;

	void print() {
		cout << read_id << ' ' << read_length << ' ' << read_start << ' ' << read_end << " | " 
			<< contig_id << ' ' << contig_length << ' ' << contig_start << ' ' << contig_end << " | "
			<< num_matches << ' ' << match_length << ' ' << strand << endl;
	}
};

vector<overlap> overlaps;
std::map<pair<string, string>, int> votes_end_start, votes_end_end, votes_start_start;

struct joint_overlap {
	string left, right;
	int votes, type;

	joint_overlap() {}
	
	joint_overlap(string _left, string _right, int _type) {
		left = _left; right = _right;
		type = _type;

		if (type == 0) votes = votes_end_start[std::make_pair(left, right)];
		else if (type == 1) votes = votes_end_end[std::make_pair(left, right)];
		else if (type == 2) votes = votes_start_start[std::make_pair(left, right)];
	}
};

inline bool cmp(const joint_overlap& x, const joint_overlap& y) {
	if (x.votes != y.votes) return x.votes > y.votes;
	if (x.left != y.left) return x.left < y.left;
	return x.right < y.right;
}

vector<Contig*> contigs;
set<Scaffold*> scaffolds;
vector<joint_overlap> joint_overlaps;
unordered_map<string, Contig*> id_to_contig;
unordered_map<string, Scaffold*> id_to_scaffold;


bool check_overlaps_end_start(overlap left, overlap right) {
	if (left.contig_id == right.contig_id) return false;
	if (left.strand != right.strand) return false;

	if (left.strand == "+" && left.read_start > right.read_start) return false;
	if (left.strand == "+" && left.read_end > right.read_end) return false;
	if (left.strand == "-" && right.read_start > left.read_start) return false;
	if (left.strand == "-" && right.read_end > left.read_end) return false;

	if (left.strand == "+" && right.contig_start + left.contig_length - left.contig_end > left.read_length) return false;
	if (left.strand == "-" && left.contig_start + right.contig_length - right.contig_end > left.read_length) return false;

	if ((double) left.contig_end / left.contig_length < 1 - CONTIG_THRESHOLD) return false;
	if ((double) right.contig_start / right.contig_length > CONTIG_THRESHOLD) return false;
	
	if (left.strand == "+" && (double) left.read_start / left.read_length > READ_THRESHOLD) return false;
	if (left.strand == "+" && (double) right.read_end / right.read_length < 1 - READ_THRESHOLD) return false;
	if (left.strand == "-" && (double) right.read_start / right.read_length > READ_THRESHOLD) return false;
	if (left.strand == "-" && (double) left.read_end / left.read_length < 1 - READ_THRESHOLD) return false;

	return true;
}

bool check_overlaps_end_end(overlap left, overlap right) {
	if (left.contig_id == right.contig_id) return false;
	if (left.strand == right.strand) return false;

	if (!(left.read_start < right.read_start && left.read_end < right.read_end || 
		right.read_start < left.read_start && right.read_end < left.read_end)) return false;
	
	if (left.contig_length - left.contig_end + right.contig_length - right.contig_end > left.read_length) return false;
	
	if ((double) left.contig_end / left.contig_length < 1 - CONTIG_THRESHOLD) return false;
	if ((double) right.contig_end / right.contig_length < 1 - CONTIG_THRESHOLD) return false;

	if (!((double) left.read_start / left.read_length <= READ_THRESHOLD && 
		(double) right.read_end / right.read_length > 1 - READ_THRESHOLD ||
		(double) right.read_start / right.read_length <= READ_THRESHOLD && 
		(double) left.read_end / left.read_length > 1 - READ_THRESHOLD)) return false;
	
	if ((double) (left.match_length + right.match_length) / left.read_length < MATCH_THRESHOLD) return false;

	return true;
}

bool check_overlaps_start_start(overlap left, overlap right) {
	if (left.contig_id == right.contig_id) return false;
	if (left.strand == right.strand) return false;

	if (!(left.read_start < right.read_start && left.read_end < right.read_end || 
		right.read_start < left.read_start && right.read_end < left.read_end)) return false;
	
	if (left.contig_start + right.contig_start > left.read_length) return false;
	
	if ((double) left.contig_start / left.contig_length > CONTIG_THRESHOLD) return false;
	if ((double) right.contig_start / right.contig_length > CONTIG_THRESHOLD) return false;

	if (!((double) left.read_start / left.read_length <= READ_THRESHOLD &&
		(double) right.read_end / right.read_length > 1 - READ_THRESHOLD ||
		(double) right.read_start / right.read_length <= READ_THRESHOLD &&
		(double) left.read_end / left.read_length > 1 - READ_THRESHOLD)) return false;

	return true;
}

void read_overlaps(string filePath) {
	cout << "Reading overlaps..." << endl;

	ifstream file(filePath);
	string line;

	while (getline(file, line)) {
		istringstream ss(line);
		overlap o;

		ss >> o.read_id >> o.read_length >> o.read_start >> o.read_end;
		ss >> o.strand;
		ss >> o.contig_id >> o.contig_length >> o.contig_start >> o.contig_end;
		ss >> o.num_matches >> o.match_length;

		overlaps.emplace_back(o);
	}
}

void read_contigs(string filePath) {
	cout << "Reading contigs..." << endl;

	ifstream file(filePath);
	string name, data;
	int ind = 0;

	while (getline(file, name)) {
		name.erase(0, 1);
		getline(file, data);

		Contig* c = new Contig(name, data);
		Scaffold* s = new Scaffold(c);

		contigs.emplace_back(c);
		scaffolds.insert(s);

		id_to_contig[name] = c;
		id_to_scaffold[name] = s;
	}
}

void rate_overlaps() {
	cout << "Checking overlaps..." << endl;

	int group_start = 0;
	
	for (int i = 1; i <= overlaps.size(); ++i) {
		if (i < overlaps.size() && overlaps[i].read_id == overlaps[i-1].read_id)
			continue;
		
		for (int j = group_start; j < i; ++j) {
			for (int k = j + 1; k < i; ++k) {
				if (check_overlaps_end_start(overlaps[j], overlaps[k])) {
					votes_end_start[std::make_pair(overlaps[j].contig_id, overlaps[k].contig_id)] += 1;
				}
				else if (check_overlaps_end_start(overlaps[k], overlaps[j])) {
					votes_end_start[std::make_pair(overlaps[k].contig_id, overlaps[j].contig_id)] += 1;
				}
				else if (check_overlaps_end_end(overlaps[j], overlaps[k])) {
					votes_end_end[std::make_pair(overlaps[j].contig_id, overlaps[k].contig_id)] += 1;
				}
				else if (check_overlaps_start_start(overlaps[j], overlaps[k])) {
					votes_start_start[std::make_pair(overlaps[j].contig_id, overlaps[k].contig_id)] += 1;
				}
			}
		}

		group_start = i;
	}

	for (auto overlap : votes_end_start)
		joint_overlaps.emplace_back(joint_overlap(overlap.first.first, overlap.first.second, 0));

	for (auto overlap : votes_end_end)
		joint_overlaps.emplace_back(joint_overlap(overlap.first.first, overlap.first.second, 1));

	for (auto overlap : votes_start_start)
		joint_overlaps.emplace_back(joint_overlap(overlap.first.first, overlap.first.second, 2));

	sort(joint_overlaps.begin(), joint_overlaps.end(), cmp);
}

void scaffold() {
	cout << "Sorting contigs..." << endl;

	for (auto o : joint_overlaps) {
		Contig* left_contig = id_to_contig[o.left];
		Contig* right_contig = id_to_contig[o.right];

		Scaffold* left_scaffold = id_to_scaffold[left_contig->get_scaffold_id()];
		Scaffold* right_scaffold = id_to_scaffold[right_contig->get_scaffold_id()];

		if (left_scaffold == right_scaffold) continue;

		if (left_contig->is_reversed())
			left_scaffold->reverse();
		
		if (right_contig->is_reversed())
			right_scaffold->reverse();

		if (o.type == 0) {
			if (!left_contig->get_end() || !right_contig->get_start()) continue;
			left_contig->set_end(false);
			right_contig->set_start(false);
		} else if (o.type == 1) {
			if (!left_contig->get_end() || !right_contig->get_end()) continue;
			left_contig->set_end(false);
			right_contig->set_end(false);
			right_scaffold->reverse();
		} else if (o.type == 2) {
			if (!left_contig->get_start() || !right_contig->get_start()) continue;
			left_contig->set_start(false);
			right_contig->set_start(false);
			left_scaffold->reverse();
		}		
		
		left_scaffold->merge(right_scaffold);
		scaffolds.erase(right_scaffold);

		if (scaffolds.size() == 1) break;
	}
}

int main(int argc, char **argv) {
	if (argc != 4) {
		cout << "ERROR: 3 arguments expected (<contigs> <overlaps> <output>)!" << endl;
		return 0;
	}

	read_contigs(string(argv[1]));
	read_overlaps(string(argv[2]));
	rate_overlaps();
	scaffold();

	ofstream output;
	output.open(string(argv[3]));

	cout << "Output to " << argv[3] << "..." << endl;

	for (auto s : scaffolds) {
		output << '>' << s->get_id() << endl;
		output << s->get_merged() << endl;
	}

	output.close();

	return 0;
}