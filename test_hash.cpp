#include <iostream> 
#include <fstream>
#include <sstream>
#include <unordered_map> 
#include <cstring>
#include <iterator>
#include <vector>
#include <algorithm>
#include <list>
#include <regex>
#include "ksw2.h"

using namespace std; 

unordered_map<int, unordered_map<int, unordered_map<string, string>>> genome_cache;

string cache_finding (int off_q, int len_q, string sub_ref){
	if(genome_cache.find(off_q)!= genome_cache.end()){
		unordered_map<int, unordered_map<string, string>> sub_q = genome_cache[off_q];
		if(sub_q.find(len_q) != sub_q.end()){
			unordered_map<string, string> sub_r = sub_q[len_q];
			if(sub_r.find(sub_ref) != sub_r.end()){
				return sub_r[sub_ref];
			}
		}
	}
	return "";
}
void add_to_cache (int off_q, int len_q, string sub_ref, string alignment){
	if(genome_cache.find(off_q) == genome_cache.end()){
		unordered_map<string, string> res;
		res.insert(make_pair(sub_ref, alignment));
		unordered_map<int, unordered_map<string, string>> sub_r;
		sub_r.insert(make_pair(len_q, res));
		genome_cache.insert(make_pair(off_q,sub_r));
	}
	else{
		unordered_map<int, unordered_map<string, string>> sub_r = genome_cache[off_q];
		if(sub_r.find(len_q) == sub_r.end()){
			unordered_map<string, string> res;
			res.insert(make_pair(sub_ref, alignment));
			sub_r.insert(make_pair(len_q,res));
		}
		else{
			unordered_map<string, string> res = sub_r[len_q];
			if(res.find(sub_ref) == res.end()){
				res.insert(make_pair(sub_ref, alignment));
				sub_r.insert(make_pair(len_q,res));
			}
		}
	}
}
string global_aligment(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape){
	void *km = 0;
	// int8_t q = 4, e = 2;
	int w = -1, rep = 1;
	ksw_extz_t ez;
	int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
	int tl = strlen(tseq), ql = strlen(qseq);
	uint8_t *ts, *qs, c[256];
	memset(&ez, 0, sizeof(ksw_extz_t));
	memset(c, 4, 256);
	c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
	c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
	ts = (uint8_t*)malloc(tl);
	qs = (uint8_t*)malloc(ql);
	for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
	for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
	ez.score = ksw_gg(km, ql, (uint8_t*)qseq, tl, (uint8_t*)tseq, 5, mat, gapo, gape, w, &ez.m_cigar, &ez.n_cigar, &ez.cigar);
	string rs = "";
	for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
		rs = rs + to_string(ez.cigar[i]>>4) + "MID"[ez.cigar[i]&0xf];
	free(ez.cigar); free(ts); free(qs);
	return rs;
}
string extension_alignment(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape){
	int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
	int tl = strlen(tseq), ql = strlen(qseq);
	uint8_t *ts, *qs, c[256];
	ksw_extz_t ez;

	memset(&ez, 0, sizeof(ksw_extz_t));
	memset(c, 4, 256);
	c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
	c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
	ts = (uint8_t*)malloc(tl);
	qs = (uint8_t*)malloc(ql);
	for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
	for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
	ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
	// printf("%d", ez.n_cigar);
	string rs = "";
	for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
		rs = rs + to_string(ez.cigar[i]>>4) + "MID"[ez.cigar[i]&0xf];
	free(ez.cigar); free(ts); free(qs);
	return rs;
}
string alignment(string ref, string q, list<int>ref_index, list<int>q_index, list<int>match_len, int sc_mch, int sc_mis, int gapo, int gape){
	string result = "";
	int begin_match_ref= ref_index.front();
	int begin_match_q = q_index.front();
	
	if(begin_match_q == 0 && begin_match_ref != 0 ){
		result += to_string(begin_match_ref) +"D";
	}
	if(begin_match_q != 0 && begin_match_ref == 0){
		result += to_string(begin_match_q) +"I";
	}
	if(begin_match_q != 0 && begin_match_ref != 0){
		string begin_ref = ref.substr(0,begin_match_ref);
		string begin_q = q.substr(0,begin_match_q);
		string align = cache_finding(0, begin_match_q, begin_ref);

		if( align == ""){
			reverse(begin_ref.begin(), begin_ref.end());
			reverse(begin_q.begin(), begin_q.end());
			string align_rs = extension_alignment(begin_ref.c_str(), begin_q.c_str(), sc_mch, sc_mis, gapo, gape);
			reverse(align_rs.begin(), align_rs.end());
			result += align_rs;
			add_to_cache(0, begin_match_q,begin_ref, align_rs);
		}
		else{
			result += align;
		}
		
	}

	std::vector<int> ref_vector{begin(ref_index), end(ref_index)};
	std::vector<int> q_vector{begin(q_index), end(q_index)};
	std::vector<int> len_vector{begin(match_len), end(match_len)};
	for(int i = 0; i< q_index.size()-1; ++i){
		int off_q = q_vector[i] + len_vector[i];
		int len_q = q_vector[i+1]-off_q;
		string sub_q = q.substr(off_q, len_q);
		string sub_ref = ref.substr(ref_vector[i]+len_vector[i], ref_vector[i+1]- ref_vector[i]-len_vector[i]);
		result += to_string(len_vector[i]) + "M";
		string find_cache_rs = cache_finding(off_q, len_q, sub_ref);
		if(find_cache_rs != ""){
			result += find_cache_rs;
		}
		else{
			// sub_ref = sub_ref.substr(0, 15);
			// sub_q = sub_q.substr(0,10);
			string glob_align = global_aligment(sub_ref.c_str(), sub_q.c_str(), sc_mch, sc_mis, gapo, gape);
			add_to_cache(off_q, len_q, sub_ref, glob_align);
			result += glob_align;
		}

	}
	int end_match_ref = ref_index.back();
	int end_match_q = q_index.back();
	int end_match_len = match_len.back(); 
	if(end_match_q + end_match_len == q.length() && end_match_ref + end_match_len < ref.length()){
		result += to_string(end_match_len)+"M";
		result += to_string(ref.length()- end_match_ref- end_match_len) +"D";
	}
	else if(end_match_q + end_match_len < q.length() && end_match_ref + end_match_len == ref.length()){
		result += to_string(end_match_len)+"M";
		result += to_string(q.length()- end_match_q- end_match_len) +"I"; 
	}
	else if (end_match_q + end_match_len == q.length() && end_match_ref + end_match_len == ref.length()){
		result += to_string(end_match_len)+"M";
	}
	else{
		result += to_string(end_match_len)+"M";
		end_match_ref = end_match_ref + end_match_len;
		end_match_q = end_match_q + end_match_len;
		int end_len_q = q.length() - end_match_q;
		string end_sub_ref = ref.substr(end_match_ref, ref.length()-end_match_ref);
		string end_align = cache_finding(end_match_q, end_len_q, end_sub_ref);
		if(end_align == ""){
			string end_sub_q = q.substr(end_match_q, end_len_q);
			// end_sub_ref = end_sub_ref.substr(0,6);
			// end_sub_q = end_sub_q.substr(0,4);
			string end_align_rs = extension_alignment(end_sub_ref.c_str(), end_sub_q.c_str(), sc_mch, sc_mis, gapo, gape);
			result += end_align_rs;
			// cout << "this is it "<< end_sub_ref.length() <<endl;
			// cout << "and " << end_sub_q.length() << endl;
		}
		else{
			result += end_align;
		}
		result += end_align;



	}
	return result;
}

string read_file(string filename){
	ifstream file(filename);
	string line, id, genome;
	while(getline(file, line)){
		// cout << line << endl;
		if(line.empty()){
			continue;
		}
		if(line[0] == '>'){
			id = line.substr(1);
		}
		else{
			genome += line;
		}
	}
	file.close();
	return genome;
}

vector<vector<string>> read_exact_match(string filename){
	ifstream file(filename);
	std::vector<vector<string>> match_pair;
	string line, id;
	while(getline(file, line)){
		// cout << line << endl;
		if(line.empty()){
			continue;
		}
		if(line[0] == '>'){
			id = line.substr(1);
		}
		else{
			line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1");
			istringstream buf(line);
    		istream_iterator<string> beg(buf), end;
   			vector<string> tokens(beg, end);
			match_pair.push_back(tokens);
		}
	}
	file.close();
	return match_pair;
}  
int main() 
{ 	
	string query = read_file("./test/query.fasta");
	// cout << "all"<<query.length() << endl;
	string ref = read_file("./test/ref.fasta");
	// cout <<"all ref "<< ref.length()<<endl;
	// list<int> ref_index {9375, 9446, 9535, 9574, 9598, 9629, 9676, 9727, 9776};
	// list<int> q_index {47, 118, 207, 246, 270, 301, 348, 399, 448};
	// list<int> len {28, 28, 26, 23, 30, 33, 29, 48, 20};
	// cout<<alignment(ref, query, ref_index, q_index, len, 1, -2, 2, 1).c_str()<<endl;
	// string s1 = "atATGCTAcCGCat";
	// string s2 = "AT";
	// printf("%s",extension_alignment(s1.c_str(),s2.c_str(), 1, -2, 2, 1).c_str());
	list<int> ref_index ;
	list<int> q_index ;
	list<int> len ;
	vector<vector<string>> matches = read_exact_match("result.txt");
	for(auto& apair:matches){
		// cout << apair[0] <<" "<< apair[1]<< " "<<apair[2]<<endl;
		ref_index.push_back(stoi(apair[0]));
		q_index.push_back(stoi(apair[1]));
		len.push_back(stoi(apair[2]));
	}

	cout<<alignment(ref, query, ref_index, q_index, len, 1, -2, 2, 1).c_str()<<endl;
	return 0;

} 