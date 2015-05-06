/***************************************************************************
 *   Copyright (C) 2006 by Federico Gherardini   *
 *   pier.federico.gherardini@uniroma2.it  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#include <list>

#include "common.h"
#include <cassert>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cctype>
#include <cstdlib>


using namespace std;

struct Match
{
	vector<string> probe_res, target_res;
	string probe_chain, target_chain;
	vector<float> probe_trasl, target_trasl;
	vector<vector<float> > rot_mat;
	string match_num;
	
	Match(const vector<string> &fields);
	void Get_residue_numbers(const string &s, string *p, string *t);
	string Get_residue_string(const vector<string> &amino, const string &chain, const string &model_num);
	void Write_chimera_file(const string &probe_file, const string &target_file);
};

class PdbWriter
{
	private:
		string Find_file(const vector<string> &dir_list, const string &file_name);
	public:
		void Parse(const string &file_name, const vector<string> &dir_list, const string &out_dir);

		static const int match_len_idx = 6;
		static const int pairing_start_idx = 7;
		static const int probe_idx = 1;
		static const int target_idx = 2;
		static const int match_num_idx = 0;
};


		
template<class T> std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b)
{
	assert(a.size() == b.size());
	std::vector<T> res(a.size());
	res = a;
	for(unsigned long i = 0; i < b.size(); i++)
		res[i] += b[i];
	return res;
}
	
template<class T> std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b)
{
	assert(a.size() == b.size());
	std::vector<T> res(a.size());
	res = a;
	for(unsigned long i = 0; i < b.size(); i++)
		res[i] -= b[i];
	
	return res;
}


template<class T, class U> std::vector<T>operator/(const std::vector<T> &v, const U &m)
{
	std::vector<T> res(v.size());
	for(unsigned long i = 0; i < v.size(); i++)
		res[i] = v[i] / m;
	return res;
}
	

void rotate(vector<float> *v, const vector<vector<float> > &m)
{
	vector<float> res(3);
	for(int i = 0; i < 3; i++)
	{
		float temp = 0;
		for(int j = 0; j < 3; j++)
			temp += m[i][j] * (*v)[j];
		res[i] = temp;
	}
	for(int i = 0; i < 3; i++)
		(*v)[i] = res[i];
}


	
void rotate_pdb_file(const string &in_name, const string &out_name, 
		     const vector<float> &trasl, const vector<vector<float> > &m, bool to_rotate)
{
	string temp;
	ifstream in(in_name.c_str());
	FILE *out = fopen(out_name.c_str(), "w");
	while(getline(in, temp) && !in.eof())
	{
		if(temp.substr(0, 4) == "ATOM" ||
				 temp.substr(0, 6) == "HETATM")
		{
			string coord = temp.substr(30, 24);
			coord = trim_copy(coord);
			float x, y, z;
			x = y = z = 0.0f;
			
			if(sscanf(coord.c_str(), "%f%f%f", &x, &y, &z) != 3)
			{
				cout << "Error while parsing x y z coordinates, aborting...\n";
				exit(1);
			}
			vector<float> v(3);
			v[0] = x; v[1] = y; v[2] = z;
			v = v + trasl;
			if(to_rotate)
				rotate(&v, m);
			fprintf(out, temp.substr(0, 30).c_str());
			for(int i = 0; i < 3; i++)
				fprintf(out, "%8.3f", v[i]);
			fprintf(out, temp.substr(54, temp.size() - 54).c_str());
			fprintf(out, "\n");
				
		}
		else
			fprintf(out, "%s\n", temp.c_str());
	}
	in.close();
	fclose(out);
}






void Match::Get_residue_numbers(const string &s, string *p, string *t)
{
	unsigned int first, last;
	first = s.find_first_of("(") + 1;
	last = s.find_first_of(")");
	*p = s.substr(first, last - first);
	first = s.find_first_of("(", last + 1) + 1;
	last = s.find_first_of(")", last + 1);
	*t = s.substr(first, last - first);
}


Match::Match(const vector<string> &fields) : probe_trasl(3), target_trasl(3), rot_mat(3)
{
	for(unsigned int k = 0; k < 3; k++)
		rot_mat[k].resize(3);
	
	match_num = fields[PdbWriter::match_num_idx];
	const string *s;
	s = &(fields[PdbWriter::probe_idx]);
	probe_chain = s->substr(s->size() - 1);
	s = &(fields[PdbWriter::target_idx]);
	target_chain = s->substr(s->size() - 1);
	
	
	unsigned int match_len = atoi(fields[PdbWriter::match_len_idx].c_str());
	unsigned int start_trans = PdbWriter::pairing_start_idx + match_len;
	unsigned int i, j;
	
	
	for(i = PdbWriter::pairing_start_idx, j = 0; i < start_trans; i++, j++)
	{
		string p, t;
		Get_residue_numbers(fields[i], &p, &t);
		probe_res.push_back(p);
		target_res.push_back(t);
	}
	
	for(i = start_trans, j = 0; i < start_trans + 3; i++, j++)
		probe_trasl[j] = (float)atof(fields[i].c_str());
	start_trans = i;
	for(i = start_trans, j = 0; i < start_trans + 3; i++, j++)
		target_trasl[j] = (float)atof(fields[i].c_str());
	start_trans = i;
	unsigned int r, c;
	for(i = start_trans, r = 0, c = 0; i < start_trans + 9; i++)
	{
		c = (i - start_trans) % 3;
		r = (i - start_trans) / 3;
		rot_mat[r][c] = (float)atof(fields[i].c_str());
	}
}

string Match::Get_residue_string(const vector<string> &amino, const string &chain, const string &model_num)
{
	string ret;
	ret = "#" + model_num;
	for(unsigned int i = 0; i < amino.size(); i++)
		ret += ":" + amino[i] + "." + chain;
	return ret;
}


void Match::Write_chimera_file(const string &probe_file, const string &target_file)
{
	string out_file = "chim_" + match_num + ".cmd";
	ofstream out(out_file.c_str());
	out << "close all\nopen " << probe_file << "\nopen " << target_file << "\n";
	out << "~show\n";
	string probe_res_s = Get_residue_string(probe_res, probe_chain, "0");
	string target_res_s = Get_residue_string(target_res, target_chain, "1");
	
	out << "display " << probe_res_s << "\n";
	out << "display " << target_res_s << "\n";
	out << "repr stick " << probe_res_s << "\n";
	out << "repr stick " << target_res_s << "\n";
	out << "focus " << probe_res_s << "\n";
	out.close();
}
	



string PdbWriter::Find_file(const vector<string> &dir_list, const string &file_name)
{
	string retval = "";
	for(unsigned int i = 0; i < dir_list.size(); i++)
	{
		string name = dir_list[i] + "/" + file_name;
		ifstream in(name.c_str());
		if(!(in.rdstate() & ifstream::failbit))
		{
			retval = name;
			in.close();
			break;
		}
		else
			in.close();
	}
	return retval;
}
	
	



void PdbWriter::Parse(const string &file_name, const vector<string> &dir_list, const string &out_dir)
{
	ifstream in(file_name.c_str());
	string temp;
		
	while(getline(in, temp) && !in.eof())
	{
		vector<string> fields;
		split(&fields, temp, '\t');
		Match m(fields);
		//We consider the last letter to be the chain name
		string probe_name = fields[probe_idx].substr(0, fields[probe_idx].size() - 1);
		string target_name = fields[target_idx].substr(0, fields[target_idx].size() - 1);
		
		string probe_in = Find_file(dir_list, probe_name);
		string target_in = Find_file(dir_list, target_name);
		if(probe_in == "")
		{
			cout << "Error: can't find file " << probe_name <<
			" in any of directories you provided\n";
			exit(1);
		}
		if(target_in == "")
		{
			cout << "Error: can't find file " << target_name <<
			" in any of directories you provided\n";
			exit(1);
		}
		
		string probe_out = out_dir + "/" + base_name(probe_in) + string("_rot") + fields[match_num_idx] + ".ent";
		string target_out = out_dir + "/" + base_name(target_in) + string("_rot") + fields[match_num_idx] + ".ent";
		
		rotate_pdb_file(probe_in, probe_out,  m.probe_trasl, m.rot_mat, false);
		rotate_pdb_file(target_in, target_out, m.target_trasl, m.rot_mat, true);
		m.Write_chimera_file(probe_out, target_out);
		
	}
	in.close();
}

void parse_options(int argc, const char *argv[], string *input_file, vector<string> *dir_list, string *out_dir)
{
    int i;
    char c;
	
    for(i = 0; i < argc; i++)
    {
		if(argv[i][0] == '-')
			switch(c = argv[i][1])
			{
				case 'i':
					*input_file = argv[i + 1]; break;
				case 'd':
					split(dir_list, argv[i + 1], ','); break;
				case 'o':
					*out_dir = argv[i + 1]; break;
				default:
					cout << "Unrecognized option!" << c << endl; break;
			}
    }
}



int main(int argc, const char *argv[])
{
	if(argv[1] == string("-h"))
	{
		cout << "superpose_pdb_files [-i input_file] [-d list of pdb file dir] [-o output_dir]\n";
		return 0;
	}
	else
	{
		string input_file, out_dir;
		out_dir = ".";
		vector<string> dir_list;
		parse_options(argc, argv, &input_file, &dir_list, &out_dir);
		PdbWriter *parser = new PdbWriter;
		parser->Parse(input_file, dir_list, out_dir);
		delete parser;
	}
	return 0;
}
	




