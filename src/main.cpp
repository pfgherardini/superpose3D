/***************************************************************************
 *   Copyright (C) 2007-2009 by Pier Federico Gherardini   *
 *   pier.federico.gherardini@uniroma2.it   *
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


#include "DefinitionCenter.h"
#include "SimilarityMatrix.h"
#include "PdbParser.h"
#include "Chain.h"
#include "Amino.h"
#include "ChainWrapper.h"
#include "common.h"
#include "Comparison.h"
#include "Match.h"
#include "Job.h"
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

void Parse_parameters_file(const string &file_name, options *opt);
void Parse_command_line(int argc, char *argv[], options *opt);
void Parse_input_file(const string &file_name, options *opt, vector<string> *probes, vector<string> *targets);

int main(int argc, char *argv[])
{
	setvbuf(stdout, NULL, _IONBF,0);
	options opt;
	vector<string> probes, targets;
	Parse_command_line(argc, argv, &opt);
	if(opt.param_file_name == "")
	{
		cout << "Error! No parameters file has been specified, aborting...\n";
		exit(1);
	}
	if(opt.input_file_name == "")
	{
		cout << "Error! No input file has been specified, aborting...\n";
		exit(1);
	}
	
	Parse_parameters_file(opt.param_file_name, &opt);
		
	
	Parse_input_file(opt.input_file_name, &opt, &probes, &targets);
	
	if(!(opt.output_file))
	{
		string temp = opt.input_file_name;
		temp = base_name(temp) + "_out.txt";
		opt.output_file = opt.trans_file = fopen(temp.c_str(), "w");
	}
	
#ifndef NDEBUG
	opt.sim_mat->Print();
	opt.sim_mat->Print_index_map();
#endif	
	Job *job = new Job(probes, targets, opt);
	job->Run();
	
	fclose(opt.output_file);
	return 0;
}


void Parse_input_file(const string &file_name, options *opt, vector<string> *probes, vector<string> *targets)
{
	ifstream in(file_name.c_str());
	if(in.rdstate() & ifstream::failbit)
	{
		cout << "Error! Can't open " << file_name << ", aborting...\n";
		exit(1);
	}
		
	string temp;
	
	
	while(getline(in, temp) && !in.eof())
	{
		if(temp.size())
		{
			vector<string> fields;
			split(&fields, temp, ' ');
			assert(fields.size() == 2);
			if(fields[0] == "probe")
				split(probes, fields[1], ';');
			else if(fields[0] == "target")
				split(targets, fields[1], ';');
			else if(fields[0] == "pdb_dir")
				opt->pdb_dir = fields[1];
			else if(fields[0] == "output_file")
			{
				if(opt->output_file)
					cout << "Warning! an output file has already been " <<
					 "specified on the command line, I'll use that instead\n";
				else
					opt->output_file = opt->trans_file = fopen(fields[1].c_str(), "w");
			}
		}
	}
	in.close();
}
	

void Parse_command_line(int argc, char *argv[], options *opt)
{
	char c;
	for(int i = 0; i < argc; i++)
		if(argv[i][0] == '-')
			switch(c = argv[i][1])
			{
				case 'i':
					opt->input_file_name = argv[i + 1]; break;
				case 'v':
					opt->verbose = true; break;
				case 'o':
					opt->trans_file = opt->output_file = fopen(argv[i + 1], "w"); break;
				case 'p':
					opt->param_file_name = argv[i + 1];
					break;
				default:
					cout << "Unrecognized option!" << c << endl; break;
			}
}
			
			
			


void Parse_parameters_file(const string &file_name, options *opt)
{
	ifstream in(file_name.c_str());
	if(in.rdstate() & ifstream::failbit)
	{
		cout << "Error! Can't open " << file_name << ", aborting...\n";
		exit(1);
	}
	
	string temp;
	vector<string> res_def, sim_matrix;
	
	
	while(getline(in, temp) && !in.eof())
	{
		if(temp.size() && temp[0] != '#')
		{
			vector<string> fields;
			split(&fields, temp, ' ');
			assert(fields.size() == 2);
			if(temp.substr(0, 3) == "def")
				res_def.push_back(fields[1]);
			else if(temp.substr(0, 5) == "equiv")
				sim_matrix.push_back(fields[1]);
			else if(temp.substr(0, 4) == "mode")
			{
				if(fields[1] == "res")
					opt->mode = options::mode_res;
				else if(fields[1] == "atom")
					opt->mode = options::mode_atom;
				else if(fields[1] == "res_multiple")
					opt->mode = options::mode_res_multiple;
				else
				{
					cout << "Unknown mode " << fields[1] << " aborting...\n";
					exit(1);
				}
			}
			else if(fields[0] == "rmsd_threshold")
				opt->rmsd_thresh = (float)atof(fields[1].c_str());
			else if(fields[0] == "neighbour_threshold")
				opt->neighbour_threshold = (float)atof(fields[1].c_str());
			else if(fields[0] == "score_max")
				opt->score_max = atoi(fields[1].c_str());
			else if(fields[0] == "score_min")
				opt->score_min = atoi(fields[1].c_str());
			else if(fields[0] == "max_longest_matches")
				opt->max_longest_matches = atoi(fields[1].c_str());
			else if(fields[0] == "pdb_res_dict")
				opt->pdb_res_dict = new PdbResidueDictionary(fields[1]);
			else
			{
				cout << "Error! Unrecognized option key " << fields[0] << "\n";
				exit(1);
			}
		}
	}
	opt->def_center = new DefinitionCenter(res_def, *opt);
	opt->sim_mat = new SimilarityMatrix(sim_matrix, *opt);

	in.close();
	
}
