/***************************************************************************
 *   Copyright (C) 2006 by Pier Federico Gherardini   *
 *   pier.federico.gherardini@uniroma2.it  
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

#include "Job.h"
#include "Match.h"
#include "PdbParser.h"
#include "Comparison.h"
#include <cstdio>
#include <iostream>
#include <map>
#include <cstdlib>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#endif

using std::list;
using std::cout;
using std::vector;
using std::string;
using std::pair;


class Match_rmsd_lt
{
	public:
		bool operator()(const Match *a, const Match *b)
		{
			return a->getRmsd() < b->getRmsd();
		}
};


void Job::Compare(ChainWrapper *c1, ChainWrapper *c2)
{
	if(opt.verbose)
	{
		cout << "comparing probe " << c1->getParent_chain().getFull_id(); 
		cout << " with target " << c2->getParent_chain().getFull_id() << "\n";
	}
	list<Match *> match_list;
	
	Comparison comparison(opt);
	comparison.Start(c1, c2, &match_list);
	list<Match *>::iterator i;
	static unsigned long match_id = 1; //Warning static variable here!
	if(opt.single_match)
	{
		if(match_list.size())
		{
			match_list.sort(Match_rmsd_lt());
			(match_list.front())->Write(match_id++);
			for(i = match_list.begin(); i != match_list.end(); i++)
				delete *i;
		}
	}
	else
		for(i = match_list.begin(); i != match_list.end(); i++)
		{
			(*i)->Write(match_id++);
			delete *i;
		}
}

void Job::Run(void)
{
	
	cout << "Running Job...\n";
	list<ChainWrapper *>::iterator c1, c2;
	
	
	list<ChainWrapper *> P, T;
	cout << "Creating ChainWrappers... ";
	if(opt.verbose)
		cout << "\n";
	for(list<Chain *>::iterator ch = probes.begin(); ch != probes.end(); ch++)
		P.push_back(new ChainWrapper(**ch, opt));
	for(list<Chain *>::iterator ch = targets.begin(); ch != targets.end(); ch++)
		T.push_back(new ChainWrapper(**ch, opt));
	
	cout << "Done\n";
	
	if(P.size() > 0 && T.size() == 0)
	{
		cout << "Only probes were provided, running all against all...\n";
		for(c1 = P.begin(); c1 != P.end(); c1++)
		{
			for(c2 = c1, c2++; c2 != P.end(); c2++)
				Compare(*c1, *c2);
			cout << "Done with probe " << ((*c1)->getParent_chain()).getFull_id() << "\n";
		}
	}
	else
	{
		for(c1 = P.begin(); c1 != P.end(); c1++)
		{
			for(c2 = T.begin(); c2 != T.end(); c2++)
				Compare(*c1, *c2);
			cout << "Done with probe " << ((*c1)->getParent_chain()).getFull_id() << "\n";
		}
	}
			
	cout << "Done\n";
	
	for(c1 = P.begin(); c1 != P.end(); c1++)
		delete *c1;
	for(c2 = T.begin(); c2 != T.end(); c2++)
		delete *c2;
	/* 
	else if(_mode == superpose_only)
		for(c1 = probes.begin(); c1 != probes.end(); c1++)
		{
			for(c2 = targets.begin(); c2 != targets.end(); c2++)
				Do_superposition_only(*c1, *c2, opt);
			cout << "Done with probe " << (*c1)->getId() << "\n";
		}
	*/
}


Chain *Job::Select_zone_amino(const Chain &ch, const string &zone_string)
{
	vector<string> fields;
	split(&fields, zone_string, 'z');
	string resnum = fields[0];
	float thresh = atof(fields[1].c_str());
	int idx = ch.getIndex(resnum);
	if(idx == -1)
	{
		cout << "Error! Can't find residue " << resnum
			<< " in chain " << ch.getFull_id() << "\n";
		exit(1);
	}
	vector<bool> sel(ch.size(), false);
	sel[idx] = true;
	
	if(opt.verbose == true)
		cout << "Selecting zone around " << resnum << " with threshold "
			<< thresh << "\n";
	
	for(unsigned int i = 0; i < ch.size(); i++)
	{
		if(i != (unsigned int)idx)
			if(ch[idx].Is_neighbour(ch[i], thresh))
			{
				if(opt.verbose)
					cout << "\t---selecting " << ch[i].getPdb_resnum() << "\n";
				sel[i] = true;
			}
	}
	return new Chain(ch, sel);
}

			
	

Chain *Job::Slice_chain(const Chain &ch, const string &range_string)
{
	vector<bool> sel(ch.size());
	list<pair<int, int> > ranges;
	process_range(range_string, &ranges);
	for(unsigned int i = 0; i < ch.size(); i++)
	{
		if(in_range(ranges, ch[i].getPdb_resnum()))
			sel[i] = true;
		else
			sel[i] = false;
	}
	return new Chain(ch, sel);
}

#ifdef _WIN32

void Job::List_directory_recursively(const string &dirname, vector<string> *files)
{
	WIN32_FIND_DATA file_info;
	HANDLE hfile;
	string pattern = dirname + "/*.*";
	hfile = ::FindFirstFile(pattern.c_str(), &file_info);
	
	if(hfile != INVALID_HANDLE_VALUE)
	{
		do
		{
			if(file_info.cFileName[0] != '.')
			{
				string filename = dirname + string("/") + file_info.cFileName; 
				if(file_info.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
					List_directory_recursively(filename, files);
				else
					files->push_back(filename);
				
			}
		} while(::FindNextFile(hfile, &file_info) == TRUE);
		
		::FindClose(hfile);
	}
}
	
	
	
#else
void Job::List_directory_recursively(const string &dirname, vector<string> *files)
{
	chdir(dirname.c_str());
	char *cwd = new char[2048];
	if(!getcwd(cwd, (2048 - 1) * sizeof(char)))
	{
		string msg = "Error opening subdir " + dirname;
		perror(msg.c_str());
		exit(1);
	}
		
	dirent entity;
	dirent *p = &entity;
	DIR *dir = opendir("./");
	if(!dir)
	{
		cout << "Error opening directory " << dirname << "\n";
		exit(1);
	}
	while(!readdir_r(dir, &entity, &p))
	{
		if(!p)
			break;
		else
		{
			if(p->d_name != string("..") && p->d_name != string("."))
			{
				struct stat s;
				stat(p->d_name, &s);
				if(S_ISDIR(s.st_mode))
					List_directory_recursively(p->d_name, files);
				else
					files->push_back(cwd + string("/") + p->d_name);
			}
		}
	}
	chdir("..");
	closedir(dir);
	delete [] cwd;
}

#endif //_WIN32
	
	
void Job::Load_chains(list<Chain *> *l, const vector<string> &v)
{
	for(unsigned int i = 0; i < v.size(); i++)
	{
		string::size_type first, last;
		first = v[i].find_first_of("(");
		last = v[i].find_last_of(")");
		
		if(!(first != string::npos && last != string::npos))
		{
			cout << "Missing chain specification from " << v[i] << ", aborting...\n";
			exit(1);
		}
		string s = v[i].substr(first + 1, last - (first + 1));
		vector<string> fields;
		split(&fields, s, ':');
		string ch_id = fields[0];
		vector<string> files_to_load;
		
		if(v[i][0] == '[')
		{
			string::size_type k = v[i].find_first_of(']');
			string dirname = v[i].substr(1, k - 1);
			List_directory_recursively(dirname, &files_to_load);
		}
		else
			files_to_load.push_back(opt.pdb_dir + "/" + v[i].substr(0, first));
		
		for(unsigned int j = 0; j < files_to_load.size(); j++)
		{
			const string &file_name = files_to_load[j];
			PdbParser parser(opt, file_name);
			vector<string> chains_to_load;
			if(ch_id == "*")
				parser.Get_chain_list(&chains_to_load);
			else
				chains_to_load.push_back(ch_id);
			for(unsigned int a = 0; a < chains_to_load.size(); a++)
			{
				Chain *ch = parser.Load_chain(chains_to_load[a]);
				if(!ch->size())
				{
					cout << "Error chain " << ch->getFull_id() << " is empty\n";
					delete ch;
				}
				else
				{
#ifndef NDEBUG
					ch->Write(ch->getFull_id() + ".conv.pdb");
#endif
					if(fields.size() == 2)
					{
						Chain *temp = 0;
						if(fields[1].find_first_of("z") != string::npos)
							temp = Select_zone_amino(*ch, fields[1]);
						else
							temp = Slice_chain(*ch, fields[1]);
						delete ch;
						ch = temp;
					}
					l->push_back(ch);
				}
			}	
		}
	}
}




Job::Job(const vector<string> &_probes, const vector<string> &_targets, const options &_opt) : opt(_opt)
{
	Load_chains(&probes, _probes);
	Load_chains(&targets, _targets);	
}
			
Job::~Job(void)
{
	list<Chain *>::iterator i;
	for(i = probes.begin(); i != probes.end(); i++)
		if(*i)
			delete *i;
	
	for(i = targets.begin(); i != targets.end(); i++)
		if(*i)
			delete *i;
}



	
/*
void Job::Do_superposition_only(Chain *c1, Chain *c2, const options &opt)
{
	static unsigned long match_id = 1; //Warning static variable here!
	assert(c1->GetSize() == c2->GetSize());
	short *amino1, *amino2;
	amino1 = new short[c1->GetSize()];
	amino2 = new short[c2->GetSize()];
	for(unsigned int i = 0; i < c1->GetSize(); i++)
		amino1[i] = amino2[i] = i;
	
	Match *m = new Match(c1, c2);
	//prof and rms are useless, set them to 0
	m->Set_values(c1->GetSize(), 0, 0, amino1, amino2);
	m->Write(opt, match_id++);
	delete m;
	delete [] amino1;
	delete [] amino2;
}
*/
	
	
	


