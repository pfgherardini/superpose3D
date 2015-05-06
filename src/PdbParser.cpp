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


#include "PdbParser.h"
#include "ResidueDefinition.h"
#include "PdbResidueDictionary.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "common.h"

using std::string;
using std::map;
using std::ifstream;
using std::cout;
using std::vector;


void PdbParser::Handle_amino_end(Amino **a, Chain *c)
{
	if(*a)
	{
		if((*a)->getLabel() == "GLY")
			(*a)->Reconstruct_cb();
		Amino *amino = Convert_pdb_amino(**a);
		if(amino)
			c->Add_residue(amino);
		delete *a;
		*a = 0;
	}
}

inline string PdbParser::Get_chain(const string &line)
{
	if(line[21] == ' ')
		return "_";
	else
		return line.substr(21, 1);
}


void PdbParser::Read_lines(const string &input_file)
{
	lines.clear();
	ifstream in(input_file.c_str());
	string temp;
	
	while(getline(in, temp) && !in.eof())
	{
		//Only load relevant lines, avoid miscellanous crap see e.g. 2gqu
		if(temp.substr(0, 5) == "MODEL" ||
		   temp.substr(0, 4) == "ATOM" ||
		   temp.substr(0, 3) == "TER" ||
		   temp.substr(0, 6) == "HETATM" ||
		   temp.substr(0, 6) == "ENDMDL" ||
		   temp.substr(0, 3) == "END")
			lines.push_back(temp);
	}
	in.close();
}

unsigned int PdbParser::Find_chain_start(const string &chain)
{
	for(unsigned int i = 0; i < lines.size(); i++)
	{
		const string &s = lines[i];
		if(s.substr(0, 4) == "ATOM" && Get_chain(s) == chain)
			return i;
	}
	cout << "Chain " << chain << " not found, aborting...\n";
	exit(1);
}




inline bool PdbParser::Chain_finished(const string &cur_chain, const string &line)
{
	if(line.substr(0, 3) == "TER" ||
		line.substr(0, 3) == "END" ||
		line.substr(12, 4) == "OXT" ||
		line.substr(0, 6) == "ENDMDL" ||
		(Get_chain(line) != cur_chain))
		return true;
	else
		return false;
}


void PdbParser::Calculate_chain_offsets(const string &chain, vector<float> *offsets, unsigned int *chain_start)
{
	offsets->resize(3);
	unsigned int i = Find_chain_start(chain);
	*chain_start = i;
	float max[3], min[3];
	const int x = 0;
	const int y = 1;
	const int z = 2;
	
	const string &s = lines[i++]; //Init values with first line and discard it
	max[x] = min[x] = (float)atof(trim_copy(s.substr(30, 8)).c_str());
	max[y] = min[y] = (float)atof(trim_copy(s.substr(38, 8)).c_str());
	max[z] = min[z] = (float)atof(trim_copy(s.substr(46, 8)).c_str());
	
	for(; i < lines.size(); i++)
	{
		const string &line = lines[i];
		if(Chain_finished(chain, line))
			break;
		float cur[3];
		cur[x] = (float)atof(trim_copy(line.substr(30, 8)).c_str());
		cur[y] = (float)atof(trim_copy(line.substr(38, 8)).c_str());
		cur[z] = (float)atof(trim_copy(line.substr(46, 8)).c_str());
		
		for(unsigned int j = 0; j < 3; j++)
		{
			if(cur[j] > max[j]) 
				max[j] = cur[j];
			if(cur[j] < min[j])
				min[j] = cur[j];
		}
	}
	for(unsigned int j = 0; j < 3; j++)
		(*offsets)[j] = (max[j] + min[j]) / 2;
}
		
Atom *PdbParser::Create_atom(const string &s, const vector<float> &offsets)
{
	Atom *ret = new Atom;
	ret->x = (int)((atof(s.substr(30, 8).c_str()) - offsets[0]) * 100);
	ret->y = (int)((atof(s.substr(38, 8).c_str()) - offsets[1]) * 100);
	ret->z = (int)((atof(s.substr(46, 8).c_str()) - offsets[2]) * 100);	
	ret->label = trim_copy(s.substr(12, 4));
	return ret;
}


bool PdbParser::Check_amino_integrity(const vector<vector<const Atom *> > &m)
{
	bool retval = false;
	for(unsigned int i = 0; i < m.size(); i++)
		if(m[i].size())
			retval = true;
		else if(opt.mode == options::mode_res) //If we are in residue mode residues need to be 100% complete
			retval = false;
	return retval;
}
	

inline int ROUND(float a)
{
	return int(a > 0.0 ? a + 0.5 : a - 0.5);
}
		
//ATTENZIONE!!!! Qui volendo forse converrebbe fare che la label dell' amminoacido
//che viene creato sia quella definita dall'utente

Amino *PdbParser::Convert_pdb_amino(const Amino &a)
{
	const ResidueDefinition &def = def_center.Get_definition(a.getLabel());
	Amino *res = new Amino(def.getNum_atoms());
	res->setPdb_resnum(a.getPdb_resnum());
	res->setLabel(a.getLabel());
	
	vector<vector<const Atom *> > m(res->getNum_atoms());
	
	for(unsigned int i = 0; i < a.getNum_atoms(); i++)
	{
		short idx = def.Atom2index((a[i])->label);
		if(idx >= 0)
		{
			assert((unsigned int)idx < m.size());
			m[idx].push_back(a[i]);
		}
	}
	
	if(Check_amino_integrity(m))
	{
		for(unsigned int i = 0; i < m.size(); i++)
		{
			if(m[i].size())
			{
				(*res)[i] = new Atom;
				for(unsigned int j = 0; j < m[i].size(); j++)
				{
					((*res)[i])->x += (m[i][j])->x;
					((*res)[i])->y += (m[i][j])->y;
					((*res)[i])->z += (m[i][j])->z;
				}
				//ATTENZIONE!!!! BISOGNA RISOLVERE LA QUESTIONE DELL'ARROTONDAMENTO!!!
				((*res)[i])->x = ROUND(((*res)[i])->x / (float)(m[i].size()));
				((*res)[i])->y = ROUND(((*res)[i])->y / (float)(m[i].size()));
				((*res)[i])->z = ROUND(((*res)[i])->z / (float)(m[i].size()));
			}
		}
	}
	else
	{
		cout << "Residue " << res->getPdb_resnum() <<
			" is incomplete according to the current definition and will be removed\n";
		delete res;
		res = 0;
	}
	return res;
}

void PdbParser::Get_chain_list(vector<string> *v)
{
	if(!lines.size())
		Read_lines(file_name);
	string temp;
	string prev_chain;
	v->clear();
	
	for(unsigned int i = 0; i < lines.size(); i++)
	{
		const string &temp = lines[i];
		if(temp.substr(0, 4) == "ATOM")
		{
			string id = Get_chain(temp);
			if(id != prev_chain)
				v->push_back(id);
			prev_chain = id;
		}
	}
}
		
	

Chain *PdbParser::Load_chain(const string &chain_id)
{
	cout << "Loading " + file_name + chain_id + "\n";
	if(!lines.size())
		Read_lines(file_name);
	string temp;
	Amino *cur_amino = 0;
	Chain *cur_chain = new Chain(base_name(file_name, false), chain_id);
	vector<float> offsets(3);
	unsigned int chain_start = 0;
	Calculate_chain_offsets(chain_id, &offsets, &chain_start);
	cur_chain->setOffsets(offsets);

	for(unsigned int i = chain_start; i < lines.size(); i++)
	{
		const string &temp = lines[i];
		if(Chain_finished(cur_chain->getChain_id(), temp))
		{ 
			Handle_amino_end(&cur_amino, cur_chain);
			break;
		}
		if(!PdbResidueDictionary::Is_amino(temp))
			continue;
		//We paste together resnum and insertion code
		string resnum = trim_copy(temp.substr(22, 5));
		if(!cur_amino || cur_amino->getPdb_resnum() != resnum)
		{
			if(cur_amino)
				Handle_amino_end(&cur_amino, cur_chain);
			cur_amino = new Amino;
			cur_amino->Load_from_pdb_line(temp);
		}
		Atom *a = Create_atom(temp, offsets);
		if(temp[16] != cur_amino->getAlt_loc())
		{
			if(cur_amino->getAlt_loc() == ' ')
				cur_amino->setAlt_loc(temp[16]);
			else
			{
				delete a;
				continue;
			}
		}
		if(!cur_amino->Add_atom(a)) delete a;
	}
	if(cur_amino)
		Handle_amino_end(&cur_amino, cur_chain);
	return cur_chain;
}


