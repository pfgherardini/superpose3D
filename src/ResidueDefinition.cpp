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


#include "ResidueDefinition.h"
#include <cassert>
#include <iostream>
#include "common.h"
#include <cstdlib>

using std::string;
using std::map;
using std::cout;
using std::vector;



short ResidueDefinition::Find_atom_with_regex(const map<string, unsigned short> &m, const string &key) const
{
	map<string, unsigned short>::const_iterator x;
	short retval = -1;
	for(x = m.begin(); x != m.end(); x++)
	{
		if((x->first)[0] == '\\')
			if(key.find(x->first.substr(1)) != string::npos)
			{
				retval = x->second;
				break;
			}
	}
	return retval;
}
		
	
	


short ResidueDefinition::Atom2index(const std::string &s) const
{
	map<string, unsigned short>::const_iterator x = atom2index.find(s);
	short retval = -1;
	if(x == atom2index.end())
	{
		//A specification for the side chain has been provided
		if(s != "N" && s != "CA" &&
		   s != "C" && s != "O" &&
		   ((x = atom2index.find("side_chain")) != atom2index.end()))
			retval = x->second;
		/*else if(label == "*")  //If this is a catch_all definition, search for partial name matches in case these have been defined
			retval = Find_atom_with_regex(atom2index, s);*/
	}
	else
		retval = x->second;
	return retval;
}




unsigned short ResidueDefinition::Pseudoatom_name2index(const std::string &s) const
{
	std::map<std::string, unsigned short>::const_iterator x = pseudoatom_name2index.find(s);
	assert(x != pseudoatom_name2index.end());
	return x->second;
/*	else
	{
		//If there is a regular expression (e.g. \N) that has not been expanded this
		//must be a catch_all definition
		assert(label == "*");
		short temp = Find_atom_with_regex(pseudoatom_name2index, s);
		assert(temp >= 0);
		retval = temp;
	}*/
}



void ResidueDefinition::Expand_wildcards_outside_avg(vector<string> *res)
{
	vector<string> v = *res;
	res->clear();
	for(unsigned int i = 0; i < v.size(); i++)
	{
		const string &s = v[i];
		if(s[0] == '\\')
		{
			vector<string> temp;
			pdb_res_dict.Expand_wildcard(this->label, s, &temp);
			vector<string>::iterator j;
			for(j = temp.begin(); j != temp.end(); j++)
				res->push_back(*j);
		}
		else
			res->push_back(s);
	}
}
	

ResidueDefinition::ResidueDefinition(const string &s, const PdbResidueDictionary &_dict) : pdb_res_dict(_dict)
{
	vector<string> fields;
	split(&fields, s, '=');
	assert(fields.size() == 2);
	label = fields[0];
	string temp = fields[1];
	split(&fields, temp, ';');
	Expand_wildcards_outside_avg(&fields);
	unsigned short i;
	num_atoms = fields.size();
	
	for(i = 0; i < fields.size(); i++)
	{
		string t = fields[i];
		vector<string> v;
		split(&v, t, ':');
		const string &def = v[0];
		if(def.substr(0, 4) == "avg(")
		{
			if(v.size() == 1) //Abort a pseudoname HAS to be provided in this case
			{
				cout << "You _have_ to provide a name for pseudotatoms that are defined " <<
					"as the geometric center of several PDB atoms. Aborting...\n";
				exit(1);
			}
			Process_avg_clause(def, i);
		}
		else
			atom2index[v[0]] = i;
		string pseudoname = v.size() == 1 ? v[0] : v[1]; //If no pseudoname has been provided use the PDB atom name
		pseudoatom_name2index[pseudoname] = i;
	}
}

void ResidueDefinition::Get_pseudoatom_names_matching_wildcard(const string &s, vector<string> *v) const
{
	v->clear();
	map<string, unsigned short>::const_iterator x;
	
	for(x = pseudoatom_name2index.begin(); x != pseudoatom_name2index.end(); x++)
		if(x->first.find(s.substr(1)) != string::npos)
			v->push_back(x->first);
}


void ResidueDefinition::Process_avg_clause(const string &def, unsigned int pseudoatom_number)
{
	assert(def[def.size() - 1] == ')');
	//This selects only the stuff inside the avg() clause
	string avg_clause = def.substr(4, def.size() - 1 - 4);
	vector<string> f;
	split(&f, avg_clause, ',');
	vector<string>::const_iterator x;
	for(x = f.begin(); x != f.end(); x++)
	{
		if((*x)[0] == '\\')
		{
			vector<string> expanded;
			/*if(label != "*")*/
			pdb_res_dict.Expand_wildcard(this->label, *x, &expanded);
			for(unsigned int k = 0; k < expanded.size(); k++)
				atom2index[expanded[k]] = pseudoatom_number;
		}
		else
			atom2index[*x] = pseudoatom_number;
	}
}	

void ResidueDefinition::Print(void)
{
	cout << "Label: " << label << "\n";
	map<string, unsigned short>::const_iterator x;
	cout << "\tAtom -> Index mapping:\n";
	for(x = atom2index.begin(); x != atom2index.end(); x++)
		cout << "\t--- " << x->first << " : " << x->second << "\n";
}
	
	
	
	
	
	
	


