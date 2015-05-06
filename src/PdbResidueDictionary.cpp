/***************************************************************************
 *   Copyright (C) 2009 by Pier Federico Gherardini   *
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



#include "PdbResidueDictionary.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "common.h"

using std::map;
using std::vector;
using std::string;
using std::ifstream;
using std::cout;


map<string, PdbResidueDictionary::PdbResidue*> PdbResidueDictionary::dictionary;

bool PdbResidueDictionary::Is_amino(const string &line)
{
	if(line.substr(0, 4) == "ATOM" ||
	   line.substr(0, 6) == "HETATM")
	{
		map<string, PdbResidue*>::const_iterator x = dictionary.find(line.substr(17, 3));
		if(x != dictionary.end())
			return true;
		else
			return false;
	}
	else
		return false;
}
	

void PdbResidueDictionary::Expand_wildcard(const string &res, const string &s, vector<string> *v) const
{
	v->resize(0);
	string _s = s.substr(1); //Discard the initial '\' character
	map<string, PdbResidue*>::const_iterator x = dictionary.find(res);
	assert(x != dictionary.end());
	const vector<string> &atom_names = x->second->atom_names;

	for(unsigned int i = 0; i < atom_names.size(); i++)
		if(atom_names[i].find(_s) != string::npos)
			v->push_back(atom_names[i]);
}

PdbResidueDictionary::PdbResidueDictionary(const string &file_name)
{
	ifstream in(file_name.c_str());
	if(in.rdstate() & ifstream::failbit)
	{
		cout << "Error! Can't open the residue dictionary in " << file_name << "\n" <<
			"please set the pdb_res_dict option in the parameters file to a correct value\n";
		exit(1);
	}
	
	string temp;
	
	while(getline(in, temp) && !in.eof())
	{
		vector<string> fields;
		split(&fields, temp, '\t');
		PdbResidue *res = new PdbResidue;
		res->parent_res = fields[1];
		res->one_letter_code = fields[2][0];
		split(&(res->atom_names), fields[3], ',');
		dictionary[fields[0]] = res;
	}
	in.close();
	
}
	
	
	
