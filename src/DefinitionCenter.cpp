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
#include <iostream>
#include <set>
#include <cstdlib>

using std::map;
using std::string;
using std::vector;
using std::cout;
using std::set;

DefinitionCenter::~DefinitionCenter(void)
{
	defmap_t::iterator x;
	for(x = definitions.begin(); x != definitions.end(); x++)
		if(x->second)
			delete x->second;
	definitions.clear();
	/*if(catch_all)
		delete catch_all;*/
}



void DefinitionCenter::Expand_wildcard_residue_names(vector<string> *lines)
{
	vector<string>::iterator i;
	bool do_expand = false;
	set<string> labels;
	string def_to_expand;
	for(i = lines->begin(); i != lines->end(); i++)
	{
		vector<string> fields;
		split(&fields, *i, '=');
		if(fields[0] != "*")
			labels.insert(fields[0]);
		else
		{
			def_to_expand = fields[1];
			do_expand = true;
			*i = "";
		}
		
	}
	if(do_expand)
	{
		PdbResidueDictionary::const_iterator x;
		for(x = opt.pdb_res_dict->begin(); x != opt.pdb_res_dict->end(); x++)
			if(labels.find(x->first) == labels.end())
				lines->push_back(x->first + "=" + def_to_expand);
	}
}
											   
											   
											   

DefinitionCenter::DefinitionCenter(const vector<string> &_lines, const options &_opt) : /*catch_all(0),*/ opt(_opt)
{
	if(!opt.pdb_res_dict)
	{
		cout << "Error! No PDB residue dictionary has been provided\n" <<
		"did you set the pdb_res_dict variable in the parameters file?\n";
		exit(1);
	}
	
	
	vector<string> lines = _lines;
	Expand_wildcard_residue_names(&lines);
	vector<string>::const_iterator i;
	for(i = lines.begin(); i != lines.end(); i++)
	{
		if(i->size())
		{
			ResidueDefinition *def = new ResidueDefinition(*i, *(opt.pdb_res_dict));
			/*if(def->getLabel() == "*")
				catch_all = def;*/
			//The definition was already there
			if(!(definitions.insert(std::make_pair(def->getLabel(), def)).second))
			{
				cout << "Multiple definitions found " <<
					"for residue " << def->getLabel() << " aborting...\n";
				exit(1);
			}
		}
	}
}

const ResidueDefinition &DefinitionCenter::Get_definition(const string &label) const
{
	defmap_t::const_iterator x = definitions.find(label);
	ResidueDefinition *retval = 0;
	if(x != definitions.end())
		retval = x->second;
	else
	{
		/*if(catch_all)
			retval = catch_all;*/
		
		cout << "No definition found " <<
			"for residue " << label << " and no \"catch all\" " <<
			"definition provided, aborting...\n";
		exit(1);
		
	}
	return *retval;
}
		
	
	

