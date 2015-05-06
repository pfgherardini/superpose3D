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



#ifndef _SUPERPOSE_RESIDUE_DEFINITION_H
#define _SUPERPOSE_RESIDUE_DEFINITION_H

#include <map>
#include <string>
#include <cassert>
#include <vector>
#include "PdbResidueDictionary.h"

class ResidueDefinition
{
	private:
		//This maps a PDB atom label to the index in the vector
		//of coordinates that define the current residue
		std::map<std::string, unsigned short> atom2index;
		
		//This maps the pseudoatom name to the index
		//in the vector of coordinates.
		//Pseudo atom names are defined as follows
		//
		//def ALA.CB:foo
		//
		//defines "foo" as the name of the pseudoatom representing the CB
		 
		std::map<std::string, unsigned short> pseudoatom_name2index;
		const PdbResidueDictionary &pdb_res_dict;
		std::string label;
		unsigned short num_atoms;
		void Expand_wildcards_outside_avg(std::vector<std::string> *res);
		void Process_avg_clause(const std::string &def, unsigned int pseudoatom_number);
		short Find_atom_with_regex(const std::map<std::string, unsigned short> &m, const std::string &key) const;
		

	
	public:
		//Residue defintions look like this 
		//def ALA=CB;avg(N,CO)
		ResidueDefinition(const std::string &s, const PdbResidueDictionary &_dict);
		unsigned short getNum_atoms(void) const { return num_atoms; }
		const std::string &getLabel(void) const { return label; }
		void Print(void);
		short Atom2index(const std::string &s) const;
		void Get_pseudoatom_names_matching_wildcard(const std::string &s, std::vector<std::string> *v) const;
		//ATTENZIONE!!! CONTROLLARE BENE COSA SUCCEDE QUI IN CASO DI ATOMI MANCANTI 
		unsigned short Pseudoatom_name2index(const std::string &s) const;
	
};



#endif
