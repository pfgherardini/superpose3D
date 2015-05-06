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




#include "ChainWrapper.h"
#include <cassert>
#include <iostream>

using std::vector;
using std::cout;
using std::string;
using std::map;



bool Fill_atom_vector_from_string(vector<Atom> *atom_v, const string &s, 
							 const Amino &amino, const ResidueDefinition &def);
inline bool Copy_all_atoms(std::vector<Atom> *dest, const Amino &amino);

inline bool Copy_all_atoms(vector<Atom> *dest, const Amino &amino)
{
	const vector<Atom *> &amino_atoms = amino.getAtoms();
	dest->resize(amino_atoms.size());
	for(unsigned int i = 0; i < amino_atoms.size(); i++)
	{
		if(amino_atoms[i])
			(*dest)[i] = *(amino_atoms[i]);
		else
			return false;
	}
	return true;
}


bool Fill_atom_vector_from_string(vector<Atom> *atom_v, const string &s, 
												 const Amino &amino, const ResidueDefinition &def)
{
	if(s == "ALL_ATOMS")
		return Copy_all_atoms(atom_v, amino);
	else
	{
		vector<string> str_v;
		split(&str_v, s, '-');
		atom_v->resize(str_v.size());
		for(unsigned int j = 0; j < str_v.size(); j++)
		{
			unsigned short idx = def.Pseudoatom_name2index(str_v[j]);
			assert(idx < amino.getAtoms().size());
			if(amino[idx])
				(*atom_v)[j] = *(amino[idx]); 
			else
				return false;
		}
	}
	return true;
}	
	

void PseudoResidue::Load(const Amino &amino, const map<string, unsigned short> *idx_map,
						 const ResidueDefinition &def, const options &opt)
{
	assert(opt.mode == options::mode_res);
	//Since we are in single residue mode the entries in the substitution matrix
	//are labelled with the residue name only
	idx_in_matrix = opt.sim_mat->Get_matrix_position_for_residue("", amino.getLabel()); 
	coord_cache.resize(idx_map->size());
	map<string, unsigned short>::const_iterator x;
	for(x = idx_map->begin(); x != idx_map->end(); x++)
	{
		assert(x->second < coord_cache.size());
		vector<Atom> &atom_v = coord_cache[x->second];
		bool filled = Fill_atom_vector_from_string(&atom_v, x->first, amino, def);
		assert(filled);
	}
}

unsigned int ChainWrapper::getEffective_size(void) const
{
	if(opt.mode == options::mode_atom)
		return pseudochain.size();
	else
		return parent_chain.size();
}

bool PseudoResidue::Is_neighbour(const PseudoResidue &b, float neighbour_threshold)
{
	for(unsigned int i = 0; i < coord_cache.size(); i++)
		for(unsigned int j = 0; j < b.coord_cache.size(); j++)
		{
			const vector<Atom> &v1 = coord_cache[i];
			const vector<Atom> &v2 = b.coord_cache[j];
			
			for(unsigned int n = 0; n < v1.size(); n++)
				for(unsigned int p = 0; p < v2.size(); p++)
					if(v1[n].sqrDistance_from(v2[p]) < 
					   (neighbour_threshold * neighbour_threshold * 100 * 100))
						return true;
		}
	return false;
}

void ChainWrapper::Do_neighbours(void)
{
	neighbours.resize(pseudochain.size());
	for(unsigned int i = 0; i < pseudochain.size() - 1; i++)
		for(unsigned int j = i + 1; j < pseudochain.size(); j++) 
			if(pseudochain[i].Is_neighbour(pseudochain[j], opt.neighbour_threshold))
			{
				neighbours[i].push_back(j);
				neighbours[j].push_back(i);
			}
}


void ChainWrapper::Create_single_pseudoresidue(const string &pseudoatom, unsigned int amino_idx, const ResidueDefinition &def)
{
	assert(opt.mode != options::mode_res);
	const Amino &amino = parent_chain[amino_idx];
	pseudochain.resize(pseudochain.size() + 1);
	PseudoResidue &p = pseudochain.back();
	p.coord_cache.resize(1);
	p.idx_in_matrix = opt.sim_mat->Get_matrix_position_for_residue(pseudoatom, amino.getLabel());
	p.idx_orig_residue = amino_idx;
	p.orig_atom_name = pseudoatom;
	
	vector<Atom> &atom_v = p[0];
	if(Fill_atom_vector_from_string(&atom_v, pseudoatom, amino, def))
	{
		if(opt.mode == options::mode_atom)
			p.idx_used = pseudochain.size() - 1;
		else if(opt.mode == options::mode_res_multiple)
			p.idx_used = amino_idx;
	}
	else //Discard the pseudoresidue
		pseudochain.resize(pseudochain.size() - 1);
}	
	

void ChainWrapper::Create_pseudoresidues()
{
	typedef map<string, unsigned short> map_t;
	for(unsigned int i = 0; i < parent_chain.size(); i++)
	{
		const Amino &amino = parent_chain[i];
		const ResidueDefinition &def = opt.def_center->Get_definition(amino.getLabel());
		const map_t *idx_map = opt.sim_mat->Get_indices_for_residue(amino.getLabel());
		//Residues which are not in the similarity matrix do not appear in the ChainWrapper
		if(idx_map)
		{
			if(opt.mode == options::mode_res)
			{
				pseudochain.resize(pseudochain.size() + 1);
				PseudoResidue &pseudores = pseudochain.back();
				pseudores.Load(amino, idx_map, def, opt);
				pseudores.idx_used = pseudochain.size() - 1;
				pseudores.idx_orig_residue = i;
			}
			else
			{
				map_t::const_iterator x;
				
				for(x = idx_map->begin(); x != idx_map->end(); x++)
				{
					assert(x->second == 0); //In this mode each pseudores corresponds to a single coordinate vector
					vector<string> v;
					if((x->first)[0] != '\\')
						v.push_back(x->first);
					else
						def.Get_pseudoatom_names_matching_wildcard(x->first, &v);
					vector<string>::const_iterator s;
					for(s = v.begin(); s != v.end(); s++)
						Create_single_pseudoresidue(*s, i, def);
				}
			}
		}
	}
}

/*
void ChainWrapper::Create_pseudoresidues(const SimilarityMatrix &sim_mat, 
										 const DefinitionCenter &def_center)
{
	typedef map<vector<string>, unsigned short> map_t;
	for(unsigned int i = 0; i < parent_chain.size(); i++)
	{
		const string &reslabel = parent_chain[i].getLabel();
		const ResidueDefinition &def = def_center.Get_definition(reslabel)
		const map_t *idx_map = sim_mat.Get_indices_for_residue(reslabel);
		if(idx_map)
		{
			PseudoResidue &pseudores = pseudochain[i];
			pseudores.Load(parent_chain[i], idx_map, sim_mat, def);
		}
	}
}	
*/

void ChainWrapper::Print(void)
{
	for(unsigned int i = 0; i < pseudochain.size(); i++)
	{
		cout <<  parent_chain[pseudochain[i].idx_orig_residue].getPdb_resnum() << " " <<
			parent_chain[pseudochain[i].idx_orig_residue].getLabel() << " " <<
			pseudochain[i].orig_atom_name << "\n";
		const PseudoResidue &res = pseudochain[i];
		for(unsigned int j = 0; j < res.size(); j++)
		{
			cout << "\t";
			for(unsigned int x = 0; x < res[j].size(); x++)
			{
				if(x)
					cout << "  ";
				cout << res[j][x].x << "," <<  res[j][x].y << "," << res[j][x].z;
			}
			cout << "\n";
		}
	}
}
		

ChainWrapper::ChainWrapper(const Chain &ch, const options &_opt) : parent_chain(ch), opt(_opt)
{
	if(opt.verbose)
		cout << "Creating Chainwrapper " << ch.getFull_id() << "\n";
	Create_pseudoresidues();
	if(opt.mode == options::mode_atom)
		residue_used.resize(pseudochain.size());
	else
		residue_used.resize(parent_chain.size());
	residue_used.assign(residue_used.size(), 0);
	Do_neighbours();
}
	

