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

//ATTENZIONE!!! in modalita res non funziona equiv *=*


#include "SimilarityMatrix.h"
#include <cassert>
#include <cstdlib>
#include <iostream>

using std::vector;
using std::pair;
using std::map;
using std::string;
using std::make_pair;
using std::cout;


void SimilarityMatrix::Init_AA_indices(void)
{
	AAlabel2idx.insert(make_pair("ALA", 0));
	AAlabel2idx.insert(make_pair("CYS", 1));
	AAlabel2idx.insert(make_pair("ASP", 2));
	AAlabel2idx.insert(make_pair("GLU", 3));
	AAlabel2idx.insert(make_pair("PHE", 4));
	AAlabel2idx.insert(make_pair("GLY", 5));
	AAlabel2idx.insert(make_pair("HIS", 6));
	AAlabel2idx.insert(make_pair("ILE", 7));
	AAlabel2idx.insert(make_pair("LYS", 8));
	AAlabel2idx.insert(make_pair("LEU", 9));
	AAlabel2idx.insert(make_pair("MET", 10));
	AAlabel2idx.insert(make_pair("ASN", 11));
	AAlabel2idx.insert(make_pair("PRO", 12));
	AAlabel2idx.insert(make_pair("GLN", 13));
	AAlabel2idx.insert(make_pair("ARG", 14));
	AAlabel2idx.insert(make_pair("SER", 15));
	AAlabel2idx.insert(make_pair("THR", 16));
	AAlabel2idx.insert(make_pair("VAL", 17));
	AAlabel2idx.insert(make_pair("TRP", 18));
	AAlabel2idx.insert(make_pair("TYR", 19));
}
	
void SimilarityMatrix::Init_sim_matrix(void)
{
	sim_matrix.resize(20);
	for(unsigned int i = 0; i < 20; i++)
	{
		vector<pair_s> &row = sim_matrix[i];
		row.resize(20);
		for(unsigned int j = 0; j < 20; j++)
			row[j] = make_pair(-1, -1);
	}
}

void SimilarityMatrix::Resize_matrix(unsigned int new_size)
{
	assert(new_size > sim_matrix.size());
	unsigned int old_size = sim_matrix.size();
	sim_matrix.resize(new_size);
	for(unsigned int i = 0; i < new_size; i++)
	{
		vector<pair_s> &row = sim_matrix[i];
		row.resize(new_size);
		unsigned int start_new_values = i >= old_size ? 0 : old_size;
		for(unsigned int j = start_new_values; j < new_size; j++)
			row[j] = make_pair(-1, -1);
	}
}	
	
unsigned short SimilarityMatrix::Get_matrix_position_for_residue(const string &atoms, const string &label) const
{
	string key = label;
	if(atoms.size() > 0)
		key += "/" + atoms;
	typedef map<string, unsigned short>::const_iterator it;
	it x = AAlabel2idx.find(key);
	if(x == AAlabel2idx.end())
	{
		vector<it> catch_all;
		for(x = AAlabel2idx.begin(); x != AAlabel2idx.end(); x++)
		{
			const string &s = x->first;
			if(s[0] == '*')
				catch_all.push_back(x); //Collect iterators for catch all (*) positions
			else if((s.substr(0, 3) == label) && (s[4] == '\\'))
				if(atoms.find(s.substr(5)) != string::npos)
					break;
		}
		//If you still haven't found it look in catch_all (*) positions
		if(x == AAlabel2idx.end())
		{
			vector<it>::const_iterator i;
			for(i = catch_all.begin(); i != catch_all.end(); i++)
			{
				const string &s = (*i)->first;
				if(s.substr(2) == atoms)
					x = *i;
				else if(s[2] == '\\' && (atoms.find(s.substr(3)) != string::npos))
					x = *i;
			}
		}
	}
	assert(x != AAlabel2idx.end());
	return x->second;
}



unsigned short SimilarityMatrix::__get_matrix_position_for_residue(const string &atoms, const string &label)
{
	string key;
	if(opt.mode == options::mode_res)
		key = label;
	else
		key = label + "/" + atoms;
	unsigned short new_val = AAlabel2idx.size();
	pair<map<string, unsigned short>::iterator, bool> p;
	p = AAlabel2idx.insert(make_pair(key, new_val));
	idx2AAlabel.insert(make_pair(p.first->second, p.first->first));
	if(p.second) //An insertion has been performed, resize the matrix
		Resize_matrix(sim_matrix.size() + 1);
	return p.first->second;
}
	



unsigned short SimilarityMatrix::Store_atoms_vector(const string &atoms, const string &residue)
{
	map<string, unsigned short> &m = index_map[residue];
	if(opt.mode == options::mode_res)		
	{
		if(atoms[0] == '\\')
		{
			cout << "Error. You can't use wildcards (e.g.: \\N) in residue mode... aborting\n";
			exit(1);
		}
		map<string, unsigned short>::const_iterator x = m.find(atoms);
		if(x != m.end())
			return x->second;
		else
		{
			unsigned short new_idx = m.size(); //It the length before insertion is X the new index is X
			m[atoms] = new_idx;
			return new_idx;
		}
	}
	else
	{
		vector<string> v;
		if(atoms[0] == '\\' && residue != "*") //For catch_all residues (*), wildcards are expanded in ChainWrapper::Create_pseudoresidues()
			opt.pdb_res_dict->Expand_wildcard(residue, atoms, &v);
		else
			v.push_back(atoms);
		vector<string>::const_iterator atom;
		for(atom = v.begin(); atom != v.end(); atom++)
		{
			m[*atom] = 0; //In this case each atom appears as a separate entity in the pseudoresidue
		}
		return 0;
	}
}



void SimilarityMatrix::Store_indices_in_matrix(const vector_str &residues, const vector<unsigned short> &indices,
											   const vector<unsigned short> &matrix_row_indices, equiv_type eq_type)
{
	assert(residues.size() == indices.size());
	assert(residues.size() == matrix_row_indices.size());
	
	for(unsigned int i = 0; i < residues.size() - 1; i++)
	{
		for(unsigned int j = i + 1; j < residues.size(); j++)
		{
			//ATTENZIONE!! Controllare che qualcosa non fosse già presente!
			unsigned int row = matrix_row_indices[i];
			unsigned int col = matrix_row_indices[j];
			sim_matrix[row][col] = make_pair(indices[i], indices[j]);
			sim_matrix[col][row] = make_pair(indices[j], indices[i]);
		}
		//If this is a single equivalence do only the first residue against all the other ones
		//otherwise go on and do all against all
		if(eq_type == equiv_single)
			break;
	}
}

SimilarityMatrix::SimilarityMatrix(const vector<string> &lines, const options &_opt) : opt(_opt)
{
	if(!opt.pdb_res_dict)
	{
		cout << "Error! No PDB residue dictionary has been provided\n" <<
		"did you set the pdb_res_dict variable in the parameters file?\n";
		exit(1);
	}
	
	vector<string>::const_iterator s;
	for(s = lines.begin(); s != lines.end(); s++)
	{
		vector<string> fields;
		const string &temp = *s;
		split(&fields, temp, '=');
		assert(fields.size() >= 2);
		vector<string> v;
		split(&v, fields[1], ';');
		equiv_type eq_type;
		
		if(v.size() >= 2)
		{
			//Remove all the elements of fields except the first and substitute them with v
			fields.erase(fields.begin() + 1, fields.end());
			fields.insert(fields.end(), v.begin(), v.end());
			eq_type = equiv_single;
		}
		else
			eq_type = equiv_group;
		
		
		vector<string> residues;
		vector<unsigned short> indices;
		vector<unsigned short> matrix_row_indices;
		for(unsigned int i = 0; i < fields.size(); i++)
		{
			vector<string> v;
			split(&v, fields[i], '.');
			residues.push_back(v[0]);
			string atoms;
			if(v.size() == 1)
				atoms = "ALL_ATOMS";
			else
				atoms = v[1];
			unsigned short new_idx = Store_atoms_vector(atoms, v[0]);
			indices.push_back(new_idx);
			matrix_row_indices.push_back(__get_matrix_position_for_residue(atoms, v[0]));
		}
		Store_indices_in_matrix(residues, indices, matrix_row_indices, eq_type);
	}
}		

void SimilarityMatrix::Print(void)
{
	vector_str residues(AAlabel2idx.size());
	map<string, unsigned short>::const_iterator x;
	for(x = AAlabel2idx.begin(); x != AAlabel2idx.end(); x++)
		residues[x->second] = x->first;
	

	for(unsigned int i = 0; i < residues.size(); i++)
		cout << "\t" << residues[i];
	cout << "\n";
	
	
	for(unsigned int i = 0; i < residues.size(); i++)
	{
		cout << residues[i];
		for(unsigned int j = 0; j < residues.size(); j++)
			cout << "\t" << sim_matrix[i][j].first << "," <<
			sim_matrix[i][j].second;
		cout << "\n";
	}
}
	

void SimilarityMatrix::Print_index_map(void)
{
	cout << "Mapping between pseudoatom names and the indeces of coordinate vectors\n";
	idx_map_t::const_iterator x;
	
	for(x = index_map.begin(); x != index_map.end(); x++)
	{
		cout << x->first << "\n";
		map<string, unsigned short>::const_iterator i;
		for(i = x->second.begin(); i != x->second.end(); i++)
		{
			cout << "\t";
			cout << i->first;
			cout << "  -->  ";
			cout << i->second << "\n";
		}
	}
}
			
			
			
		
		
		
