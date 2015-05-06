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

#ifndef _SUPERPOSE_SIMILARITYMATRIX_H
#define _SUPERPOSE_SIMILARITYMATRIX_H


#include <vector>
#include <map>
#include <string>
#include "common.h"


//The matrix is as follows. Row and columns are residues.
//Each element of the matrix is a pair of indices which
//refers to the vector of pseudoatoms for each residues.
//In case of custom pairing rules
//e.g
//ALA=TYR
//ALA.CB=SER.CB
//
//you would have
//
//		TYR												SER
//ALA [pair of indeces of array of coordinates] [pair of indices of vectors containing ALA.CB and SER.CB]
//Therefore in this case ALA must be made of two vectors
//1) Coordinates used when pairing with any residue except SER
//2) Coordinates used when pairing with SER
//
//These definitions are stored in ChainWrappers which are objects that wrap a Chain,
//preparing it for the comparison algorithm

struct options; //ATTENZIONE!! QUESTO SI POTRA' TOGLIERE UNA VOLTA RISOLTO IL PROBLEMA DEI MODE

class SimilarityMatrix
{
	private:
		typedef std::pair<short, short> pair_s; //-1 signals that the index hasn't been defined
		typedef std::vector<std::string> vector_str;
	
		//idx_map_t a structure that describes to ChainWrappers how the coordinates of each residue
		//type should be arranged. The organizations is as follows
		//
		//	ALA -> [ALL_ATOMS] : 0
		//		-> [CB] : 1
		//		-> [CA,CB] : 2
		//
		//Means that ALA residue should have at index 0 a vector containing all the atoms, 
		//at index 1 a vector containing the CB and at index 2 a vector containing CA and CB
		//The map is indexed by reside name and the value is a map indexed by vector of strings.
		//The value of this last map is the index in the ChainWrapper coordinates while the vector
		//of strings gives the names of the atoms that correspond to that index.
	
		enum equiv_type {equiv_group, equiv_single};
		typedef std::map<std::string, std::map<std::string, unsigned short> > idx_map_t;
		idx_map_t index_map;
		std::vector<std::vector<pair_s> > sim_matrix;
		std::map<std::string, unsigned short> AAlabel2idx;
		std::map<unsigned short, std::string> idx2AAlabel;
		const options &opt;
	
		void Init_AA_indices(void);
		void Init_sim_matrix(void);
		unsigned short Store_atoms_vector(const std::string &atoms, const std::string &residue);
		void Store_indices_in_matrix(const vector_str &residues, const std::vector<unsigned short> &indices,
									 const std::vector<unsigned short> &matrix_row_indices, equiv_type eq_type);
		void Resize_matrix(unsigned int new_size);
		unsigned short __get_matrix_position_for_residue(const std::string &atoms, const std::string &label);
	
	
	public:
		void Print(void);
		void Print_index_map(void);
		const std::vector<pair_s> &operator[](unsigned int i) const {return sim_matrix[i]; }
		SimilarityMatrix(const vector_str &lines, const options &_opt);
		unsigned short Get_matrix_position_for_residue(const std::string &atoms, const std::string &label) const;
		
		
		inline const std::map<std::string, unsigned short> *Get_indices_for_residue(const std::string &s) const
		{
			const std::map<std::string, unsigned short> *ret;
			idx_map_t::const_iterator x = index_map.find(s);
			if(x != index_map.end())
				ret = &(x->second);
			else if((x = index_map.find("*")) != index_map.end())
				ret = &(x->second);
			else
				ret = 0;
			return ret;
		}
		const std::string &Idx2AAlabel(unsigned short idx) const
		{
			std::map<unsigned short, std::string>::const_iterator x;
			x = idx2AAlabel.find(idx);
			assert(x != idx2AAlabel.end());
			return x->second;
		}
		
};

#endif




