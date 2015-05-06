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


#ifndef _SUPERPOSE_CHAINWRAPPER_H
#define _SUPERPOSE_CHAINWRAPPER_H

#include <vector>
#include "Amino.h"
#include "SimilarityMatrix.h"
#include "DefinitionCenter.h"
#include "Chain.h"
#include "common.h"
#include <iostream>
		
struct PseudoResidue
{
	std::vector<std::vector<Atom> > coord_cache;
	unsigned short idx_in_matrix;
	unsigned short idx_used; //Index in the residue_used_vector
	unsigned short idx_orig_residue;
	std::string orig_atom_name;
	
	inline std::vector<Atom> &operator[](unsigned short i) {return coord_cache[i];}
	inline const std::vector<Atom> &operator[](unsigned short i) const {return coord_cache[i];}
	unsigned int size(void) const { return coord_cache.size(); }
	void Load(const Amino &amino, const std::map<std::string, unsigned short> *idx_map, 
			  const ResidueDefinition &def, const options &opt);
	inline unsigned short getIdx_in_matrix(void) const {return idx_in_matrix; }

	bool Is_neighbour(const PseudoResidue &b, float neighbour_threshold);
};


class ChainWrapper
{
	private:

		std::vector<PseudoResidue> pseudochain;
		std::vector<std::vector<unsigned short> >neighbours;
		std::vector<char> residue_used;
	
		
		const Chain &parent_chain;
		const options &opt;
		void Create_pseudoresidues(void);
		void Create_single_pseudoresidue(const std::string &pseudoatom, unsigned int amino_idx, const ResidueDefinition &def);
		void Do_neighbours(void);
	
		
	
	public:
		ChainWrapper(const Chain &ch, const options &_opt);
		void Print(void);
		inline char &Residue_used(unsigned short pseudores_idx)
			{
				return residue_used[pseudochain[pseudores_idx].idx_used];
			}
		
		inline const std::vector<std::vector<unsigned short> > &Neighbours(void) const {return neighbours; }
		unsigned int size(void) const {return pseudochain.size(); }
		PseudoResidue &operator[](unsigned int i) {return pseudochain[i]; }
		const PseudoResidue &operator[](unsigned int i) const {return pseudochain[i]; }
		const Chain &getParent_chain(void) const {return parent_chain; }
		unsigned int getEffective_size(void) const;
		void Print_used(void) const
		{
			std::cout << parent_chain.getFull_id() << " --- ";
			for(unsigned int i = 0; i < residue_used.size(); i++)
				std::cout << (int)(residue_used[i]);
			std::cout << "\n";
		}
	
		void Print_neighbours(void) const
		{
			for(unsigned int i = 0; i < neighbours.size(); i++)
			{
				std::cout << parent_chain.getFull_id() << " -- " << i << " /";
				for(unsigned int j = 0; j < neighbours[i].size(); j++)
					std::cout << " " << neighbours[i][j];
				std::cout << "\n";
			}
		}
		
};



#endif

