/***************************************************************************
 *   Copyright (C) 2006 by Gabriele Ausiello, Pier Federico Gherardini   *
 *   gabriele.ausiello@uniroma2.it
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

#ifndef _SUPERPOSE_CHAIN_H
#define _SUPERPOSE_CHAIN_H
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "Amino.h"



class Chain  
{
	private:
		std::string file_name;
		std::string chain_id;
		std::vector<float> offsets;
		
		void Do_neighbours(void);
	
	public:

	
	
		std::vector<Amino *> residues;
	
				
		Chain(const Chain &ch, const std::vector<bool> &sel);
		Chain(const Chain &ch, const std::vector<unsigned int> &indices);
		Chain(const char *_id, unsigned int size);
		Chain(const char *_id, unsigned int size, const std::string &annot);
		Chain(const std::string &_file_name, const std::string &_chain_id) : file_name(_file_name), chain_id(_chain_id) {}
		virtual ~Chain();
		void Clear(void);
		
		
		bool Load(const Chain &ch, const std::vector<bool> &sel);

	
		const Amino &operator[](int n) const {return *(residues[n]); }
		int getIndex(const std::string &resnum) const;
		inline unsigned int size() const {return residues.size(); }
		const std::string &getChain_id(void) const {return chain_id; }
		const std::string getFull_id(void) const {return file_name + chain_id; }
		void setOffsets(const std::vector<float> &off) {offsets = off; }
		const std::vector<float> &getOffsets(void) const {return offsets; }
		void Add_residue(Amino *amino) {residues.push_back(amino); }
		void Write(const std::string &file_name);
};


#endif 
