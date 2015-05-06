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

#ifndef _SUPERPOSE_PDBPARSER_H
#define _SUPERPOSE_PDBPARSER_H
#include "Chain.h"
#include "ResidueDefinition.h"
#include "DefinitionCenter.h"
#include "Amino.h"
#include "common.h"
#include <map>
#include <string>
#include <cstdio>

class PdbParser
{
	private:
		const options &opt;
		const DefinitionCenter &def_center;
		const std::string &file_name;
		std::vector<std::string> lines;
		
		inline std::string Get_chain(const std::string &line);
		void Handle_amino_end(Amino **a, Chain *c);
		unsigned int Find_chain_start(const std::string &chain);
		void Calculate_chain_offsets(const std::string &chain, std::vector<float> *offsets, unsigned int *chain_start);
		Atom *Create_atom(const std::string &s, const std::vector<float> &offsets);
		bool Check_amino_integrity(const std::vector<std::vector<const Atom *> > &m);
		Amino *Convert_pdb_amino(const Amino &a);
		inline bool Chain_finished(const std::string &cur_chain, const std::string &line);
		void Read_lines(const std::string &input_file);
	
	public:
	
		Chain *Load_chain(const std::string &chain_id);
		void Get_chain_list(std::vector<std::string> *v);	
		PdbParser(const options &_opt, const std::string &_file_name) : 
			opt(_opt), def_center(*(opt.def_center)), file_name(_file_name) {}
		
		~PdbParser(void) {}
		void Flush_lines(void) { lines.clear(); }
};

#endif
