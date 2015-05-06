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



#ifndef _SUPERPOSE_DEFINITIONCENTER_H
#define _SUPERPOSE_DEFINITIONCENTER_H


#include "ResidueDefinition.h"
#include <map>
#include <string>
#include <vector>
#include "common.h"


struct options;

class DefinitionCenter
{
	private:
		typedef std::map<std::string, ResidueDefinition*> defmap_t;
		defmap_t definitions;
		//ResidueDefinition *catch_all;
		const options &opt;
		void Expand_wildcard_residue_names(std::vector<std::string> *lines);
	
	public:
		virtual ~DefinitionCenter(void);
		DefinitionCenter(const std::vector<std::string> &_lines, const options &_opt);
		const ResidueDefinition &Get_definition(const std::string &label) const;
	
};
	






#endif
