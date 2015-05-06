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


#ifndef _SUPERPOSE_PDBRESIDUEDICTIONARY_H
#define _SUPERPOSE_PDBRESIDUEDICTIONARY_H

#include <map>
#include <string>
#include <vector>


class PdbResidueDictionary
{
	private:
		
		struct PdbResidue
		{
			std::string parent_res;
			char one_letter_code;
			std::vector<std::string> atom_names;
		};
			
	
		static std::map<std::string, PdbResidue*> dictionary;
	public:
		static bool Is_amino(const std::string &label);
		PdbResidueDictionary(const std::string &file_name);
		void Expand_wildcard(const std::string &res, const std::string &s, std::vector<std::string> *v) const;
		typedef std::map<std::string, PdbResidue *>::const_iterator const_iterator;
		const_iterator begin(void) const {return dictionary.begin(); }
		const_iterator end(void) const {return dictionary.end(); }
};




#endif




