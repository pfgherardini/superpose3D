/***************************************************************************
 *   Copyright (C) 2006 by Pier Federico Gherardini   *
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

#ifndef _SUPERPOSE_JOB_H_
#define _SUPERPOSE_JOB_H_
#include "Chain.h"
#include "common.h"
#include "ChainWrapper.h"
#include <list>
#include <vector>
#include <string>

	
class Job
{
	public:
		virtual ~Job(void);
		void Run(void);
		enum mode {normal, all_against_all, clustering, superpose_only};
		void setMode(const enum mode &m) {_mode = m;}
		Job(const std::vector<std::string> &_probes, const std::vector<std::string> &_targets, const options &_opt);

private:
		std::list<Chain *> probes;
		std::list<Chain *> targets;
		const options &opt;
		enum mode _mode;
		void Compare(ChainWrapper *c1, ChainWrapper *c2);
		void Do_superposition_only(Chain *c1, Chain *c2);
		void Load_chains(std::list<Chain *> *l, const std::vector<std::string> &v);
		void List_directory_recursively(const std::string &dirname, std::vector<std::string> *files);
		Chain *Slice_chain(const Chain &ch, const std::string &range_string);
		Chain *Select_zone_amino(const Chain &ch, const std::string &zone_string);
};


#endif

