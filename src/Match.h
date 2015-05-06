/***************************************************************************
 *   Copyright (C) 2006 - 2009 by Gabriele Ausiello, Pier Federico Gherardini   *
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

#ifndef _SUPERPOSE_MATCH_H
#define _SUPERPOSE_MATCH_H
#include "Chain.h"
#include "common.h"
#include "superpose.h"
#include "ChainWrapper.h"
#include <vector>

	
class Match  
{
	private:		
		float rmsd;
		std::vector<unsigned short> amino1, amino2;
		const ChainWrapper &C1, &C2;
		const options &opt;
		std::vector<double> rot_matrix;
		std::vector<double> probe_trasl, target_trasl;
		
		
				
		//void Write_set(const options &opt, const Set &s, const Chain &ch, short *amino);

	public:
		Match(const ChainWrapper &c1, const ChainWrapper &c2, const options &_opt);
		inline void Set_values(unsigned int superposition_length,
							   float _rmsd, const std::vector<unsigned short> &a1, const std::vector<unsigned short> &a2)
		{
			//The lenght of the superposition must be passed as a parameter because a1 and a2 are always
			//score_max * 2 long
			rmsd = _rmsd;
			amino1.resize(superposition_length); amino2.resize(superposition_length);
			amino1.assign(a1.begin(), a1.begin() + superposition_length);
			amino2.assign(a2.begin(), a2.begin() + superposition_length);
		}
		
	
		unsigned int size(void) const {return amino1.size(); }
		const std::vector<unsigned short> &getProbe_res(void) const {return amino1; }
		const std::vector<unsigned short> &getTarget_res(void) const {return amino2; }
		float getRmsd(void) const {return rmsd; }
		void Write(unsigned long id);
		virtual ~Match() {};
		void Superpose(void);
};


#endif

