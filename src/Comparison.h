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

#ifndef _SUPERPOSE_COMPARISON_H
#define _SUPERPOSE_COMPARISON_H


#include "Superposition.h"
#include "Chain.h"
#include "ChainWrapper.h"
#include "Match.h"
#include "Tree.h"
#include "common.h"
#include <list>
#include "SimilarityMatrix.h"

class Comparison  
{
	private:
		
		unsigned int score_min, score_max;
		
		float rmsd_max;

		Superposition super;
		ChainWrapper *C1;
		ChainWrapper *C2;
		std::list<Match *> *match_list;
		unsigned int superposition_length;
		std::vector<unsigned short> amino1, amino2; //amino1 = probe; amino2 = target
		const options &opt;
		Tree tree;


		void Compare(short n1,short n2, int *max_longest_matches);
		static void (Comparison::* Add_match_to_list)(float rms);
		void Add_match_keeping_all(float rms);
		void Add_match_keeping_best(float rms);
	public:
		
		Comparison(const options &opt);
		static void setKeep_best_matches(void) {Add_match_to_list = &Comparison::Add_match_keeping_best; }
		static void setKeep_all_matches(void) {Add_match_to_list = &Comparison::Add_match_keeping_all; }
		virtual ~Comparison();
		void Start(ChainWrapper *_c1, ChainWrapper *_c2, std::list<Match *> *m_list);
};

#endif

	
