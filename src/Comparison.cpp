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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <iostream>

#include "Comparison.h"
#include "Superposition.h"
#include <map>


using std::list;
using std::pair;
using std::vector;
using std::cout;

void (Comparison::* Comparison::Add_match_to_list)(float rms) = &Comparison::Add_match_keeping_best;

/*
Comparison::Comparison(int _score_max, float _rmsd_max, const SimilarityMatrix &_sim_mat) : amino1(_score_max), 
						amino2(_score_max), sim_mat(_sim_mat)
{
	score_max = _score_max;
	rmsd_max = _rmsd_max;	
	tree.Init(score_max);
}
 */

Comparison::Comparison(const options &_opt) : amino1(_opt.score_max), 
						amino2(_opt.score_max), opt(_opt)
{
	score_max = opt.score_max;
	score_min = opt.score_min;
	rmsd_max = opt.rmsd_thresh;
	tree.Init(score_max);
}

Comparison::~Comparison() 
{
	amino1.clear();
	amino2.clear();
	tree.Clear();
}

void Comparison::Start (ChainWrapper *_c1, ChainWrapper *_c2, list<Match *> *m_list)
{	
	C1 = _c1; C2 = _c2; 
	superposition_length = 0;
	match_list = m_list;
	int max_longest_matches = opt.max_longest_matches;
	
	unsigned short a1,a2;
	for (a1 = 0; a1 < C1->size(); a1++) 
		if (!(C1->Residue_used(a1))) 
			for (a2 = 0; a2 < C2->size(); a2++) 
				if (!(C2->Residue_used(a2))) 
					Compare(a1, a2 , &max_longest_matches);

	if(!(max_longest_matches)) 
		cout << "Warning! the maximum number of longest matches has been reached\n" <<
		"when comparing " << _c1->getParent_chain().getFull_id() << " with " <<
		_c2->getParent_chain().getFull_id() << "\n" <<
		"The output you are seeing is probably not what you want.\n" <<
		"Please refer to the documentation for additional details.\n";
	
	tree.Clear();
	amino1.clear();
	amino2.clear();
	for(list<Match *>::iterator i = m_list->begin(); i != m_list->end(); i++)
		(*i)->Superpose();
}

void Comparison::Add_match_keeping_best(float rmsd)
{
	Match *temp = new Match(*C1, *C2, opt);
	temp->Set_values(superposition_length, rmsd, amino1, amino2);
	if(!match_list->size() || match_list->front()->size() == superposition_length)
		match_list->push_back(temp);
	else if(superposition_length > match_list->front()->size())
	{
		list<Match *>::iterator i;
		for(i = match_list->begin(); i != match_list->end(); i++)
			delete *i;
		match_list->clear();
		match_list->push_back(temp);
	}
	else
		delete temp;
}

void Comparison::Add_match_keeping_all(float rmsd)
{
	Match *temp = new Match(*C1, *C2, opt);
	temp->Set_values(superposition_length, rmsd, amino1, amino2);
	match_list->push_back(temp);
}


void Comparison::Compare(short n1,short n2, int *max_longest_matches)
{	
	if(!(*max_longest_matches)) return;
	unsigned short idx1 = (*C1)[n1].getIdx_in_matrix();
	unsigned short idx2 = (*C2)[n2].getIdx_in_matrix();
	pair<short, short> p = (*opt.sim_mat)[idx1][idx2];
	if(p.first == -1) //Residues are not comparable
		return;
	//Get the relevant atom vector
	const vector<Atom> &atoms1 = (*C1)[n1][p.first];
	const vector<Atom> &atoms2 = (*C2)[n2][p.second];
	assert(atoms1.size() == atoms2.size());
	Superposition c(atoms1, atoms2);
	super.Add(&c);
	float rmsd = 0.0;
	if (superposition_length > 0) 
	{
		rmsd = super.Rms();	
		if(rmsd > rmsd_max * 100) 
		{
			super.Sub(&c);
			return;
		}
	}
	amino1[superposition_length] = n1;
	amino2[superposition_length] = n2;
	C1->Residue_used(n1) = 1;
	C2->Residue_used(n2) = 1;
	int to_remove = tree.Add(n1,n2);
	superposition_length++;	

	if (!tree.Done()) 
	{
		if(superposition_length >= score_min) (this->*Add_match_to_list)(rmsd);
		if(superposition_length < score_max)
			for (unsigned int j = 0; j < superposition_length; j++) 
			{
				const vector<unsigned short> &neigh1 = C1->Neighbours()[amino1[j]];
				vector<unsigned short>::const_iterator v1;
				
				for (v1 = neigh1.begin(); v1 != neigh1.end(); v1++) 
					if (!(C1->Residue_used(*v1)))
					{
						//Get the neighbours of the amino acid in C2 matching amino1[j]
						//const vector<unsigned short> &neigh2 = C2->Neighbours()[amino2[S1->used[amino1[j]]]];
						const vector<unsigned short> &neigh2 = C2->Neighbours()[amino2[j]];
						vector<unsigned short>::const_iterator v2;
						for(v2 = neigh2.begin(); v2 != neigh2.end(); v2++) 
							if (!(C2->Residue_used(*v2))) 
								Compare(*v1, *v2, max_longest_matches);
					}
			}
		else
			(*max_longest_matches)--;
	}
	superposition_length--;
	C1->Residue_used(n1) = 0;
	C2->Residue_used(n2) = 0;
	tree.Remove(to_remove);
	super.Sub(&c);
}

