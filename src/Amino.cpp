/***************************************************************************
 *   Copyright (C) 2007-2009 by Pier Federico Gherardini   *
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


#include "Amino.h"
#include "Quaternion.h"
#include "common.h"
#include <iostream>


using std::cout;
using std::string;
using std::vector;


Amino::Amino(const Amino &a) : atoms(a.atoms.size()), pdb_resnum(a.pdb_resnum),
		one_letter_code(a.one_letter_code), alt_loc(a.alt_loc), label(a.label)
{
	for(unsigned int i = 0; i < atoms.size(); i++)
		if(a.atoms[i])
			atoms[i] = new Atom(*(a.atoms[i]));
}


void Amino::Reconstruct_cb(void)
{
	//ATTENZIONE!!!!! QUI STIAMO USANDO INTERI INVECE CHE FLOAT!!!!
	
	//L'angolo dovrebbe essere +120 in quanto un angolo positivo descrive una rotazione in
	//senso orario ponendosi sull'origine e guardando lungo l'asse. Siccome il nostro asse
	//e' CA-->N la rotazione che ci serve a noi e' in senso orario
	
	const Atom *n, *ca, *c;
	n = ca = c = 0;
	for(unsigned int i = 0; i < atoms.size(); i++)
	{
		const string &lab = (atoms[i])->label;
		if(lab == "N")
			n = atoms[i];
		else if(lab == "CA")
			ca = atoms[i];
		else if(lab == "C")
			c = atoms[i];
	}
	
	if(n && ca && c)
	{
		Quaternion q((float)n->x - (float)ca->x, (float)n->y - (float)ca->y,
					 (float)n->z - (float)ca->z, 120);
		float _c[3] = {(float)c->x, (float)c->y, (float)c->z};
		float _ca[3] = {(float)ca->x, (float)ca->y, (float)ca->z};
		float res[3];
		q.Rotate(_c, _ca, res);
		atoms.push_back(new Atom("CB", (int)res[0], (int)res[1], (int)res[2]));
	}
}


	
bool Amino::IsAmino(const string &line)
{
	if(line.substr(0, 4) == "ATOM")
	{
		if(trim_copy(line.substr(17, 3)).size() == 3)
			return true;
		else
			return false;
	}
	//ATTENZIONE!! Bisogna risolvere il caso degli amminoacidi che iniziano con HETATM
/*	else if(line.substr(0, 6) == "HETATM")
		return aa_index.find(line.substr(17, 3))
		!= aa_index.end();*/
	else
		return false;
}

bool Amino::Is_neighbour(const Amino &b, float neigh_thresh) const
{
	for(unsigned int i = 0; i < atoms.size(); i++)
		for(unsigned int j = 0; j < b.atoms.size(); j++)
			if((atoms[i])->sqrDistance_from(*(b.atoms[j])) < (neigh_thresh * neigh_thresh * 100 * 100))
				return true;
	return false;
}



void Amino::Load_from_pdb_line(const string &s)
{
	label = s.substr(17, 3);
	pdb_resnum = trim_copy(s.substr(22, 5));//Paste together resnum and insertion code
	//ATTENZIONE QUI VA RISOLTA LA COSA DEL ONE_LETTER_CODE
	//one_letter_code = aa_index[label];
	alt_loc = ' ';
}

bool Amino::Add_atom(Atom *a)
{
	if(a->label != "OXT")
	{
		atoms.push_back(a);
		return true;
	}
	else
		return false;
}

Amino::~Amino(void)
{
	vector<Atom *>::iterator i;
	for(i = atoms.begin(); i != atoms.end(); i++)
	{
		if(*i)
			delete *i;
	}
}


