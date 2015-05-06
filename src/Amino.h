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

#ifndef _SUPERPOSE_AMINO_H
#define _SUPERPOSE_AMINO_H


#include <vector>
#include <string>
#include <cassert>



struct Atom
{
	std::string label;
	int x, y, z;
	Atom(void) : x(0), y(0), z(0) {}
	Atom(const std::string &lab, int _x, int _y, int _z) :
		label(lab), x(_x), y(_y), z(_z) {}
	float sqrDistance_from(const Atom &b) const
	{
		return (float)((x - b.x) * (x - b.x) +
		(y - b.y) * (y - b.y) +
		(z - b.z) * (z - b.z));
	}
	
	Atom &operator=(const Atom &a)
	{
		this->x = a.x; this->y = a.y; this->z = a.z;
		this->label = a.label;
		return *this;
	}
};

class Amino 
{
	
	private:
		std::vector<Atom *> atoms;
		std::string pdb_resnum;
		char one_letter_code;
		char alt_loc;
		std::string label;
		
public:
		const std::vector<Atom *> &getAtoms(void) const { return atoms; }
		Atom *const &operator[](unsigned int i) const { return atoms[i]; }
		Atom *&operator[](unsigned int i) { return atoms[i]; }
		~Amino(void);
		Amino(const Amino &a);
		Amino(void) : one_letter_code(' ') {}
		Amino(unsigned short _num_atoms) : atoms(_num_atoms), one_letter_code(' ') 
		{
			/*for(unsigned int i = 0; i < atoms.size(); i++)
				atoms[i] = new Atom;*/
		}
		unsigned short getNum_atoms(void) const { return atoms.size(); }
		void setPdb_resnum(const std::string &s) { pdb_resnum = s; }
		void setLabel(const std::string &s) { label = s; }
		const std::string &getPdb_resnum(void) const { return pdb_resnum; }
		const std::string &getLabel(void) const { return label; }
		char getAlt_loc(void) const { return alt_loc; }
		void setAlt_loc(char c) { alt_loc = c; }
		char getOne_letter_code(void) const { return one_letter_code; }
		bool Add_atom(Atom *a);
		void Load_from_pdb_line(const std::string &s);
		static bool IsAmino(const std::string &line);
		bool Is_neighbour(const Amino &b, float neigh_thresh) const;
		void Reconstruct_cb(void);
};



#endif
