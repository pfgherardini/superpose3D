/***************************************************************************
 *   Copyright (C) 2006 by Gabriele Ausiello *
 *   gabriele.ausiello@uniroma2.it  
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

#ifndef _SUPERPOSE_SUPERPOSE_H
#define _SUPERPOSE_SUPERPOSE_H

#include <vector>
#include "Amino.h"


struct Point
{
	double x,y,z;
	double sqrDistance_from(const Point &p2);
};

struct Set 
{
	std::vector<Point> p;
	Set(void) {}
	void Bring_to_zero(std::vector<double> *trasl);
	unsigned int size(void) const {return p.size(); }
	void Rotate(const std::vector<double> &r);
	double Rmsd(const Set &s2);
	void Add_atoms(const std::vector<Atom> &v);
};


void set_superpose(Set *target, Set *probe, std::vector<double> *r);
void jacobi (double a[5][5], int n,double  *d, double vv[5][5]);


#endif
