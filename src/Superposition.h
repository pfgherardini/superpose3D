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
 
 
#ifndef _SUPERPOSE_SUPERPOSITION_H
#define _SUPERPOSE_SUPERPOSITION_H

#include "Amino.h"
#include <vector>
	
class Superposition 
{
	private:
		int na;
		long x1,y1,z1,x1x1,y1y1,z1z1,
			x2,y2,z2,x2x2,y2y2,z2z2,
			x1x2,x1y2,x1z2,
			y1x2,y1y2,y1z2,
			z1x2,z1y2,z1z2;
	public:
		Superposition(void);
		Superposition(const std::vector<Atom> &atoms1, const std::vector<Atom> &atoms2);
		void Clean(void);
		void Print(void);
		virtual ~Superposition() {}
		Superposition(const Atom *a1, const Atom *a2);
		float Rms(void);
	
		inline void Add(const Superposition *b)
		{
			na += b->na;	
			x1+=b->x1;y1+=b->y1;z1+=b->z1;
			x2+=b->x2;y2+=b->y2;z2+=b->z2;
			x1x1+=b->x1x1;y1y1+=b->y1y1;z1z1+=b->z1z1;
			x2x2+=b->x2x2;y2y2+=b->y2y2;z2z2+=b->z2z2;
			x1x2+=b->x1x2;x1y2+=b->x1y2;x1z2+=b->x1z2;
			y1x2+=b->y1x2;y1y2+=b->y1y2;y1z2+=b->y1z2;
			z1x2+=b->z1x2;z1y2+=b->z1y2;z1z2+=b->z1z2;
		}
		
		inline void Sub(const Superposition *b)
		{
			na -= b->na;	
			x1-=b->x1;y1-=b->y1;z1-=b->z1;
			x2-=b->x2;y2-=b->y2;z2-=b->z2;
			x1x1-=b->x1x1;y1y1-=b->y1y1;z1z1-=b->z1z1;
			x2x2-=b->x2x2;y2y2-=b->y2y2;z2z2-=b->z2z2;
			x1x2-=b->x1x2;x1y2-=b->x1y2;x1z2-=b->x1z2;
			y1x2-=b->y1x2;y1y2-=b->y1y2;y1z2-=b->y1z2;
			z1x2-=b->z1x2;z1y2-=b->z1y2;z1z2-=b->z1z2;
		}
};
	
#endif

