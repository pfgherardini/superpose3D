/***************************************************************************
 *   Copyright (C) 2006 - 2009 by Gabriele Ausiello, Pier Federico Gherardini *
 *   gabriele.ausiello@uniroma2.it 
 *	 pier.federico.gherardini@uniroma2.it								   *
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

#ifndef _SUPERPOSE_TREE_H
#define _SUPERPOSE_TREE_H
#include <set>
#include <cstring>

class Tree 
{
		class Code
		{
			public:
				short *_c;
				int size;
				Code(int _size) : size(_size)
				{
					_c = new short[size * 2];
				}
				Code(const Code &other) : size(other.size)
				{
					_c = new short[size * 2];
					memcpy(_c, other._c, size*2*sizeof(short));
				}
				~Code(void)
				{
					if(_c)
						delete [] _c;
				}
		};
		
		class Code_lt
		{
			public:
				inline bool operator()(const Code *a, const Code *b)
				{
					return memcmp(a->_c, b->_c, a->size*2*sizeof(short)) < 0;
				}
		};
					
		std::set<Code *, Code_lt> _tree;
		Code *cur_code;
		int size;
	
	public:
		Tree (void);
		const short *getCurrent_code(void) const {return cur_code->_c; }
		virtual ~Tree();
		void Init(int);
		int Add(short,short);
		void Remove(int to_remove);
		bool Done();
		void Clear();		
};


#endif



