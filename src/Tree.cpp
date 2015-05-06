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

#include "stdlib.h"
#include "string.h"
#include <cstdio> 
#include "Tree.h"


using std::set;
using std::pair;

Tree::~Tree() {
	if(cur_code)
		delete cur_code;
}

Tree::Tree() : cur_code(0) {}

void Tree::Init(int n)//This should go into the constructor
{
	size=n;
	cur_code = new Code(size);
	for(int i = 0; i < cur_code->size * 2; i++)
		cur_code->_c[i] = 0;
}

bool Tree::Done()
{
	Code *val = new Code(*cur_code);
	pair<set<Code *, Code_lt>::iterator, bool> p = _tree.insert(val);
	if(p.second) //It wasn't already there
		return false;
	else
	{
		delete val;
		return true;
	}
}

void Tree::Clear()//Dovrebbe stare nel destructor
{
	set<Code *, Code_lt>::iterator i;
	for(i = _tree.begin(); i != _tree.end(); i++)
		delete *i;
	_tree.clear();
}





int Tree::Add(short a1, short a2)
{
	int i;
	for(i = 0; i < size - 1; i++)
		if(cur_code->_c[i * 2] < a1)
			break;
	if(i < size - 1)//i.e. this is not the last position
		memmove(cur_code->_c + (i * 2) + 2, cur_code->_c + (i * 2), (size * 2 - (i * 2) - 2) * sizeof(short));
	cur_code->_c[i * 2] = a1;
	cur_code->_c[(i * 2) + 1] = a2;
	/*printf("Add    "); //DEBUG
	for(int a = 0; a < size; a++) printf("%hd%hd", code[a * 2], code[a*2 + 1]);//DEBUG
	printf("\n"); //DEBUG*/
	return i;
}	

void Tree::Remove(int to_remove)
{
	int i = to_remove;
	if(i < size - 1)
		memmove(cur_code->_c + (i * 2), cur_code->_c + (i * 2) + 2, (size * 2 - (i * 2) - 2) * sizeof(short));
	cur_code->_c[(size-1)*2]=0;
	cur_code->_c[(size-1)*2+1]=0;
	/*printf("Remove "); //DEBUG
	for(int a = 0; a < size; a++) printf("%hd%hd", code[a * 2], code[a*2 + 1]);//DEBUG
	printf("\n"); //DEBUG*/
}



