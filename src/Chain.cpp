/***************************************************************************
 *   Copyright (C) 2006 by Gabriele Ausiello, Pier Federico Gherardini   *
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

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <stddef.h>
#include <iostream>
#include "common.h"
#include "Chain.h"

using std::cout;
using std::string;
using std::vector;
using std::ofstream;

Chain::Chain(const Chain &ch, const vector<bool> &sel) : file_name(ch.file_name), 
				chain_id(ch.chain_id), offsets(ch.offsets)
{
	vector<bool>::const_iterator x;
	unsigned int n = 0;
	for(x = sel.begin(); x != sel.end(); x++)
		if(*x)
			n++;
	residues.resize(n);
	unsigned int i, j;
	for(i = j = 0; i < ch.size(); i++)
		if(sel[i])
		{			
			residues[j] = new Amino(ch[i]);
			j++;
		}
}               
/*
Chain *Chain::Append(const Chain &a)
{
	Chain *ret = new Chain(this->file_name, this->getChain_id());
	const vector<float> &off_this = this->getOffsets();
	const vector<float> &off_a = a->getOffsets();
	vector<float> off_new(3);
	for(unsigned int i = 0; i < off_new.size(); i++)
		off_new[i] = off_this[i] - off_a[i];
	ret->setOffsets(off_new);
*/	


int Chain::getIndex(const std::string &resnum) const
{
	for(unsigned int i = 0; i < residues.size(); i++)
		if((residues[i])->getPdb_resnum() == resnum)
			return i;
	return -1;
}



void Chain::Write(const string &file_name)
{
	ofstream out(file_name.c_str());
	vector<Amino *>::const_iterator i;
	
	out << file_name << "\t" << residues.size();
	for(unsigned int k = 0; k < offsets.size(); k++)
		out << "\t" << offsets[k];
	out << "\n";
	
	for(i = residues.begin(); i != residues.end(); i++)
	{
		const Amino &a = **i;
		out << a.getPdb_resnum() << "\t" << a.getLabel();
		for(unsigned int j = 0; j < a.getNum_atoms(); j++)
			if(a[j])
				out << "\t" << (a[j])->x << "\t" << (a[j])->y << "\t" << (a[j])->z;
		out << "\n";
	}
	out.close();
}




Chain::~Chain()
{
	vector<Amino *>::iterator i;
	
	for(i = residues.begin(); i != residues.end(); i++)
	{
		assert(*i);
		delete *i;
	}
}



