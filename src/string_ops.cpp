/***************************************************************************
 *   Copyright (C) 2006 by Pier Federico Gherardini   *
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

#include <string>
#include <vector>
#include "common.h"
#include <limits.h>
#include <cstdlib>

using std::string;
using std::vector;
using std::list;
using std::pair;


void split(vector<string> *fields, const string &s, char delim)
{
	unsigned int i, cur;
	fields->clear();
	fields->resize(1);
	cur = 0;
	for(i = 0; i < s.size(); i++)
	{
		if(s[i] != delim)
			(*fields)[cur] += s[i];
		else
		{
			fields->resize(fields->size() + 1);
			cur++;
		}
	}
}	

/*	
vector<string> *split(const char *s, char delim)
{
    vector<string> *res = new vector<string>(1);
    unsigned int i, cur;
    cur = 0;
    for(i = 0; i < strlen(s); i++)
    {
        if(s[i] != delim)
            (*res)[cur] += s[i];
        else
        {
            res->resize(res->size() + 1);
            cur++;
        }
    }
    return res;
}
*/

string trim_copy(const string &_s)
{
	string::size_type start = _s.find_first_not_of(" \n\t\r");
	string::size_type end = _s.find_last_not_of(" \n\t\r");
	if(start == string::npos || end == string::npos) 
		return "";
	else
		return _s.substr(start, end - start + 1);
}

void trim(string *_s)
{
	string &s = *_s;
	string::size_type start = s.find_first_not_of(" \n\t\r");
	string::size_type end = s.find_last_not_of(" \n\t\r");
	if(start == string::npos || end == string::npos) 
		s.clear();
	else
		s = s.substr(start, end - start + 1);
}


void process_range(const string &s, list<pair<int, int> > *ranges)
{
	ranges->clear();
	vector<string> fields;
	split(&fields, s, ',');
	vector<string>::iterator i = fields.begin();
	for(; i != fields.end(); i++)
	{
		pair<int, int> p;
		vector<string> temp;
		split(&temp, (*i), '-');
		p.first = atoi(temp[0].c_str());
		if(temp.size() == 1)
			p.second = INT_MIN;
		else if(temp.size() == 2)
			p.second = atoi(temp[1].c_str());
		ranges->push_back(p);
	}
}


string base_name(const string &s, bool remove_suffix)
{
	//ATTENZIONE!! SU WINDOWS QUESTA FUNZIONE VA CAMBIATA!!!
	string::size_type first = s.find_last_of("/\\");
	if(first == string::npos)
		first = 0;
	else 
		first++;
	string::size_type last;
	if(remove_suffix)
	{
		last = s.find_last_of(".");
		if(last == string::npos)
			last = s.size();
	}
	else
		last = s.size();
	return s.substr(first, last - first);
}
	
	
string remove_spaces(const string &_s)
{
	string s;
	string::const_iterator x;
	for(x = _s.begin(); x != _s.end(); x++)
		if(!isspace(*x))
			s += *x;
	return s;
}

string remove_insertion_code(const string &_s)
{
	string::const_iterator x;
	string s;
    //Rimuoviamo insertion code
    for(x = _s.begin(); x != _s.end(); x++)
        if(isdigit(*x))
			s += *x;
	return s;
}

/*
bool in_range(const pair<int, int> p, const string &_s)
{
	string s = remove_insertion_code(_s);
	if(!s.size())
		return false;
	else
	{
		int val = atoi(s.c_str());
		if(p.second == -1) //It is a single residue
		{
			if(val == p.first)
				return true;
			else
				return false;
		}
		else
		{
			if(val >= p.first && val <= p.second)
				return true;
			else
				return false;
		}
	}
}
*/

bool in_range(const list<pair<int, int> > &ranges, const string &_s)
{
	string::const_iterator x;
	string s = remove_insertion_code(_s);
			
	if(!s.size())
		return false;
	else
	{
		int val = atoi(s.c_str());
		list<pair<int, int> >::const_iterator i;
		for(i = ranges.begin(); i != ranges.end(); i++)
		{
			bool in_this_range = false;
			if(val < i->first)
				in_this_range = false;
			else if(val == i->first)
				in_this_range = true;
			//If the upper limit is not defined
			//p.second == INT_MIN and therefore val will never be greater
			else if(val <= i->second)
				in_this_range = true;
			if(in_this_range)
				return true;
		}
	}
	return false;
}





