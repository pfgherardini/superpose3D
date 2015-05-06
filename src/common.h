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

#ifndef _SUPERPOSE_COMMON_H_
#define _SUPERPOSE_COMMON_H_
#include <vector>
#include "DefinitionCenter.h"
#include "SimilarityMatrix.h"
#include "PdbResidueDictionary.h"
#include <string>
#include <cstdio>
#include <map>
#include <list>


class SimilarityMatrix;
class DefinitionCenter;

struct options
{
	FILE *output_file, *trans_file;
	DefinitionCenter *def_center;
	SimilarityMatrix *sim_mat;
	PdbResidueDictionary *pdb_res_dict;
	std::string pdb_dir, param_file_name, input_file_name;
	
	float rmsd_thresh, neighbour_threshold;
	unsigned int score_min, score_max, max_longest_matches;
	bool verbose, write_coordinates, print_transformations,
			single_match, print_match_id, use_old_transformations;
	enum {mode_res, mode_res_multiple, mode_atom} mode;
	options(void)
	{
		pdb_dir = "./"; //ATTENZIONE questo andrà cambiato su Windows!!
		mode = mode_res;
		trans_file = output_file = 0;
		sim_mat = 0;
		pdb_res_dict = 0;
		def_center = 0;
		rmsd_thresh = 0.7f;
		neighbour_threshold = 7.5f;
		max_longest_matches = 10;
		score_max = 10;
		score_min = 3;
		verbose = false;
		write_coordinates = false;
		single_match = false;
		print_transformations = true;
		print_match_id = true;
		use_old_transformations = false;
	}
	~options(void)
	{}
};

std::string base_name(const std::string &s, bool remove_suffix = true);
void trim(std::string *_s);
std::string trim_copy(const std::string &_s);
//std::vector<std::string> *split(const char *s, char delim);
void split(std::vector<std::string> *fields, const std::string &s, char delim);
void process_range(const std::string &s, std::list<std::pair<int, int> > *ranges);
bool in_range(const std::list<std::pair<int, int> > &ranges, const std::string &_s);
std::string remove_spaces(const std::string &_s);



#endif //_COMMON_H_
