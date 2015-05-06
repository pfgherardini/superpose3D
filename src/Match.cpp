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

#include "Match.h"
#include <cstdlib>
#include <cstdio>
#include "superpose.h"
#include <map>

using std::pair;
using std::vector;
using std::string;

	
Match::Match(const ChainWrapper &c1, const ChainWrapper &c2, const options &_opt) : C1(c1), C2(c2), opt(_opt),
				rot_matrix(9), probe_trasl(3), target_trasl(3)
{
	rmsd = 0.0f;
}



void Match::Write(unsigned long id)
{
	FILE *f = opt.output_file;
	FILE *trans_f = opt.trans_file;
	const Chain &ch1 = C1.getParent_chain();
	const Chain &ch2 = C2.getParent_chain();
	
	if(opt.print_match_id)
		fprintf(f, "%ld\t", id);
	fprintf(f, "%s\t%s\t", ch1.getFull_id().c_str(), ch2.getFull_id().c_str());
	fprintf(f, "%d\t%d\t", C1.getEffective_size(), C2.getEffective_size());
	fprintf(f, "%.3f\t", rmsd);
	fprintf(f, "%d\t", this->size());
	
	for(unsigned int i = 0; i < this->size(); i++)
	{
		const PseudoResidue &pseudo1 = C1[amino1[i]];
		const PseudoResidue &pseudo2 = C2[amino2[i]];
		const Amino &res1 = ch1[pseudo1.idx_orig_residue];
		const Amino &res2 = ch2[pseudo2.idx_orig_residue];
		string lab1, lab2;
		if(opt.mode != options::mode_res)
		{	
			lab1 = pseudo1.orig_atom_name;    //opt.sim_mat->Idx2AAlabel(pseudo1.getIdx_in_matrix());
			lab2 = pseudo2.orig_atom_name;	  //opt.sim_mat->Idx2AAlabel(pseudo2.getIdx_in_matrix());
			lab1 = lab1.substr(lab1.find_first_of("/") + 1);
			lab2 = lab2.substr(lab2.find_first_of("/") + 1);
			lab1 = lab1 == "ALL_ATOMS" ? "" : lab1;
			lab2 = lab2 == "ALL_ATOMS" ? "" : lab2;
		}
		fprintf(f, "%s", res1.getLabel().c_str());
		if(lab1.size()) fprintf(f, ".%s", lab1.c_str());
		fprintf(f, "(%s);%s", res1.getPdb_resnum().c_str(), res2.getLabel().c_str());
		if(lab2.size()) fprintf(f, ".%s", lab2.c_str());
		fprintf(f, "(%s)", res2.getPdb_resnum().c_str());
		if(i != this->size() - 1)
			fprintf(f, "\t");
	}
	if(opt.print_transformations)
	{
		/*if(f == trans_f)
			fprintf(f, "\t");*/
	
		for(unsigned int i = 0; i < 3; i++)
			fprintf(trans_f, "\t%8.3f", probe_trasl[i] / 100 - ch1.getOffsets()[i]);

		for(unsigned int i = 0; i < 3; i++)
			fprintf(trans_f, "\t%8.3f", target_trasl[i] / 100 - ch2.getOffsets()[i]);
	
		for(unsigned int i = 0; i < 9; i++)
			fprintf(trans_f, "\t%8.3f", rot_matrix[i]);
	}
	if(f != trans_f)
	{
		fprintf(trans_f, "\n");
		fprintf(f, "\n");
	}
	else
		fprintf(f, "\n");
		
	fflush(f);
	fflush(trans_f);
}
	
/*
void Match::Write_set(const options &opt, const Set &s, const Chain &ch, short *amino)
{

	FILE *f = opt.output_file;
	for(int i = 0; i < n; i++)
	{
		fprintf(f, "%c\t%s\t", ch[amino[i]].one_letter_code,
										ch[amino[i]].pdb_resnum.c_str());
		//We write the C-alpha first and then the bariatom
		fprintf(f, "%d\t%d\t%d\t", (int)(s.p[i + n].x), (int)(s.p[i + n].y), (int)(s.p[i + n].z));
		fprintf(f, "%d\t%d\t%d\n", (int)(s.p[i].x), (int)(s.p[i].y), (int)(s.p[i].z));
	}
}
*/




void Match::Superpose(void)
{
	Set set1, set2;
	for(unsigned int i = 0; i < amino1.size(); i++)
	{
		unsigned int n1 = amino1[i];
		unsigned int n2 = amino2[i];
		
		unsigned short idx1 = C1[n1].getIdx_in_matrix();
		unsigned short idx2 = C2[n2].getIdx_in_matrix();
		pair<short, short> p = (*(opt.sim_mat))[idx1][idx2];
		const vector<Atom> &atoms1 = C1[n1][p.first];
		const vector<Atom> &atoms2 = C2[n2][p.second];
		set1.Add_atoms(atoms1);
		set2.Add_atoms(atoms2);
	}
	//Get the translation for probe and target;
	set1.Bring_to_zero(&probe_trasl);
	set2.Bring_to_zero(&target_trasl);	
	set_superpose(&set1, &set2, &rot_matrix);
	set2.Rotate(rot_matrix);
	rmsd = (float)(set1.Rmsd(set2) / 100.0);
}




/*
void Match::Superpose(Set &set1, Set &set2)
{
	
	struct Traslazione t1,t2;	
	set_init(&set1,n*2);set_init(&set2,n*2);

	float x,y,z;
	int k;
	for (k=0;k<n;k++) {
		x=S1->chain->amino[amino1[k]].x;
		y=S1->chain->amino[amino1[k]].y;
		z=S1->chain->amino[amino1[k]].z;
		set1.p[k+n].x=(double)x;set1.p[k+n].y=(double)y;set1.p[k+n].z=(double)z;
		
		x=S1->chain->amino[amino1[k]].cx;
		y=S1->chain->amino[amino1[k]].cy;
		z=S1->chain->amino[amino1[k]].cz;
		set1.p[k].x=(double)x;set1.p[k].y=(double)y;set1.p[k].z=(double)z;

		x=S2->chain->amino[amino2[k]].x;
		y=S2->chain->amino[amino2[k]].y;
		z=S2->chain->amino[amino2[k]].z;
		set2.p[k+n].x=(double)x;set2.p[k+n].y=(double)y;set2.p[k+n].z=(double)z;

		x=S2->chain->amino[amino2[k]].cx;
		y=S2->chain->amino[amino2[k]].cy;
		z=S2->chain->amino[amino2[k]].cz;
		set2.p[k].x=(double)x,set2.p[k].y=(double)y,set2.p[k].z=(double)z;
	}
	
	t1=set_zero(&set1);t2=set_zero(&set2);	
	set_trasla(&set1,&t1);set_trasla(&set2,&t2);
	set_superpose(&set1,&set2,r);
	tx1=t1.x;ty1=t1.y;tz1=t1.z;tx2=t2.x;ty2=t2.y;tz2=t2.z;
	set_ruota(&set2,r);//Questo direi proprio che serve...
	rms = (float)set_rms(&set1, &set2) / 100;
	//set_delete(&set1);
	//set_delete(&set2);
}
 */


