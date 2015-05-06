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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "Superposition.h"
#include <cassert>

using std::vector;

const double PI23 = 2*3.1415926535897932/3.0;
const double PI = 3.1415926535897932;
const double THIRD = 1./3.;

int signR(double Z)
{
	int ret = 0;
	if (Z > 0.0)	ret = 1;
	if (Z < 0.0)	ret = -1;
	if (Z == 0.0)	ret =0;

	return ret;
}

double CBRT(double Z)
{
	double ret = 0;
	const double THIRD = 1./3.;
	ret = fabs(pow(fabs(Z),THIRD)) * signR(Z);
	return ret;
}

double SQRS(double Z)
{
	double ret = 0;
	const double HALF = 1./2.;
	ret = fabs(pow(fabs(Z),HALF))* signR(Z);
	return ret;
}

float Superposition::Rms(void)
{
	double n = na;

	double xx=x1x2-x1*x2/n,xy=x1y2-x1*y2/n,xz=x1z2-x1*z2/n;
	double yx=y1x2-y1*x2/n,yy=y1y2-y1*y2/n,yz=y1z2-y1*z2/n; 
	double zx=z1x2-z1*x2/n,zy=z1y2-z1*y2/n,zz=z1z2-z1*z2/n; 
	double d1= (double)x1x1+(double)x2x2+(double)y1y1+(double)y2y2+(double)z1z1+(double)z2z2;
	double d2= ((double)x1*x1+(double)x2*x2+(double)y1*y1+(double)y2*y2+(double)z1*z1+(double)z2*z2)/n;	
	
	double XX=xx*xx,XY=xy*xy, XZ=xz*xz;
	double YX=yx*yx,YY=yy*yy, YZ=yz*yz;
	double ZX=zx*zx,ZY=zy*zy, ZZ=zz*zz;
	
	double t1=XX+XY+XZ+YX+YY+YZ+ZX+ZY+ZZ;
	double t2=xz*yy*zx -xy*yz*zx -xz*yx*zy +xx*yz*zy +xy*yx*zz -xx*yy*zz;
	double t3=XX*YY +XX*YZ +XX*ZY +XX*ZZ +XY*YX +XY*YZ +XY*ZX +XY*ZZ +XZ*YX +XZ*YY +XZ*ZX +XZ*ZY +YX*ZY +YX*ZZ +YY*ZX +YY*ZZ +YZ*ZX +YZ*ZY; 	
	double t4=xx*xy*yx*yy+ xx*xz*yx*yz+ xy*xz*yy*yz+ xx*xy*zx*zy+ yx*yy*zx*zy+ xx*xz*zx*zz+ yx*yz*zx*zz+ xy*xz*zy*zz+ yy*yz*zy*zz; 
	
	double best2 = 0;

	double T1 = -t1/3.0;
	double T5 = t3-2.0*t4;
	double P = pow((T5/3.0 -T1*T1),3);
	double Q = -T1*T1*T1 +T5*T1*.5 +t2*t2*.5;
	double DIS = pow(Q,2)+P;
	if ( DIS < 0.0 )
	{
		double PHI=acos(Q/sqrt(-P))/3.0;
		double T=2.0*pow((-P),1.0/6.0);
		
		double s0= SQRS(T*cos(PHI)-T1);
		if((t2 > -1e-8) && (t2 < 1e-8)) best2 = s0;
		else { 
			double s1= SQRS(T*cos(PHI-PI23)-T1);
			double s2= SQRS(T*cos(PHI+PI23)-T1);
			if (t2>0) best2=s0+s1-s2;
			if (t2<0) best2=s0+s1+s2;
		}	
	}
	else best2=SQRS(CBRT(Q+sqrt(DIS))+CBRT(Q-sqrt(DIS))-T1);

	double rms2=(d1-d2-2.0*best2)/ n;
	if((rms2 > -1e-6) && (rms2 < 1e-6)) rms2=0.0;
	
	return (float)sqrt(rms2);
}



Superposition::Superposition(void)
{	
	na=0;
	x1=0;y1=0;z1=0;
	x2=0;y2=0;z2=0;
	x1x1=0;y1y1=0;z1z1=0;
	x2x2=0;y2y2=0;z2z2=0;
	x1x2=0;x1y2=0;x1z2=0;
	y1x2=0;y1y2=0;y1z2=0;
	z1x2=0;z1y2=0;z1z2=0;
}

Superposition::Superposition(const vector<Atom> &atoms1, const vector<Atom> &atoms2)
{
	na=0;
	x1=0;y1=0;z1=0;
	x2=0;y2=0;z2=0;
	x1x1=0;y1y1=0;z1z1=0;
	x2x2=0;y2y2=0;z2z2=0;
	x1x2=0;x1y2=0;x1z2=0;
	y1x2=0;y1y2=0;y1z2=0;
	z1x2=0;z1y2=0;z1z2=0;
	
	assert(atoms1.size() == atoms2.size());
	for(unsigned int i = 0; i < atoms1.size(); i++)
	{
		Superposition temp(&(atoms1[i]), &(atoms2[i]));
		this->Add(&temp);
	}
}
	

Superposition::Superposition(const Atom *a1, const Atom *a2)
{
	na=1;	
	x1x2=a1->x*a2->x;x1y2=a1->x*a2->y;x1z2=a1->x*a2->z;
	y1x2=a1->y*a2->x;y1y2=a1->y*a2->y;y1z2=a1->y*a2->z;
	z1x2=a1->z*a2->x;z1y2=a1->z*a2->y;z1z2=a1->z*a2->z;
	x1x1=a1->x*a1->x;y1y1=a1->y*a1->y;z1z1=a1->z*a1->z;
	x2x2=a2->x*a2->x;y2y2=a2->y*a2->y;z2z2=a2->z*a2->z;
	x1=a1->x;y1=a1->y;z1=a1->z;
	x2=a2->x;y2=a2->y;z2=a2->z;
}

void Superposition::Print(void)
{
	printf("super(%d %ld %ld %ld %ld %ld %ld)\n",na,x1,x2,y1,y2,z1,z2); 
	printf("__________(%d %ld %ld %ld %ld %ld %ld)\n",na,x1x1,x2x2,y1y1,y2y2,z1z1,z2z2); 
	printf("__________(%d %ld %ld %ld %ld %ld %ld %ld %ld %ld)\n",na,x1x2,x1y2,x1z2,y1x2,y1y2,y1z2,z1x2,z1y2,z1z2); 
}

void Superposition::Clean(void)
{
	na=0;
	x1=0;y1=0;z1=0;
	x2=0;y2=0;z2=0;
	x1x1=0;y1y1=0;z1z1=0;
	x2x2=0;y2y2=0;z2z2=0;
	x1x2=0;x1y2=0;x1z2=0;
	y1x2=0;y1y2=0;y1z2=0;
	z1x2=0;z1y2=0;z1z2=0;
}



/*
 Superposition::Superposition(struct Amino *a1, struct Amino *a2)
 {
 na=1;	
 x1x2=a1->x*a2->x+a1->cx*a2->cx;x1y2=a1->x*a2->y+a1->cx*a2->cy;x1z2=a1->x*a2->z+a1->cx*a2->cz;
 y1x2=a1->y*a2->x+a1->cy*a2->cx;y1y2=a1->y*a2->y+a1->cy*a2->cy;y1z2=a1->y*a2->z+a1->cy*a2->cz;
 z1x2=a1->z*a2->x+a1->cz*a2->cx;z1y2=a1->z*a2->y+a1->cz*a2->cy;z1z2=a1->z*a2->z+a1->cz*a2->cz;
 x1x1=a1->x*a1->x+a1->cx*a1->cx;y1y1=a1->y*a1->y+a1->cy*a1->cy;z1z1=a1->z*a1->z+a1->cz*a1->cz;
 x2x2=a2->x*a2->x+a2->cx*a2->cx;y2y2=a2->y*a2->y+a2->cy*a2->cy;z2z2=a2->z*a2->z+a2->cz*a2->cz;
 x1=a1->x+a1->cx;y1=a1->y+a1->cy;z1=a1->z+a1->cz;
 x2=a2->x+a2->cx;y2=a2->y+a2->cy;z2=a2->z+a2->cz;
 }*/



/*
 
 
 Superposition::Superposition(struct Amino *a1, struct Amino *a2)
 {
 na=1;	
 x1x2=a1->x*a2->x;x1y2=a1->x*a2->y;x1z2=a1->x*a2->z;
 y1x2=a1->y*a2->x;y1y2=a1->y*a2->y;y1z2=a1->y*a2->z;
 z1x2=a1->z*a2->x;z1y2=a1->z*a2->y;z1z2=a1->z*a2->z;
 x1x1=a1->x*a1->x;y1y1=a1->y*a1->y;z1z1=a1->z*a1->z;
 x2x2=a2->x*a2->x;y2y2=a2->y*a2->y;z2z2=a2->z*a2->z;
 x1=a1->x;y1=a1->y;z1=a1->z;
 x2=a2->x;y2=a2->y;z2=a2->z;
 
 #ifndef SUPERPOSE_CA_ONLY	
 x1x2+=a1->cx*a2->cx;x1y2+=a1->cx*a2->cy;x1z2+=a1->cx*a2->cz;
 y1x2+=a1->cy*a2->cx;y1y2+=a1->cy*a2->cy;y1z2+=a1->cy*a2->cz;
 z1x2+=a1->cz*a2->cx;z1y2+=a1->cz*a2->cy;z1z2+=a1->cz*a2->cz;
 x1x1+=a1->cx*a1->cx;y1y1+=a1->cy*a1->cy;z1z1+=a1->cz*a1->cz;
 x2x2+=a2->cx*a2->cx;y2y2+=a2->cy*a2->cy;z2z2+=a2->cz*a2->cz;
 x1+=a1->cx;y1+=a1->cy;z1+=a1->cz;
 x2+=a2->cx;y2+=a2->cy;z2+=a2->cz;
 #endif
 }
 
 */

/*
 void Superposition::Init(struct Amino *a1, struct Amino *a2)
 {
 na=1;	
 x1x2=a1->x*a2->x+a1->cx*a2->cx;x1y2=a1->x*a2->y+a1->cx*a2->cy;x1z2=a1->x*a2->z+a1->cx*a2->cz;
 y1x2=a1->y*a2->x+a1->cy*a2->cx;y1y2=a1->y*a2->y+a1->cy*a2->cy;y1z2=a1->y*a2->z+a1->cy*a2->cz;
 z1x2=a1->z*a2->x+a1->cz*a2->cx;z1y2=a1->z*a2->y+a1->cz*a2->cy;z1z2=a1->z*a2->z+a1->cz*a2->cz;
 x1x1=a1->x*a1->x+a1->cx*a1->cx;y1y1=a1->y*a1->y+a1->cy*a1->cy;z1z1=a1->z*a1->z+a1->cz*a1->cz;
 x2x2=a2->x*a2->x+a2->cx*a2->cx;y2y2=a2->y*a2->y+a2->cy*a2->cy;z2z2=a2->z*a2->z+a2->cz*a2->cz;
 x1=a1->x+a1->cx;y1=a1->y+a1->cy;z1=a1->z+a1->cz;
 x2=a2->x+a2->cx;y2=a2->y+a2->cy;z2=a2->z+a2->cz;
 }*/


