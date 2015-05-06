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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "superpose.h"


using std::vector;


//ATTENZIONE!!! PROBE E TARGET SONO SCAMBIATI????

void set_superpose(struct Set *target, struct Set *probe, vector<double> *r)
{
	double qmatr[5][5];
	double dd[5];
	double vv[5][5];
	unsigned int j;
	int i,indro;
	double weight = 1.0/(double)target->size();
	double xrj,yrj,zrj;
	double xpj,ypj,zpj;
	double q1,q2,q3,q4;

	
	assert(target->size() == probe->size());
	
	for (i=1;i<=4;i++) for (j=1;j<=4;j++) qmatr[i][j]=0.0;
	
	for (j=0;j<target->size();j++)
	{
		
		xpj=target->p[j].x;ypj=target->p[j].y;zpj=target->p[j].z;
		xrj=probe->p[j].x; yrj=probe->p[j].y; zrj=probe->p[j].z;

		qmatr[1][1]+=weight*(xpj*xpj+ypj*ypj+zpj*zpj+
			xrj*xrj+yrj*yrj+zrj*zrj
			-2.0*xpj*xrj-2.0*ypj*yrj-2.0*zpj*zrj);
		qmatr[2][1]+=weight*2.0*(ypj*zrj-zpj*yrj);
		qmatr[3][1]+=weight*2.0*(-xpj*zrj+zpj*xrj);
		qmatr[4][1]+=weight*2.0*(xpj*yrj-ypj*xrj);
		qmatr[2][2]+=weight*(xpj*xpj + ypj*ypj + zpj*zpj +
			xrj*xrj + yrj*yrj + zrj*zrj -
			2.0*xpj*xrj + 2.0*ypj*yrj + 2.0*zpj*zrj);
		qmatr[3][2]+=weight*(-2.0*(xpj*yrj+ypj*xrj));
		qmatr[4][2]+=weight*(-2.0*(xpj*zrj+zpj*xrj));
		qmatr[3][3]+=weight*(xpj*xpj + ypj*ypj + zpj*zpj +
			xrj*xrj + yrj*yrj + zrj*zrj +
			2.0*xpj*xrj - 2.0*ypj*yrj + 2.0*zpj*zrj);
		qmatr[4][3]+=weight*(-2.0*(ypj*zrj+zpj*yrj));
		qmatr[4][4]+=weight*(xpj*xpj + ypj*ypj + zpj*zpj +
			xrj*xrj + yrj*yrj + zrj*zrj +
			2.0*xpj*xrj + 2.0*ypj*yrj - 2.0*zpj*zrj);
        }
	qmatr[1][2]=qmatr[2][1];qmatr[1][3]=qmatr[3][1];qmatr[1][4]=qmatr[4][1];
	qmatr[2][3]=qmatr[3][2];qmatr[2][4]=qmatr[4][2];qmatr[3][4]=qmatr[4][3];

	jacobi(qmatr,4,dd,vv);
	

	if (dd[1]<dd[2]) indro=1 ; else indro=2;
	if( dd[3]<dd[indro]) indro=3;
	if( dd[4]<dd[indro]) indro=4;

	q1=vv[indro][1];q2=vv[indro][2];q3=vv[indro][3];q4=vv[indro][4];
	
	r->resize(9);
	(*r)[0]= q1*q1 + q2*q2 - q3*q3 - q4*q4;
	(*r)[1]= 2.0*( - q1*q4 + q2*q3  );
	(*r)[2]= 2.0*(   q1*q3 + q2*q4  );
	(*r)[3]= 2.0*(   q1*q4 + q2*q3  );
	(*r)[4]= q1*q1 - q2*q2 + q3*q3 - q4*q4;
	(*r)[5]= 2.0*( - q1*q2 + q3*q4 );
	(*r)[6]= 2.0*( - q1*q3 + q2*q4 );
	(*r)[7]= 2.0*(   q1*q2 + q3*q4 );
	(*r)[8]= q1*q1 - q2*q2 - q3*q3 + q4*q4;
}

void jacobi (double a[5][5], int n, double  *d, double vv[5][5])
{
    double b[100],z[100];
	int ip,iq,i,j;
	double tresh,c,s,h,g,tau,sm,t,theta;
	int nrot;

	for (ip=1;ip<=n;ip++)
	{
		for (iq=1;iq<=n;iq++) vv[iq][ip]=0.0;	     			
			vv[ip][ip]=1.0;
	}

	for (ip=1;ip<=n;ip++)
	{
        	b[ip]=a[ip][ip];
      	  	d[ip]=b[ip];
        	z[ip]=0.0;
	}

	nrot=0;
	for (i=1;i<=50;i++)
	{
		sm=0.0;
		for (ip=1;ip<=n-1;ip++)
			for (iq=ip+1;iq<=n;iq++) sm+=fabs(a[iq][ip]);
		if(sm==0.0) return;
		if(i<4) 
			tresh=(0.2*sm)/(n*n);
		else 
			tresh=0.0;
		for(ip=1;ip<=n-1;ip++)
		{
			for (iq=ip+1;iq<=n;iq++)
			{
				g=100.0*fabs(a[iq][ip]);
				if ((i>4)&&(fabs(d[ip])+g==fabs(d[ip]))
					&& (fabs(d[iq])+g==fabs(d[iq]))) a[iq][ip]=0.0;
				else if(fabs(a[iq][ip])>tresh)
				{
					h=d[iq]-d[ip];
					if(fabs(h)+g==fabs(h)) t=a[iq][ip]/h;
					else
					{
						theta=0.5*h/a[iq][ip];
						t=1./(fabs(theta)+sqrt(1.0+theta*theta));
						if(theta<0.0)t=-t;
					}
              
					c=1./sqrt(1+t*t);
					s=t*c;
					tau=s/(1.+c);
					h=t*a[iq][ip];
					z[ip]=z[ip]-h;
					z[iq]=z[iq]+h;
					d[ip]=d[ip]-h;
					d[iq]=d[iq]+h;
					a[iq][ip]=0.0;
					
					for (j=1;j<=ip-1;j++)
					{
						g=a[ip][j];
						h=a[iq][j];
						a[ip][j]=g-s*(h+g*tau);
						a[iq][j]=h+s*(g-h*tau);
					}
					
					for(j=ip+1;j<=iq-1;j++)
					{
						g=a[j][ip];
						h=a[iq][j];
						a[j][ip]=g-s*(h+g*tau);
						a[iq][j]=h+s*(g-h*tau);
					}
					
					for (j=iq+1;j<=n;j++)
					{
						g=a[j][ip];
						h=a[j][iq];
						a[j][ip]=g-s*(h+g*tau);
						a[j][iq]=h+s*(g-h*tau);
					}
					for (j=1;j<=n;j++)
					{
						g=vv[ip][j];
						h=vv[iq][j];
						vv[ip][j]=g-s*(h+g*tau);
						vv[iq][j]=h+s*(g-h*tau);
					}
					nrot++;
				}
			}
		}
		
		for (ip=1;ip<=n;ip++)
		{
          		b[ip]=b[ip]+z[ip];
          		d[ip]=b[ip];
          		z[ip]=0.;
		}
	}
	printf("50 iterations should never happen\n");
}



void Set::Add_atoms(const vector<Atom> &v)
{
	for(unsigned int i = 0; i < v.size(); i++)
	{
		p.resize(p.size() + 1);
		Point &point = p.back();
		point.x = v[i].x;
		point.y = v[i].y;
		point.z = v[i].z;
	}
}
	

//Determina l'offset del centro di massa dallo zero
//sum(x)/N - t = 0
//quindi t = -sum(x)/N.
//D'altra parte w = 1 / N, quindi sum(x)/N = sum(x * w);

void Set::Bring_to_zero(vector<double> *trasl)
{
	trasl->resize(3);
	trasl->assign(trasl->size(), 0.0);
	double w = (double)(1.0) / p.size();
	
	for(unsigned int i = 0; i < p.size(); i++)
	{
		(*trasl)[0] -= p[i].x * w;
		(*trasl)[1] -= p[i].y * w;
		(*trasl)[2] -= p[i].z * w;
	}
	
	for(unsigned int i = 0; i < p.size(); i++)
	{
		p[i].x += (*trasl)[0];
		p[i].y += (*trasl)[1];
		p[i].z += (*trasl)[2];
	}
}

void Set::Rotate(const vector<double> &r)
{
	assert(r.size() == 9);
	double xt,yt,zt;
	
	for(unsigned int i = 0; i < this->size(); i++)
	{	
		xt = r[0]*p[i].x + r[1]*p[i].y + r[2]*p[i].z;
		yt = r[3]*p[i].x + r[4]*p[i].y + r[5]*p[i].z;
		zt = r[6]*p[i].x + r[7]*p[i].y + r[8]*p[i].z;
		p[i].x = xt; p[i].y = yt; p[i].z = zt;
	}
}

double Set::Rmsd(const Set &s2)
{	
	double rmsd = 0.0;
	assert(this->size() == s2.size());
	
	for(unsigned int i = 0; i < this->size(); i++)
		rmsd += this->p[i].sqrDistance_from(s2.p[i]);
	return sqrt(rmsd / (double)(this->size()));
}

double Point::sqrDistance_from(const Point &p2)
{
	return (double)((this->x-p2.x)*(this->x-p2.x)+(this->y-p2.y)*(this->y-p2.y)+(this->z-p2.z)*(this->z-p2.z));	
}



