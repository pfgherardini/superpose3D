/***************************************************************************
 *   Copyright (C) 2007 by Pier Federico Gherardini   *
 *   pier.federico.gherardini@uniroma2.it   *
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

#include "Quaternion.h"
#include <cmath>

const double PI = 3.14159;

Quaternion::Quaternion(float X, float Y, float Z, float angle)
{
	float vec[3] = {X, Y, Z};
	Normalize_vector(vec);
	float rad = float((angle / 180.0f) * PI);
	w = cos(rad / 2.0f);
	float temp = sin(rad / 2.0f);
	x = vec[0] * temp;
	y = vec[1] * temp;
	z = vec[2] * temp;
}

void Quaternion::Normalize_vector(float v[3])
{
	float norm = sqrt(Sqr_len(v));
	v[0] /= norm;
	v[1] /= norm;
	v[2] /= norm;
}

void Quaternion::Rotate(float c[3], float ca[3], float r[3])
{
	float M[3][3];
	Create_rotation_matrix(M);
	float v[3] = {c[0] - ca[0], c[1] - ca[1], c[2] - ca[2]};
	for(int i = 0; i < 3; i++)
	{
		r[i] = 0.0;
		for(int j = 0; j < 3; j++)
			r[i] += v[j] * M[j][i];
	}
	
	r[0] += ca[0];
	r[1] += ca[1];
	r[2] += ca[2];
}
		
void Quaternion::Create_rotation_matrix(float matrix[3][3])
{
	matrix[0][0] = 1 - 2 * square(y) - 2 * square(z);
	matrix[0][1] = 2 * x * y - 2 * w * z;
	matrix[0][2] = 2 * x * z + 2 * w * y;
	
	matrix[1][0] = 2 * x * y + 2 * w * z;
	matrix[1][1] = 1 - 2 * square(x) - 2 * square(z);
	matrix[1][2] = 2 * y * z - 2 * w * x;
	
	matrix[2][0] = 2 * x * z - 2 * w * y;
	matrix[2][1] = 2 * y * z + 2 * w * x;
	matrix[2][2] = 1 - 2 * square(x) - 2 * square(y);
}



