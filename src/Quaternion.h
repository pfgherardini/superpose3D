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
#ifndef _QUATERNION_H_
#define _QUATERNION_H_

class Quaternion
{
	private:
		float w, x, y, z;
		float square(float a) {return a * a;}
		void Normalize_vector(float v[3]);
		float Sqr_len(float v[3]) {return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];}
	public:
		Quaternion(void) {};
		Quaternion(float X, float Y, float Z, float angle);
		void Create_rotation_matrix(float matrix[3][3]);
		void Rotate(float c[3], float ca[3], float r[3]);
};
		
#endif //_QUATERNION_H_

