//Supplementary functions for the healpix algorithm
#include <stdio.h>
#include <omp.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "general_funcs.cc"
using namespace std;

struct rotation_type  {
	float theta;
	vec3 basis_set[3];
};

//random angles
int const Nalpha = 35; //don't change this EVER
int const Ntheta = pow(Nalpha, 3);
float random_angles[Ntheta][3];

//create list of Ntheta = Nalpha^3 angles on a grid from [0, 2*pi)
void get_random_angles(void)  {
	#pragma omp parallel
	{
	#pragma omp for
	for (int i = 0; i < Nalpha; i++)  {
		for (int j = 0; j < Nalpha; j++)  {
			for (int k = 0; k < Nalpha; k++)  {
				int ind = i*pow(Nalpha, 2) + j*Nalpha + k;
				float alpha0 = 2.*pi/Nalpha*i;
				float beta0  = 2.*pi/Nalpha*j;
				float gamma0 = 2.*pi/Nalpha*k;
				random_angles[ind][0] = alpha0;
				random_angles[ind][1] = beta0;
				random_angles[ind][2] = gamma0;
			}
		}
	}
	}
}

//Get rotation angle/unit vectors
struct rotation_type get_rotation(vec3 v1, vec3 v2)  {
	struct rotation_type rotation;
	
	//angle between v1 and v2
	float theta = acos(dot(v1, v2));
	
	//get the unit vector ortogonal to the plane of v1, v2
	vec3 u3 = cross(v2, v1);
	u3 = u3/length(u3);
	
	//find a unique basis
	vec3 u1 = v1/length(v1);
	vec3 u2 = cross(u1, u3);
	u2 = u2/length(u2);
	
	//if in the lower half plane, correct theta
	if (dot(v2, u2) < 0.)  {
		theta = 2.*pi - theta;
	}
	
	rotation.theta = theta;
	rotation.basis_set[0] = u1;
	rotation.basis_set[1] = u2;
	rotation.basis_set[2] = u3;
	
	return rotation;
}

//Solve for the rotation operation that maps v1 onto v2, then apply it to v3 to get v4. 
//This is so I don't have to store the random rotation angles for each unit vector.  
//Acknowlegements: Robert Dawson

vec3 rotate(struct rotation_type rotation, vec3 v3)  {
	vec3 v4;
	
	float theta = rotation.theta;
	vec3 u1 = rotation.basis_set[0];
	vec3 u2 = rotation.basis_set[1];
	vec3 u3 = rotation.basis_set[2];
	
	//find the component of v3 in the new basis
	float v3u1 = dot(v3, u1);
	float v3u2 = dot(v3, u2);
	float v4u3 = dot(v3, u3);
	
	//rotate around u3
	float v4u1 = cosf(theta)*v3u1 - sinf(theta)*v3u2;
	float v4u2 = sinf(theta)*v3u1 + cosf(theta)*v3u2;
	
	//return to cartesian basis
	v4.x = v4u1*u1.x + v4u2*u2.x + v4u3*u3.x;
	v4.y = v4u1*u1.y + v4u2*u2.y + v4u3*u3.y;
	v4.z = v4u1*u1.z + v4u2*u2.z + v4u3*u3.z;
	
	return v4;
}

//Rotate a unit vector around the x, y, and z axes through a set of randonly chosen angles.  
//The direction tag indicates the order of multiplication of the rotation matrices.  
//direction = 0 - 5 are the forward rotations, 6-11 are the inverse rotations.  
vec3 get_unit_vector(vec3 v, float alpha0, float beta0, float gamma0, unsigned short int direction)  {
	vec3 vec_out;
	float xhat, yhat, zhat;
	float xhatp, yhatp, zhatp;
	
	xhat = (float) v.x;
	yhat = (float) v.y;
	zhat = (float) v.z;
	
	switch (direction)  {
		case 0 : 
			xhatp = (cos(beta0)*cos(gamma0))*xhat 
				  + (-cos(beta0)*sin(gamma0))*yhat
				  + (sin(beta0))*zhat;
			yhatp = (sin(alpha0)*sin(beta0)*cos(gamma0) + cos(alpha0)*sin(gamma0))*xhat 
				  + (-sin(alpha0)*sin(beta0)*sin(gamma0) + cos(alpha0)*cos(gamma0))*yhat
				  + (-sin(alpha0)*cos(beta0))*zhat;
			zhatp = (-cos(alpha0)*sin(beta0)*cos(gamma0) + sin(alpha0)*sin(gamma0))*xhat
				  + (cos(alpha0)*sin(beta0)*sin(gamma0) + sin(alpha0)*cos(gamma0))*yhat
				  + (cos(alpha0)*cos(beta0))*zhat;
			break;
		case 6 :
			xhatp = (cos(gamma0)*cos(beta0))*xhat
				  + (sin(gamma0)*cos(alpha0) + cos(gamma0)*sin(beta0)*sin(alpha0))*yhat
				  + (sin(gamma0)*sin(alpha0) - cos(gamma0)*sin(beta0)*cos(alpha0))*zhat;
			yhatp = (-sin(gamma0)*cos(beta0))*xhat
				  + (cos(gamma0)*cos(alpha0) - sin(gamma0)*sin(beta0)*sin(alpha0))*yhat
				  + (cos(gamma0)*sin(alpha0) + sin(gamma0)*sin(beta0)*cos(alpha0))*zhat;
			zhatp = (sin(beta0))*xhat
				  + (-cos(beta0)*sin(alpha0))*yhat
				  + (cos(beta0)*cos(alpha0))*zhat;
			break;
		default : 
			printf("Invalid direction\n");
	}
	
	vec_out.x = xhatp;
	vec_out.y = yhatp;
	vec_out.z = zhatp;
	
	return vec_out; 
}