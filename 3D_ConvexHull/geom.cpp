#include "geom.h"
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstdlib>

#include <vector>

using namespace std;

double d = 0;



/* returns 2 times the signed area of triangle abc. The area is positive if c
 is to the left of ab, and negative if c is to the right of ab
 */
int signed_area3D(point3d a, point3d b, point3d c, point3d d) {
	int det = determinant(a,b,c,d);
	return det;
}

//calcualte the determinant of a 4X4 matrix based on the four points
//more info on this on page 26 of the textbook
int determinant(point3d a, point3d b, point3d c, point3d d) {
	return -(a.z-d.z)*(b.y-d.y)*(c.x-d.x)+(a.y-d.y)*(b.z-d.z)*(c.x-d.x)
			  +(a.z-d.z)*(b.x-d.x)*(c.y-d.y)-(a.x-d.x)*(b.z-d.z)*(c.y-d.y)
			  -(a.y-d.y)*(b.x-d.x)*(c.z-d.z)+(a.x-d.x)*(b.y-d.y)*(c.z-d.z);
}

/* return 1 if p,q,r, t on same plane, and 0 otherwise */
int coplanar(point3d p, point3d q, point3d r, point3d t) {
	return (signed_area3D(p, q, r, t) == 0);
}


/* return 1 if d is  strictly left of abc; 0 otherwise */
int left (point3d a, point3d b, point3d c, point3d d) {
	//note here to deal with co-planar points we added in the
	//= part of the >=
	return (signed_area3D(a, b, c, d) >= 0);
}


/* compute and return the convex hull of the points */
vector<triangle3d> brute_force_hull(vector<point3d> points) {
	
	//signed_area3D(points[0], points[0], points[0], points[0]);
	
	vector<triangle3d> result;
	
	//loop through all possible triplets of points
	for (int a = 0; a < points.size(); a++) {
		for (int b = 0; b < points.size(); b++) {
			//make sure not the same point
			if (a != b) {
				for (int c = 0; c < points.size(); c++) {
					//make sure not the same point
					if (a != c && b != c ) {
						//set variable equal to true, to be set to false
						//later on if the plane is not extreme
						bool allLeftOf = true;
						//loop through all other possible points to make a
						//tetrahedron with a,b,c
						for (int d = 0; d < points.size(); d++) {
							//make sure d isn't any of the other points
							if (a != d && b != d && c != d) {
								//if at any point a point d is right of a,b,c
								//break the loop and set "allLeftOf" to false
								if (!left(points.at(a), points.at(b), points.at(c), points.at(d))) {
									allLeftOf = false;
									break;
								}
							}
						}
						//if any points were right of don't add the three
						//else add the three points a,b,c to result as a triangle
						if (allLeftOf) {
							triangle3d triangle;
							triangle.a = points.at(a);
							triangle.b = points.at(b);
							triangle.c = points.at(c);
							result.push_back(triangle);
						}
					}
				}
			}
		}
	}
	
	return result;
}

