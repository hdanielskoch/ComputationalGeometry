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

/* **************************************** */
/* returns the signed area of triangle abc. The area is positive if c
   is to the left of ab, and negative if c is to the right of ab
 */
int signed_area2D(point2D a, point2D b, point2D c) {
	return ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}

/* **************************************** */
/* return 1 if p,q,r collinear, and 0 otherwise */
int collinear(point2D p, point2D q, point2D r) {
 	return (signed_area2D(p, q, r) == 0);
}

/* **************************************** */
/* return 1 if c is  strictly left of ab; 0 otherwise */
int left (point2D a, point2D b, point2D c) {
  	return (signed_area2D(a, b, c) > 0);
}

/* **************************************** */
/* return 1 if same point*/
int samePoint(point2D a, point2D b){
	return ((a.x == b.x) && (a.y == b.y));
}

/* **************************************** */
/* Determines if two segments intersect */
bool intersect(lineSegment2D seg1, lineSegment2D seg2) {
    bool a1left = left(seg2.p1, seg2.p2, seg1.p1);
    bool a2left = left(seg2.p1, seg2.p2, seg1.p2);
    bool b1left = left(seg1.p1, seg1.p2, seg2.p1);
    bool b2left = left(seg1.p1, seg1.p2, seg2.p2);

    return ((a1left != a2left) && (b1left != b2left));
}
