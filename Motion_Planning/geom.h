/*
   Laura Toma
*/


#ifndef __geom_h
#define __geom_h


typedef struct _point2d {
    int x,y;
    int index;
    double weight;
    int prev;
} point2D;


typedef struct _lineSegment2D {
    point2D p1, p2;
} lineSegment2D;


typedef struct _rect2D  {
    point2D origin;
    float width, height;
} rect2D;


//add any functions you might need to operate on these basic types
/* returns the signed area of triangle abc. The area is positive if c
   is to the left of ab, and negative if c is to the right of ab
 */
int signed_area2D(point2D a, point2D b, point2D c);


/* return 1 if p,q,r collinear, and 0 otherwise */
int collinear(point2D p, point2D q, point2D r);


/* return 1 if c is  strictly left of ab; 0 otherwise */
int left (point2D a, point2D b, point2D c);

/* Returns 1 if a and b are the same points */
int samePoint(point2D a, point2D b);

/* Determines if two segments intersect */
bool intersect(lineSegment2D seg1, lineSegment2D seg2);


/* Determines if a ray and segment properly intersect (endpoint of segment on ray) */
bool properIntersect(lineSegment2D ray, lineSegment2D seg, double mRay, double bRay);

#endif
