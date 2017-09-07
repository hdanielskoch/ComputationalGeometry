#include "geom.h"
#include "pointStack.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

/* **************************************** */
/* returns the signed area of triangle abc. The area is positive if c
 is to the left of ab, and negative if c is to the right of ab
 */
double signed_area2D(point2D a, point2D b, point2D c) {
	return ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}



/* **************************************** */
/* return 1 if p,q,r collinear, and 0 otherwise */
int collinear(point2D p, point2D q, point2D r) {
	return (signed_area2D(p, q, r) == 0);
}

/* **************************************** */
/* return 1 if same point*/
int samePoint(point2D a, point2D b){
	return ((a.x == b.x) && (a.y == b.y));
}

/* **************************************** */
/* return 1 if c is  strictly left of ab; 0 otherwise */
int left (point2D a, point2D b, point2D c) {
	return (signed_area2D(a, b, c) > 0);
}

/* **************************************** */
/* returns the index of the lowest point (smallest y value)
   If same y value, take furthest left point (smallest x value) */
int lowestInd(point2D* p, int n){
	int lowInd = INT_MAX;
	int lowPointX = INT_MAX; 
	int lowPointY = INT_MAX;
	int i;
	for(i=0; i<n; i++){
		
		//lowest point so far
		if(p[i].y < lowPointY){
			lowInd = i;
			lowPointY = p[i].y;
			lowPointX = p[i].x;
		}
		//Same y coordinate => choose farthest left x coord
		else if((p[i].y == lowPointY) && (p[i].x < lowPointX)){
				lowInd = i;
				lowPointY = p[i].y;
				lowPointX = p[i].x;
		}
	}
	return lowInd;
}

/* compute the convex hull of the points in p; the points on the CH are returned as a list
 */
pointNode* graham_scan(point2D* p, int n) {
	pointStack* stack = createPointStack();
	//Find lowest point
	int lowInd = lowestInd(p, n);
	int xLow = p[lowInd].x;
	int yLow = p[lowInd].y;
	
	double PI = acos(-1.0);
	double r = 0.0; //radius from lowest point
	double dx = 0.0; //Difference in x coords between the lowest point and some other point
	double dy = 0.0;
	//Stores angle in array
	point2D* angArr = (point2D*)malloc(n*sizeof(point2D));
	//point2D angArr[n];
	
	//Add angle to p array
	int i;
	for(i=0; i<n; i++){
		//Skip lowest point and all points equal to lowest point
		if((p[i].x == xLow) && (p[i].y == yLow)){

			//This will automatically be the lowest angle
			angArr[i].angle = -1.0; 
			angArr[i].x = p[i].x;
			angArr[i].y = p[i].y;
			continue;
		}

		dx = (double)(p[i].x - p[lowInd].x);
		dy = (double)(p[i].y - p[lowInd].y);
		r = sqrt(dx*dx + dy*dy);
		
		//Set angle in degrees
		angArr[i].angle = acos(dx/r) * 360.0 / (2.0*PI);
		angArr[i].x = p[i].x;
		angArr[i].y = p[i].y;
	}
	
	//Sort array by angle, then by x if same angle, then by y if same y
	qsort(angArr, n, sizeof(point2D), cmpfunc);
	
	//Initialize two pointers to move around stack
	pointNode *p1 = (pointNode *) malloc(sizeof(pointNode));
	p1->p = angArr[0];

	//Iterate until a second DIFFERENT point is found
	pointNode *p2 = (pointNode *) malloc(sizeof(pointNode));
	p2->p = angArr[1];

	i = 2;
	while(samePoint(p1->p, p2->p) == 1){
		p2->p = angArr[i];
		i++;
	}

	//First two points in angle array must be the first two points in CH
	push(stack, p1);
	push(stack, p2);

	int j = i; //Start iterating from i
	
	//Do single pass over to construct convex hull
	for (i = j; i < n; i++){

		//p2 is top of stack and p1 is below it
		p1 = stack->firstPoint->next;
		p2 = stack->firstPoint;
		
		//Next point in consideration
		pointNode *p3 = (pointNode *) malloc(sizeof(pointNode));
		p3->p = angArr[i];

		// if((samePoint(p1->p, p2->p) == 1) || (samePoint(p2->p, p3->p) == 1) || (samePoint(p1->p, p3->p) == 1)){
		// 	//printf("here");
		// 	continue;
		// }

		if((angArr[i].x == angArr[i-1].x) && (angArr[i].y == angArr[i-1].y) && (angArr[i].angle == angArr[i-1].angle)){
			//printf("here");
			continue;
		}

		//Point is left of previous two points or on same line and p2 and p3 are not the same points
		if (left(p1->p, p2->p, p3->p) || collinear(p1->p, p2->p, p3->p)) {
			push(stack, p3);
			printf("%d, %d\n",p3->p.x, p3->p.y);
		//Point is to the right of previous two points
		} else {
			pop(stack);
			i--;
		}
	}
	return stack->firstPoint;
}

//adapted from http://stackoverflow.com/questions/6103636/c-qsort-not-working-correctly
//this compare function is
/* Returns 1 if A should appear before B,
   returns -1 otherwise */
int cmpfunc (const void * a, const void * b) {
	point2D *A = (point2D *)a;
	point2D *B = (point2D *)b;
	
	//Compare angles
	if ((A->angle - B->angle) > 0.0)
		return 1;
	else if ((A->angle - B->angle) < 0.0)
		return -1;
	//Angles are the same, compare x values
	else if (A->x < B->x)
		return -1;
	else if (A->x > B->x)
		return 1;
	//angles and x values are the same, compare y values
	else if (A->y < B->y)
		return 1;
	else
		return -1;
}

