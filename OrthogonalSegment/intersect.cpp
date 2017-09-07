/*
 
 intersect.cpp
 
 Laura Toma
 
 What it does:
 
 Draws a set of horizontal and vertical line segments in the default 2D
 projection. Then it pretends to compute their intersections using the
 line sweep algorithm, and simulates the sweep line moving from left to
 right.
 
 */

#include "geom.h"
#include "rtimer.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <assert.h>
#include <algorithm>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif



#include <vector>
#include <map>
using namespace std;



GLfloat red[3] = {1.0, 0.0, 0.0};
GLfloat green[3] = {0.0, 1.0, 0.0};
GLfloat blue[3] = {0.0, 0.0, 1.0};
GLfloat black[3] = {0.0, 0.0, 0.0};
GLfloat white[3] = {1.0, 1.0, 1.0};
GLfloat gray[3] = {0.5, 0.5, 0.5};
GLfloat yellow[3] = {1.0, 1.0, 0.0};
GLfloat magenta[3] = {1.0, 0.0, 1.0};
GLfloat cyan[3] = {0.0, 1.0, 1.0};

//struct that will be sorted, containing x-values and their associated segments
typedef struct _xSeg {
	segment2D seg;
	int x;
} xSeg;

//adapted from http://stackoverflow.com/questions/1380463/sorting-a-vector-of-custom-objects
struct less_than_key {
	inline bool operator() (const xSeg& xSeg1, const xSeg& xSeg2) {
		if (xSeg1.x < xSeg2.x){
			return true;
		} else if(xSeg1.x > xSeg2.x){
			return false;
			//x values are the same
		} else if (xSeg1.seg.start.y > xSeg2.seg.start.y){
			return true;
		} else if (xSeg1.seg.start.y < xSeg2.seg.start.y){
			return false;
		} else if (xSeg1.seg.end.y > xSeg2.seg.end.y){
			return true;
		} else if (xSeg1.seg.end.y < xSeg2.seg.end.y){
			return false;
			//Segments are the same
		} else
			return false;
	}
};


/* forward declarations of functions */
void display(void);
void keypress(unsigned char key, int x, int y);
void timerfunc();

void initialize_segments_random();
void initialize_segments_horizontal();
void initialize_segments_helix();
void initialize_segments_zero();
void initialize_segments_worst_case();
void init_segs_hash();
void print_segments();

//renders the sweep line
void draw_sweep_line();

//renders the active structure
void draw_active_structure();

//renders the intersection points
void draw_intersection_points();

vector<xSeg> initializeListSegements();
void computeIntersection(xSeg x);
bool vertical(segment2D l);
void searchAS(int x, int y1, int y2);




/* global variables */
const int WINDOWSIZE = 500;

//we got two test cases so far: random and horizontal; add more!
int init_case = 0;
const int NB_TEST_CASES = 6;

//NOTE: all the structures below need to be global so that they can be rendered

//number of segments requested by user
int n;

//the array of  segments
vector<segment2D>  segments;

//the active structure that stores the segments intersecting the sweep line
map <int,segment2D> as;

vector<xSeg> X;

vector<xSeg>::iterator nextX;

//the intersections points of the segments
vector<point2D> intpoints;


//the events. do we need this??
//vector<sweepEvent> events;


//current position of sweep line; this is used to animate the sweep line moving
int sweep_line_x = 0;

/* ************************************************** */
//fills global variable "segments" with n segments
void initialize_segments() {
	switch (init_case)  {
			
		case 0:
			initialize_segments_random();
			break;
			
		case 1:
			initialize_segments_horizontal();
			break;
		case 2:
			initialize_segments_zero();
			break;
		case 3:
			initialize_segments_helix();
			break;
		case 4:
			init_segs_hash();
			break;
		case 5:
			initialize_segments_worst_case();
			break;
			
		default:
			initialize_segments_random();
	}
	
	sweep_line_x = -1;
	X.clear();
	X = initializeListSegements();
	nextX = X.begin();
	as.clear();
	intpoints.clear();
	
	init_case = (init_case+1) % NB_TEST_CASES;
	return;
}




/* ************************************************** */
//fills global variable "segments" with n horizontal segments
void initialize_segments_horizontal() {
	int i;
	point2D a,b;
	segment2D s;
	
	//clear the vector of segments
	segments.clear();
	
	//a long horizontal segment
	a.x = 1;
	a.y = WINDOWSIZE/2;
	b.x = WINDOWSIZE - 10;
	b.y = a.y;
	
	s.start = a;
	s.end = b;
	segments.push_back(s);
	
	//n-1 vertical segments
	for (i=0; i<n-1; i++) {
		
		a.x = i*WINDOWSIZE/n;
		a.y = WINDOWSIZE/2 - random() % ((int)(.4*WINDOWSIZE));
		b.x = a.x;
		b.y = WINDOWSIZE/2 + random() % ((int)(.4*WINDOWSIZE));
		s.start = a;
		s.end = b;
		segments.push_back(s);
	}
	
}


/* ****************************** */
//fills global variable "segments" with n random segments, half horizontal and half vertical
void initialize_segments_random() {
	
	//clear the vector
	segments.clear();
	
	int i;
	point2D a, b;
	segment2D s;
	for (i=0; i<n; i++) {
		if (random()%2 == 0) {
			//horizontal segment
			a.x = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE));
			a.y =  (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE));
			b.y = a.y;
			b.x = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE));
			
			if (a.x < b.x) {
				s.start = a; s.end = b;
			} else {
				s.start = b; s.end = a;
			}
			
		} else {
			//vertical segment
			a.x = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE));
			b.x = a.x;
			a.y = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE));
			b.y = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE));
			
			if (a.y < b.y) {
				s.start = a; s.end = b;
			} else {
				s.start = b; s.end = a;
			}
		}
		
		//insert the segment in the array of segments
		segments.push_back (s);
	} //for i
}

/* ************************************************** */

//fills global variable "segments" with n segments, half horizontal and half vertical,
//  so that a horizontal line intersects all the other vertical lines
void initialize_segments_worst_case() {
	//clear the vector
	segments.clear();
	
	int i;
	int offset = 100;
	int length = 300;
	int pad = length / (n/2);
	point2D a, b;
	segment2D s;
	
	// printf("pad = %d\n", pad);
	
	int count1 = 0, count2 = 0;
	
	for (i = 0; i < n; i++) {
		
		if (i % 2 == 0) {
			// horizontal segment
			a.x = offset;
			a.y = offset + pad * count1;
			b.x = offset + length;
			b.y = offset + pad * count1;
			s.start = a;
			s.end = b;
			count1++;
		}
		else {
			// vertical segment
			a.x = offset + pad * count2;
			a.y = offset;
			b.x = offset + pad * count2;
			b.y = offset + length;
			s.start = a;
			s.end = b;
			count2++;
		}
		
		// insert the segment in the array of segments
		segments.push_back(s);
	}
}

/* Initializes the segments in a square that has zero intersections.
 * Jack Ward
 */
void initialize_segments_zero() {
	// Clear the vector
	segments.clear();
	int i;
	const int MARGIN = 80;
	const int PADDING = 4;
	point2D a, b;
	segment2D s;
	
	// Horizontal segments
	for (i = PADDING; i < n/2; i++) {
		a.x = MARGIN;
		a.y = MARGIN + PADDING * i + PADDING;
		b.x = MARGIN + PADDING * i;
		b.y = MARGIN + PADDING * i + PADDING;
		
		s.start = a;
		s.end = b;
		
		segments.push_back(s);
	}
	
	// Vertical segments
	for (i = PADDING; i < n/2; i++) {
		a.x = MARGIN + PADDING * i + PADDING;
		a.y = MARGIN;
		b.x = MARGIN + PADDING * i + PADDING;
		b.y = MARGIN + PADDING * i;
		
		s.start = a;
		s.end = b;
		
		segments.push_back(s);
	}
	
}

/* Jack and Duncan's test case*/
void init_segs_hash() {
	segments.clear();
	int i;
	int offset = 100;
	int len = 200;
	int pad = len/2;
	point2D a, b;
	segment2D s;
	int j = 0;
	for (i = 0; i < n; i++) {
		if (i % 3 == 0) {
			a.x = 40 + j;
			a.y = 20;
			b.x = 40 + j;
			b.y = 480;
			j = j + 1;
		}
		if (i % 3 == 1) {
			a.x = 38;
			b.x = 462;
			a.y = 200 + j;
			b.y = 200 + j;
			j++;
		}
		if (i % 3 == 2) {
			if (260 - j > 0) {
				a.x = 260 - j;
				a.y = 20;
				b.y = 480;
				b.x = 260 - j;
				j++;
			}
		}
		s.start = a;
		s.end = b;
		segments.push_back(s);
	}
}

void initialize_segments_helix(){
	//clear the vector
	segments.clear();
	int pad = WINDOWSIZE/(n/2 + 2);
	int length = pad * 2;
	point2D a, b;
	segment2D s;
	// printf("pad = %d\n", pad);
	int count1 = 1, count2 = 1;
	for (int i = 0; i < n; i++) {
		if (i % 2 == 0) {
			// horizontal segment
			a.y = pad * (count1 + 1);
			b.y = pad * (count1 + 1);
			a.x = pad * count1;
			b.x = a.x + length;
			count1++;
		}
		else {
			// vertical segment
			a.x = pad * (count2 + 1);
			b.x = pad * (count2 + 1);
			a.y = pad * count2;
			b.y = a.y + length;
			count2++;
		}
		s.start = a;
		s.end = b;
		// insert the segment in the array of segments
		segments.push_back(s);
	}
}


/* ************************************************** */
void print_segments() {
	
	for (int i=0; i<segments.size(); i++) {
		printf("segment %d: [(%d,%d), (%d,%d)]\n",
			   i, segments[i].start.x, segments[i].start.y, segments[i].end.x, segments[i].end.y);
		
	}
}


vector<xSeg> initializeListSegements() {
	vector<xSeg> vect;
	vect.clear();
	
	for (int i = 0; i < segments.size(); i++) {
		//check if vertical
		if (vertical(segments[i])) {
			xSeg a;
			a.x = segments[i].start.x;
			a.seg = segments[i];
			vect.push_back(a);
			//else must be a horizontal, thus add start and end points
		} else {
			xSeg a;
			a.x = segments[i].start.x;
			a.seg = segments[i];
			vect.push_back(a);
			xSeg b;
			b.x = segments[i].end.x;
			b.seg = segments[i];
			vect.push_back(b);
		}
	}
	
	//sort the vector by X value
	sort(vect.begin(), vect.end(), less_than_key());
	
	return vect;
}

/* ****************************** */
int main(int argc, char** argv) {
	
	//read number of points from user
	if (argc!=2) {
		printf("usage: viewPoints <nbPoints>\n");
		exit(1);
	}
	n = atoi(argv[1]);
	printf("you entered n=%d\n", n);
	assert(n >0);
	
	//the default is to initialize the segments randomly
	initialize_segments_random();
	
	//Rtimer rt1;
	//rt_start(rt1);
	
	//initialize the list of all events
	X.clear();
	as.clear();
	X = initializeListSegements();
	nextX = X.begin();
	
	//rt_stop(rt1);
 
	
	//print the timing
	//  char buf [1024];
	//rt_sprint(buf,rt1);
	//printf("run time:  %s\n\n", buf);
	//fflush(stdout);
	
 
	
	
	/* initialize GLUT  */
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
	glutInitWindowPosition(100,100);
	glutCreateWindow(argv[0]);
	
	/* register callback functions */
	glutDisplayFunc(display);
	glutKeyboardFunc(keypress);
	glutIdleFunc(timerfunc);  //<--------note we got an idle function, we'll use it to animate
	
	/* init GL */
	/* set background color black*/
	glClearColor(0, 0, 0, 0);
	
	
	/* give control to event handler */
	glutMainLoop();
	return 0;
}

bool vertical(segment2D l) {
	return l.start.x == l.end.x;
}

void computeIntersection(xSeg x) {
	//if vertical
	if (vertical(x.seg)) {
		int y1 = x.seg.start.y;
		int y2 = x.seg.end.y;
		
		if (y1 > y2) {
			int temp = y1;
			y1 = y2;
			y2 = temp;
		}
		
		searchAS(x.x, y1, y2);
		//if start of horizontal
	} else if (x.x == x.seg.start.x) {
		as.insert(pair<int,segment2D>(x.seg.start.y, x.seg));
	} else {
		as.erase(x.seg.start.y);
	}
}

void searchAS(int x, int y1, int y2){
	map<int,segment2D>::iterator low = as.lower_bound(y1);
	map<int,segment2D>::iterator high = as.upper_bound(y2);
	
	//Active structure contains at least one horizontal line in y range
	if(low != as.end()){
		//Iteratate through y range
		for(map<int,segment2D>::iterator it=low; it!= high; it++){
			point2D intSect;
			intSect.x = x;
			intSect.y = it->second.start.y;
			intpoints.push_back(intSect);
		}
	}
	
}

/* draw the segments stored in global variable segments
 NOTE: We draw in the local coordinate system (0,0) to (WINSIZE,WINSIZE)
 */
void draw_segments(){
	
	//set color
	glColor3fv(yellow);
	
	int i;
	for (i=0; i<segments.size(); i++) {
		glBegin(GL_LINES);
		glVertex2f(segments[i].start.x, segments[i].start.y);
		glVertex2f(segments[i].end.x, segments[i].end.y);
		glEnd();
	}
}

/*
 draw the sweep line
 NOTE: We draw in the local coordinate system (0,0) to (WINSIZE,WINSIZE)
 */
void draw_sweep_line() {
	
	//sweep line color
	glColor3fv(red);
	
	//the current position of sweep line is sweep_line_x; assume it's a
	//segment from y=0 to y=windowsize;
	glBegin(GL_LINES);
	glVertex2f(sweep_line_x, 0);
	glVertex2f(sweep_line_x, WINDOWSIZE);
	glEnd();
}


/* draw a segment with current color
 NOTE: We draw in the local coordinate system (0,0) to (WINSIZE,WINSIZE)
 */
void draw_segment(segment2D s) {
	glBegin(GL_LINES);
	glVertex2f(s.start.x, s.start.y);
	glVertex2f(s.end.x, s.end.y);
	glEnd();
}


/* draw all the elements in the active structure
 NOTE: we draw in the local coordinate system (0,0) to (WINSIZE,WINSIZE)
 */
void draw_active_structure() {
	//active structure lines color
	glColor3fv(cyan);
	
	for (map<int,segment2D>::iterator it = as.begin(); it != as.end(); it++) {
		draw_segment(it->second);
	}
	
}



/* draw all the elements in intpoints
 NOTE: we draw in the local coordinate system  (0,0) to (WINSIZE,WINSIZE)
 */
void draw_intersection_points() {
	
	const int R= 2;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	
	glColor3fv(green);
	//Iterate over all interection points
	for(int i=0; i < intpoints.size(); i++){
		//draw a small square centered at (points[i].x, points[i].y)
		glBegin(GL_POLYGON);
		glVertex2f(intpoints[i].x -R,intpoints[i].y-R);
		glVertex2f(intpoints[i].x +R,intpoints[i].y-R);
		glVertex2f(intpoints[i].x +R,intpoints[i].y+R);
		glVertex2f(intpoints[i].x -R,intpoints[i].y+R);
		glEnd();
	}
}




/* ****************************** */
void display(void) {
	
	//clear everything, start from scratch
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	
	/* The default GL window is [-1,1]x[-1,1] with the origin in the center.
	 
	 The points are in the range (0,0) to (WINSIZE,WINSIZE). This is our local coordinate system.
	 
	 We first transform our local coordinate system to [-1,1] x [-1,1]
	 */
	glScalef(2.0/WINDOWSIZE, 2.0/WINDOWSIZE, 1.0);
	glTranslatef(-WINDOWSIZE/2, -WINDOWSIZE/2, 0);
	
	
	while(sweep_line_x == nextX->x) {
		//cout << sweep_line_x << ", " << nextX->x << endl;
		computeIntersection(*nextX);
		nextX++;
	}
	
	//draw my scene in the local coordinate system (0,0) to (WINSIZE,WINSIZE)
	draw_segments();
	draw_active_structure();
	draw_intersection_points();
	draw_sweep_line();
	
	
	/* execute the drawing commands */
	glFlush();
}



/* ****************************** */
void keypress(unsigned char key, int x, int y) {
	switch(key) {
		case 'q':
			exit(0);
			break;
			
		case 'i':
			initialize_segments();
			glutPostRedisplay();
			break;
	}
}


void timerfunc() {
	
	sweep_line_x++;
	
	glutPostRedisplay();
}


/* Handler for window re-size event. Called back when the window first appears and
 whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
	
	// Set the viewport to cover the new window
	glViewport(0, 0, width, height);
 
	glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
	glLoadIdentity();             // Reset
	gluOrtho2D(0.0, (GLdouble) width, 0.0, (GLdouble) height);
}




