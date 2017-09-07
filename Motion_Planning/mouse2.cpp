/* mouse2.c

Laura Toma

What it does:  The user can enter a polygon by clicking on the mouse.


*/

#include "geom.h"
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstdlib>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <cfloat>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#include <vector>

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

GLint fillmode = 0;



/* forward declarations of functions */
void display(void);
void keypress(unsigned char key, int x, int y);
void SpecialInput(int key, int x, int y);
void mousepress(int button, int state, int x, int y);
void timerfunc();

void initialize_polygon();
void print_polygon(vector<point2D> poly);
void setDirection();
bool is_inside(point2D point);
bool simple();

void draw_red(vector<point2D> poly, double guard_x, double guard_y);
void setReflexAngles();
void printReflex();
void findAllGreen(int mouse_x, int mouse_y);
void getGreen(int prevRed, int curRed, int nextRed, int guard_x, int guard_y, int redPointsIndex);
void getClosest(int guard_x, int guard_y, bool direction, int redPointsIndex);
void move_guard();


/* our coordinate system is (0,0) x (WINDOWSIZE,WINDOWSIZE) where the
   origin is the lower left corner */


/* global variables */
const int WINDOWSIZE = 700;

//the current polygon
vector<point2D>  poly;

//coordinates of last mouse click
double mouse_x=-10, mouse_y=-10;
//initialized to a point outside the window

//for moving the guard
double directions[2] = {-0.1, 0.1};
double dx = 1.0;
double dy = 0.0;
double r = sqrt(1.0);

//when this is 1, then clicking the mouse results in those points being stored in poly
int poly_init_mode = 0;

//global variable to actiavte the movement of the guard
//this can be toggled by hitting b
bool guard_moving = false;

//Number of green points added to VP
int numGreen = 0;

//True if points were made in ccw order
bool ccw = true;

//Stores red points, green points, and visible points
vector<point2D> red_points;
vector<point2D> VP;
vector<point2D> green_points;


void draw_circle(double x, double y){
  int precision = 100;
  double r = 4;
  double theta = 0;
  glBegin(GL_POLYGON);
  for(int i = 0; i < precision; i++){
    theta = i * 2 * M_PI/precision;
    glVertex2f(x + r*cos(theta), y + r*sin(theta));
  }
  glEnd();
}

/*
Usage

void glutMouseFunc(void (*func)(int button, int state, int x, int y));

Description

glutMouseFunc sets the mouse callback for the current window. When a
user presses and releases mouse buttons in the window, each press and
each release generates a mouse callback. The button parameter is one
of GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, or GLUT_RIGHT_BUTTON. For
systems with only two mouse buttons, it may not be possible to
generate GLUT_MIDDLE_BUTTON callback. For systems with a single mouse
button, it may be possible to generate only a GLUT_LEFT_BUTTON
callback. The state parameter is either GLUT_UP or GLUT_DOWN
indicating whether the callback was due to a release or press
respectively. The x and y callback parameters indicate the window
relative coordinates when the mouse button state changed. If a
GLUT_DOWN callback for a specific button is triggered, the program can
assume a GLUT_UP callback for the same button will be generated
(assuming the window still has a mouse callback registered) when the
mouse button is released even if the mouse has moved outside the
window.
*/
/* ************************************************************ */
/* Allows user to click anywhere and store the point  */
void mousepress(int button, int state, int x, int y) {

    if(state == GLUT_DOWN) {
        point2D p;
        p.x = x;
        p.y = WINDOWSIZE - y;
        //make sure the user clicked inside the polygon
        //if the program is in poly_init_mode == 1, the is_inside will return
        //true no matter what because there is no polygon yet to look inside
        if (is_inside(p)) {
            mouse_x = x;
            mouse_y = y;
            //(x,y) are in wndow coordinates, where the origin is in the upper
            //left corner; our reference system has the origin in lower left
            //corner, this means we have to reflect y
            mouse_y = WINDOWSIZE - mouse_y;

            printf("mouse click at (x=%d, y=%d)\n", (int)mouse_x, (int)mouse_y);

            if (poly_init_mode == 1) {
                point2D p = {mouse_x, mouse_y};
                poly.push_back(p);
            }
        }
    }
    glutPostRedisplay();
}

/* ************************************************************ */
/* initialize  polygon stored in global variable poly  */
void initialize_polygon() {

  //clear the vector, in case something was there
  poly.clear();

  int n = 10; //size of polygon
  double rad = 100;
  double step = 2 * M_PI / n;
  point2D p;
  for (int i=0; i<n; i++) {

    p.x = WINDOWSIZE/2 + rad * cos (i * step);
    p.y = WINDOWSIZE/2 + rad * sin (i * step);


    //insert the segment in the array of segments
    poly.push_back(p);
  } //for i
}

/* ************************************************** */
void print_polygon(vector<point2D> poly) {

  for (unsigned int i=0; i<poly.size()-1; i++) {
    printf("edge %d: [(%d,%d), (%d,%d)]\n",
	   i, poly[i].x, poly[i].y, poly[i+1].x, poly[i+1].y);
  }
  //print last edge from last point to first point
  int last = poly.size()-1;
    printf("edge %d: [(%d,%d), (%d,%d)]\n",
	   last, poly[last].x, poly[last].y, poly[0].x, poly[0].y);

}

/* ****************************** */
int main(int argc, char** argv) {

  initialize_polygon();
  print_polygon(poly);


  /* initialize GLUT  */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
  glutInitWindowPosition(100,100);
  glutCreateWindow(argv[0]);

  /* register callback functions */
  glutDisplayFunc(display);
  glutKeyboardFunc(keypress); //for all normal keys
  glutSpecialFunc(SpecialInput); //for arrow keys
  glutMouseFunc(mousepress);
  glutIdleFunc(timerfunc);

  /* init GL */
  /* set background color black*/
  glClearColor(0, 0, 0, 0);
  /* here we can enable depth testing and double buffering and so
     on */

  /* give control to event handler */
  glutMainLoop();
  return 0;
}


/* ****************************** */
/* draw the polygon */
/* has a special case that draws the polygon filled in if 'fill' is set to true */
void draw_polygon(vector<point2D> poly, bool fill){

    if (poly.size() == 0) return;

    //check if the this polygon should be filled
    if (fill) {
        //if so, set color to cyan and set transparency
        glColor4f(cyan[0], cyan[1], cyan[2], 0.1);
        //create a point for the guard
        point2D guard;
        guard.x = mouse_x;
        guard.y = mouse_y;
        //run through the polygon (visible polygon in this case)
        for (int i = 0; i < poly.size(); i++) {
            //do wraparound for next point
            int j = i + 1;
            if (j == poly.size()) {
                j = 0;
            }

            //draw a triangle from the guard to the two neighboring points
            //in the vector
            glBegin(GL_TRIANGLES);
            glVertex2f(guard.x, guard.y);
            glVertex2f(poly[i].x, poly[i].y);
            glVertex2f(poly[j].x, poly[j].y);
            glEnd();
        }
    } else {
        //draw the outline of the polygon
        //set color to yellow
        glColor3fv(yellow);
        int i;
        for (i=0; i<poly.size()-1; i++) {
            glBegin(GL_LINES);
            glVertex2f(poly[i].x, poly[i].y);
            glVertex2f(poly[i+1].x, poly[i+1].y);
            glEnd();
        }
        //render last segment between last point and forst point
        int last=poly.size()-1;
        glBegin(GL_LINES);
        glVertex2f(poly[last].x, poly[last].y);
        glVertex2f(poly[0].x, poly[0].y);
        glEnd();
    }
}





/* ****************************** */
void display(void) {

  glClear(GL_COLOR_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity(); //clear the matrix


  /* The default GL window is [-1,1]x[-1,1] with the origin in the
     center.

     Our system of coordinates (in which we generate our points) is
     (0,0) to (WINSIZE,WINSIZE), with the origin in the lower left
     corner.

     We need to map the points to [-1,1] x [-1,1]

     Assume we are the local coordinate system.

     First we scale down to [0,2] x [0,2] */
  glScalef(2.0/WINDOWSIZE, 2.0/WINDOWSIZE, 1.0);
   /* Then we translate so the local origin goes in the middle of teh
     window to (-WINDOWSIZE/2, -WINDOWSIZE/2) */
  glTranslatef(-WINDOWSIZE/2, -WINDOWSIZE/2, 0);

  //now we draw in our local coordinate system (0,0) to
  //(WINSIZE,WINSIZE), with the origin in the lower left corner.
  draw_polygon(poly, false);
  simple();

  point2D point;
  point.x = mouse_x;
  point.y = mouse_y;

  //Not drawing and guard is inside polygon
  if (poly_init_mode == 0 && is_inside(point)) {
      //Clear all vectors
      red_points.clear();
      VP.clear();
      green_points.clear();

      //Draw red points,
      draw_red(poly, mouse_x, mouse_y);

      //Determine CCW vs CW and determine reflex angles
      setDirection();
      setReflexAngles();

      //Find all green points from rays through reflex angles
      findAllGreen(mouse_x, mouse_y);

      //Draw visible polygon
      draw_polygon(VP, true);
      numGreen = 0;
      //draw a circle where the mouse was last clicked. Note that this
      //point is stored as a global variable and is modified by the mouse
      //handler function
      glColor3fv(blue);
      draw_circle(mouse_x, mouse_y);
  }

  /* execute the drawing commands */
  glFlush();
}

//Using the ray-casting algorithm, determine if a point is inside the polygon
//This technique "casts" a ray from the point, to the right, and then counts
//the intersections of that ray with the polygon. If the number of intersections
//is even, the point is outside the polygon. if its odd, the point is inside.
bool is_inside(point2D point) {
    //this is just to allow drawing when in drawing mode (affects only the call
    //to is_inside within the mouse press handler)
    if (poly_init_mode == 1) {
      return true;
    }

    //create the ray
    lineSegment2D ray;
    ray.p1.x = point.x;
    ray.p1.y = point.y;
    ray.p2.x = WINDOWSIZE;
    ray.p2.y = point.y;

    //make a new segment that will be built from the vertices of the polygon
    lineSegment2D seg;
    int numIntersections = 0;
    for (int i = 0; i < (poly.size() - 1); i++) {
        seg.p1 = poly[i];
        seg.p2 = poly[i+1];
        //count intersections
        if (intersect(ray, seg)) {
            numIntersections++;
        }
    }
    //deal with wrap-around out here by checking the last vertex with the first
    //vertex for intersection
    seg.p1 = poly[poly.size() - 1];
    seg.p2 = poly[0];
    if (intersect(ray, seg)) {
        numIntersections++;
    }

    //if number of intersections was even, return false; the point is outside.
    //else return true; the point is inside
    if (numIntersections % 2 == 0) {
        return false;
    } else {
        return true;
    }
}

/* ************************************************************* */
/* Determines if polygon is simple*/
bool simple(){
    bool isSimple = true;
    for(int i=0; i < poly.size(); i++){
        lineSegment2D seg1;
        seg1.p1 = poly[i];

        //Wrap around
        int next1 = i+1;
        if(next1 == (poly.size())){
            next1 = 0;
        }
        seg1.p2 = poly[next1];

        for(int j=0; j < poly.size(); j++){
            //Create segment
            lineSegment2D seg2;
            seg2.p1 = poly[j];

            //Wrap around
            int next2 = j+1;
            if(next2 == (poly.size())){
                next2 = 0;
            }
            seg2.p2 = poly[next2];

            //Same point or neighboring segments, skip
            if(i == j || i == next2 || j == next1){
                continue;
            }

            if (intersect(seg1, seg2)){
                glColor3fv(green);
                glBegin(GL_LINES);
                glVertex2f(seg1.p1.x, seg1.p1.y);
                glVertex2f(seg1.p2.x, seg1.p2.y);
                glEnd();
                glBegin(GL_LINES);
                glVertex2f(seg2.p1.x, seg2.p1.y);
                glVertex2f(seg2.p2.x, seg2.p2.y);
                glEnd();
                isSimple = false;
            }
        }
    }
    return isSimple;
}

/* ************************************************************* */
/* Determines if points were added clockwise or counterclockwise */
void setDirection(){
  int lowestInd = 0;
  //Find lowest point (yvalue)
  for(int i=1; i<poly.size(); i++){
    if(poly[i].y < poly[lowestInd].y){
      lowestInd = i;
    }
  }
  //Hanlde wrap around
  int nextInd = lowestInd + 1;
  int prevInd = lowestInd - 1;
  if(nextInd >= poly.size()){
    nextInd = 0;
  }
  if(prevInd < 0){
    prevInd = poly.size()-1;
  }
  //Previous point is to the left of cur->next (CCW)
  if(left(poly[lowestInd], poly[nextInd], poly[prevInd])){
    ccw = true;
  //Previous point is to the right of cur->next (CW)
  } else{
    ccw = false;
  }
}
/* ************************************************************* */
// Determine all reflex angles (based off of cw of ccw)
// Store whether angle is reflex in point struct within poly
void setReflexAngles(){
    int prev, next;
    for(int i=0; i< poly.size(); i++){

      //Allow wrap around
      prev = i-1;
      next = i+1;
      if(i==0){
        prev = poly.size()-1;
      }else if(i==poly.size()-1){
        next = 0;
      }

      //Reflex angle if points are to the left and clockwise
      if(left(poly[prev], poly[i], poly[next])){
        poly[i].reflex = !ccw;
      //Reflex angle if points are to the right and counterclockwise
      }else{
        poly[i].reflex = ccw;
      }
    }
}
/* ************************************************************* */
/* Print boolean values for whether angle is reflex or not */
/* Used to debug mostly */
void printReflex(){
  for(int i=0; i< poly.size(); i++){
    cout << poly[i].reflex << endl;
  }
}

/* Draws all red points by checking if rays from the guard intersect points */
/* the red points are defined by the vertices (reflex and non-reflex) that are */
/* visible to the guard. the green-points (found later) are the points that */
/* are visible past reflex angle vertices. */
void draw_red(vector<point2D> poly, double guard_x, double guard_y) {
    //For each point in polygon, create segment from guard to point
    for (int i = 0; i < poly.size(); i++) {
        lineSegment2D seg1;
        seg1.p1.x = guard_x;
        seg1.p1.y = guard_y;
        seg1.p2 = poly[i];
        bool intersected = false;

        //Check intersection points with every other segment in the polygon
        for (int j = 0; j < poly.size(); j++) {
            //wrap around to find end point if neccessary
            double endPt;
            if (j+1 < poly.size()) {
                endPt = j+1;
            } else {
                endPt = 0;
            }
            //make sure it's not the same segment
            if (i != j && i != endPt) {
                //make segment
                lineSegment2D seg2;
                seg2.p1 = poly[j];
                seg2.p2 = poly[endPt];
                //check intersection
                if(intersect(seg1, seg2)){
                    //once we've found an intersection,
                    //set the bool to true and get out of the loop
                    intersected = true;
                    break;
                }
            }
        }

        //make sure that the ray intersects at least one vertex
        if(!intersected) {
            point2D p = poly[i];
            //Set index of point to link this to index in poly array
            p.index = i;
            //Add red points to a red points vector and VP
            red_points.push_back(p);
            VP.push_back(p);
            glColor3fv(red);
            draw_circle(p.x, p.y);
        }
    }
}
/* ******************************************************** */
/* Find all green points associated with red reflex angles*/
void findAllGreen(int mouse_x, int mouse_y){
  for(int i=0; i < red_points.size(); i++){
      //Allow wrap around
      int prev = i-1;
      int next = i+1;
      if(i==0){
        prev = red_points.size()-1;
      }else if(i==red_points.size()-1){
        next = 0;
      }
    //Only look at reflex angles
      if(poly[red_points[i].index].reflex){
          //Pass in index of point in poly array
          getGreen(red_points[prev].index, red_points[i].index, red_points[next].index, mouse_x, mouse_y, i);
      }
  }
}

/* ****************************************************************************************** */
/* Finds the closest intersection (green point), draws it, and adds it to the VP */
void getClosest(int guard_x, int guard_y, bool direction, int redPointsIndex){
  //Get closest point
  double distance = 0.0;
  double minDistance = DBL_MAX;
  int greenIndex = 0;
  for(int i=0; i < green_points.size(); i++){
    distance = sqrt((guard_x - green_points[i].x)*(guard_x - green_points[i].x)
                    + (guard_y - green_points[i].y)*(guard_y - green_points[i].y));
    //Shortest distance so far
    if(distance < minDistance){
      minDistance = distance;
      greenIndex = i;
    }
  }

  //Insert point after red point in VP
  if(direction == true){
    VP.insert(VP.begin() + redPointsIndex + 1 + numGreen, green_points[greenIndex]);
  //Insert point before red point
  }else{
    VP.insert(VP.begin() + redPointsIndex + numGreen, green_points[greenIndex]);
  }

  //Draw green point and then clear all green points in vector
  numGreen++;
  glColor3fv(cyan);
  draw_circle(green_points[greenIndex].x, green_points[greenIndex].y);
  green_points.clear();
  return;
}
/* ****************************************************************************************** */
/* Find the green point from the ray associated with a red reflex angle*/
/* curRed is the index in poly associated with the red reflex point, prev is index of previous RED point, */
/* next is index of next RED point. redPointsIndex is the index of the current red point in the red_points vector */
void getGreen(int prevRed, int curRed, int nextRed, int guard_x, int guard_y, int redPointsIndex){
  //signifies whether green intersection point has been found
  bool greenFound = false;

  //determines if the green point was found iterating forward or backward
  bool forward = true;

  //skip red reflex point if the guard is on that point
    if((guard_x == poly[curRed].x) && (guard_y == poly[curRed].y)){
      return;
    }

    //Create ray from guard through red reflex angle to edge of window
    lineSegment2D ray;
    ray.p1.x = guard_x;
    ray.p1.y = guard_y;

    //Guard and red reflex point have ame x coords, different y (draw ray vertically up or down)
    if(guard_x == poly[curRed].x){
      //Draw ray vertically upward
      if(guard_y < poly[curRed].y){
        ray.p2.x = guard_x;
        ray.p2.y = WINDOWSIZE;
      //Draw ray vertically downward
      }else{
        ray.p2.x = guard_x;
        ray.p2.y = 0;
      }
    //Different x coords
    } else{
        //Compute slope m and intercept b of ray (goes to edge of window)
        double mRay = (double)(poly[curRed].y - guard_y)/ (double)(poly[curRed].x - guard_x);
        double bRay = (double)((double)guard_y - mRay*(double)guard_x);
        //Ray to top edge
        if(((mRay < (-1.0))||(mRay > 1.0)) && (poly[curRed].y > guard_y)){
          ray.p2.y = WINDOWSIZE;
          ray.p2.x = (((double)ray.p2.y - bRay) / mRay);
        //Ray to bottom edge
        } else if(((mRay < (-1.0))||(mRay > 1.0)) && (poly[curRed].y < guard_y)){
          ray.p2.y = 0;
          ray.p2.x = (double)((-1.0*bRay) / mRay);
        //Ray to right edge
        } else if(poly[curRed].x > guard_x){
          ray.p2.x = WINDOWSIZE;
          ray.p2.y = (mRay * (double)ray.p2.x + bRay);
        //Ray to left edge
        } else{
          ray.p2.x = 0;
          ray.p2.y = (int)bRay;
        }
    }

  //Iterate forward to next red point and check intersections with all segments in poly
  int i = curRed;
  while(i != nextRed){
    //Skip first segment
    i++;

    //Next point in poly
    int next = i+1;

    //Handle wrap around
    if(i == poly.size()){
      i = 0;
      next = i+1;
    }else if(next == (poly.size())){
      next = 0;
    }

    //Create segment from adjacent points in poly
    lineSegment2D seg;
    seg.p1 = poly[i];
    seg.p2 = poly[next];

    //Ray and segment intersect
    if(intersect(ray, seg)){
      //Create intersection point (green point)
      point2D green;
      double mSeg = (double)(seg.p2.y - seg.p1.y) / (double)(seg.p2.x - seg.p1.x);
      double bSeg = (double)seg.p2.y - mSeg*(double)seg.p2.x;

      //Recompute slope and intercept of ray since they were created in an if statement
      double mRay = (double)((poly[curRed].y - guard_y)/ (double)(poly[curRed].x - guard_x));
      double bRay = (double)guard_y - mRay*(double)guard_x;

      //Check if segments PROPERLY intersect such that ray doesn't leave polygon
      int preCur = curRed - 1;
      if(preCur == -1){
        preCur = poly.size()-1;
      }
      int nextCur = curRed + 1;
      if(nextCur == poly.size()){
        nextCur = 0;
      }
      //Ray leaves polygon (goes between two segments on either side of current red reflex point)
      if(left(ray.p1, ray.p2, poly[preCur]) != left(ray.p1, ray.p2, poly[nextCur])){
        return;
      }

      //ray is vertical
      if(ray.p1.x == ray.p2.x){
        green.x = ray.p1.x;
        green.y = mSeg * (double)green.x + bSeg;
      //segment is vertical
      }else if(seg.p1.x == seg.p2.x){
        green.x = seg.p1.x;
        green.y = mRay * (double)green.x + bRay;
      //Neither ray nor segment is vertical
      } else{
        green.x = (bRay - bSeg) / (mSeg - mRay);
        green.y = mSeg * (double)green.x + bSeg;
      }

      //Add all green points (intersections inside the polygon) to an array (will then find closest intersection)
      green_points.push_back(green);
      greenFound = true;
      forward = true;
    }
  }

   //Iterate backward to previous red point and check intersections with all segments in poly
  i = curRed;
  while(i != prevRed){
    //Skip first point
    i--;
    //previous point in poly
    int prev = i-1;

    //Handle wrap around
    if(i == -1){
      i = poly.size() - 1;
      prev = i-1;
    }else if(prev == -1){
      prev = poly.size() - 1;
    }

    //Create segment from points in array
    lineSegment2D seg;
    seg.p1 = poly[i];
    seg.p2 = poly[prev];

    //Determine if segment and ray intersect
    if(intersect(ray, seg)){
      //Create intersection point (green point)
      point2D green;

      //Compute slope of segment
      double mSeg = (double)(seg.p2.y - seg.p1.y) / (double)(seg.p2.x - seg.p1.x);
      double bSeg = (double)seg.p2.y - mSeg*(double)seg.p2.x;

      //Recompute slope of ray
      double mRay = (double)((poly[curRed].y - guard_y)/ (double)(poly[curRed].x - guard_x));
      double bRay = (double)guard_y - mRay*(double)guard_x;

      //Check if segments PROPERLY intersect such that ray doesn't leave polygon
      int preCur = curRed - 1;
      if(preCur == -1){
        preCur = poly.size()-1;
      }
      int nextCur = curRed + 1;
      if(nextCur == poly.size()){
        nextCur = 0;
      }

      //Ray leaves polygon (goes between two segments on either side of current red reflex point)
      if(left(ray.p1, ray.p2, poly[preCur]) != left(ray.p1, ray.p2, poly[nextCur])){
        return;
      }

      //ray is vertical
      if(ray.p1.x == ray.p2.x){
        green.x = ray.p1.x;
        green.y = mSeg * green.x + bSeg;
      //segment is vertical
      }else if(seg.p1.x == seg.p2.x){
        green.x = seg.p1.x;
        green.y = mRay * green.x + bRay;
      //Neither ray nor segment is vertical
      } else{
        green.x = (bRay - bSeg) / (mSeg - mRay);
        green.y = mSeg * green.x + bSeg;
      }

      //Add all green points (intersections) to an array (will then find closest intersection)
      green_points.push_back(green);
      greenFound = true;
      forward = false;
    }
  }

  //Green point was found, find closest intersection and add to VP
  if(greenFound == true){
    getClosest(guard_x, guard_y, forward, redPointsIndex);
    return;
  }
}

/*****************************/
//A special funciton do deal with arrow key inputs for moving
//the guard.
void SpecialInput(int key, int x, int y) {
    point2D guard;
    //figure out which key is pressed, and act accordingly
    //in each case set the guard-moving to false, and make
    //sure that moving it will not put it outside the poly
    switch (key) {
        case GLUT_KEY_LEFT:
            //move guard to the left
            guard_moving = false;
            mouse_x += -3;
            guard.x = mouse_x;
            guard.y = mouse_y;
            //make sure not outside with previous move
            if (!is_inside(guard)) {
                mouse_x -= -3;
            }
            break;
        case GLUT_KEY_RIGHT:
            //move guard to the right
            guard_moving = false;
            mouse_x += 3;
            guard.x = mouse_x;
            guard.y = mouse_y;
            //make sure not outside with previous move
            if (!is_inside(guard)) {
                mouse_x -= 3;
            }
            break;
        case GLUT_KEY_UP:
            //move guard up
            guard_moving = false;
            mouse_y += 3;
            guard.x = mouse_x;
            guard.y = mouse_y;
            //make sure not outside with previous move
            if (!is_inside(guard)) {
                mouse_y -= 3;
            }
            break;
        case GLUT_KEY_DOWN:
            //move guard down
            guard_moving = false;
            mouse_y += -3;
            guard.x = mouse_x;
            guard.y = mouse_y;
            //make sure not outside with previous move
            if (!is_inside(guard)) {
                mouse_y -= -3;
            }
            break;
    }

    //re-draw it all
    glutPostRedisplay();
}

/* ****************************** */
void keypress(unsigned char key, int x, int y) {
    switch(key) {
        case 'q':
        exit(0);
        break;

    //expected behaviour: press 's', then click on the points you
    //want, and press 'e' when you're done. the points will be saved
    //in global 'poly'

    case 's':
        //start re-initializeing polygon
        //turn-off guard moving
        guard_moving = false;
        //re-set all vectors
        poly.clear();
        VP.clear();
        red_points.clear();
        green_points.clear();
        //re-set the position of the guard
        mouse_x = mouse_y=0;
        //change to drawing mode
        poly_init_mode = 1;
        glutPostRedisplay();
        break;

    case 'b':
        //this is the case that allows the guard to move
        //but it is a key that toggles that ability so we
        //negate the previous value
        guard_moving = !guard_moving;

    case 'e':
        if(!simple()){
          cout << "Polygon is not simple. Please erase, ('s') and draw another one." << endl;
          break;
        }
        poly_init_mode=0;
        glutPostRedisplay();
        break;
    }
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

void move_guard() {
    //create a checker point that will be the candidate move for the guard.
    point2D checker;
    checker.x = mouse_x + dx;
    checker.y = mouse_y + dy;
    //make sure that the move does not put the guard outside the poly
    while(!is_inside(checker)) {
        //if it does put the guard outside, re-set the checker
        //and get a new random direction
        checker.x = mouse_x;
        checker.y = mouse_y;
        //also the randomness is minimized by dividing by two
        //dx = -(directions[(rand() % 2)] - ((double)rand()/RAND_MAX)/2);
        //dy = -(directions[(rand() % 2)] - ((double)rand()/RAND_MAX)/2);
        if(dx > 0.0){
          dx = -(directions[(rand() % 2)] + ((double)rand()/RAND_MAX)/10);
        } else{
          dx = -(directions[(rand() % 2)] - ((double)rand()/RAND_MAX)/10);
        }
        //Last dy was positive
        if(dy > 0.0){
          dy = -sqrt(r*r - dx*dx);
        //Last dy was negative
        } else {
          dy = sqrt(r*r - dx*dx);
        }

        //update the checker, and then check if it's inside now
        checker.x += dx;
        checker.y += dy;
    }
    //once the direction allows us to stay inside the poly,
    //we can update the position of the guard
    mouse_x += dx;
    mouse_y += dy;
}


void timerfunc() {
  //check if it's appropriate to move the guard
  //if so, move the guard
  if (poly_init_mode == 0 && guard_moving) {
      move_guard();
      glutPostRedisplay();
  }

}
