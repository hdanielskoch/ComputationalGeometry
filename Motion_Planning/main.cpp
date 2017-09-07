/* main.cpp

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
#include <climits>
#include <assert.h>
#include <cfloat>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#include <vector>
//#include <array>

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
void mousepress(int button, int state, int x, int y);
void timerfunc();

void initialize_polygon();
bool is_inside(point2D point, vector<point2D> poly);
bool simple(vector<point2D> poly);
bool check_seg_intersections(lineSegment2D seg1);
void check_all_points(int poly, int pt);
void check_local_points(int poly, int pt);
void calculate_VG();
void printVG();
void addBoundingBox();
bool check_seg_points_not_equal(lineSegment2D seg1, lineSegment2D seg2);
void makeAdjacency();
void printAdjList();
// void polygon_error(GLenum errno);
// void combine(GLdouble coords[3], void* vertex_data[4], GLfloat weight[4], void** outData);

//For Dijkstra's
void runDijkstra();
double get_distance(point2D a, point2D b);
int get_index_min_distance(vector<point2D> Q, vector<double> weights);
void printShortestPath();


/* our coordinate system is (0,0) x (WINDOWSIZE,WINDOWSIZE) where the
   origin is the lower left corner */


/* global variables */
const int WINDOWSIZE = 700;

//coordinates of last mouse click
double mouse_x=-10, mouse_y=-10;


//when this is 1, then clicking the mouse results in those points being stored in poly
int poly_init_mode = 0;

//When this is one, the user is done drawing polygones
int startInit = 0;
int endInit = 0;

//Stores all polygons, temporary polygon, bounding box, VG, adjacency list, and shortest path
vector<vector<point2D> > polys;
vector<point2D> polyTemp;
vector<point2D> boundBox;
vector<lineSegment2D> VG;
vector<vector<point2D> > adjList;
vector<lineSegment2D> shortestPath;

//Start and end points of VG
point2D start;
point2D endPt;

/* Draws small cirlce to indicate point */
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
        mouse_x = x;
        mouse_y = y;
        //(x,y) are in wndow coordinates, where the origin is in the upper
        //left corner; our reference system has the origin in lower left
        //corner, this means we have to reflect y
        mouse_y = WINDOWSIZE - mouse_y;
        printf("mouse click at (x=%d, y=%d)\n", (int)mouse_x, (int)mouse_y);

        //User is now entering start point
        if(startInit){
          start.x = mouse_x;
          start.y = mouse_y;
          cout << "Start: " << start.x << " " << start.y << endl;
        //User is now entering end point
        }else if(endInit){
          endPt.x = mouse_x;
          endPt.y = mouse_y;
          cout << "End: " << endPt.x << " " << endPt.y << endl;
        }
        //User is initializing polygon barriers
        //make sure the user clicked inside the polygon
        //if the program is in poly_init_mode == 1, the is_inside will return
        //true no matter what because there is no polygon yet to look inside
        else if (is_inside(p, polyTemp)) {

            printf("mouse click at (x=%d, y=%d)\n", (int)mouse_x, (int)mouse_y);

            if (poly_init_mode == 1) {
                point2D p = {mouse_x, mouse_y};

                //Add point to temporary polygon currently being constructed
                polyTemp.push_back(p);
            }
        }
    }
    glutPostRedisplay();
}

/* ****************************** */
int main(int argc, char** argv) {
  /* initialize GLUT  */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
  glutInitWindowPosition(100,100);
  glutCreateWindow(argv[0]);

  /* register callback functions */
  glutDisplayFunc(display);
  glutKeyboardFunc(keypress); //for all normal keys
  glutMouseFunc(mousepress);

  /* init GL */
  /* set background color black*/
  glClearColor(0.9, 0.9, 0.9, 0.5);
  /* here we can enable depth testing and double buffering and so
     on */

  /* give control to event handler */
  glutMainLoop();
  return 0;
}

//Might to try get polygons to be filled later using Jack's code 
//(right now it doesn't compile)


// // Called for errors in tessellation
// void polygon_error(GLenum errno) {
//   printf("ERROR!!!\n");
//   printf("The error is %s\n", gluErrorString(errno));
// }

// // Called for intersections in tessellation
// void combine(GLdouble coords[3], void* vertex_data[4], GLfloat weight[4],
//              void** outData) {
//   array<double, 3> point = {coords[0], coords[1], coords[2]};
//   *outData = &point[0];
// }

// // Draw the polygon!
// void draw_polygon(vector<point2D> poly) {
//   if (poly.size() == 0) return;

//   // convert the structs to arrays for the OpenGL C API
//   // std::vector<double*> points;
//   vector<array<double, 3>> points;
//   for (point2D p : poly) {
//     array<double, 3> point = {(double)p.x, (double)p.y, 0};
//     points.push_back(point);
//   }

//   // Tessellate and fill the polygon
//   auto tess = gluNewTess();
//   gluTessCallback(tess, GLU_BEGIN, (GLvoid (*) ())glBegin);
//   gluTessCallback(tess, GLU_VERTEX, (GLvoid (*) ())glVertex3dv);
//   gluTessCallback(tess, GLU_END, (GLvoid (*) ())glEnd);
//   gluTessCallback(tess, GLU_ERROR, (GLvoid (*) ())polygon_error);
//   gluTessCallback(tess, GLU_TESS_COMBINE, (GLvoid (*) ())combine);

//   gluBeginPolygon(tess);
//   for (size_t i = 0; i < points.size(); i++) {
//     gluTessVertex(tess, &points[i][0], &points[i][0]);
//   }
//   gluEndPolygon(tess);
// }

/* ****************************** */
/* Draws a polygon from a 1D vector */
/* has a special case that draws the polygon filled in if 'fill' is set to true */
void draw_polygon(vector<point2D> poly, bool fill){
    if (poly.size() == 0) return;

    //check if the this polygon should be filled (Doesn't work)
    if (fill) {
        //if so, set color to cyan and set transparency
        glColor4f(cyan[0], cyan[1], cyan[2], 0.1);

        //run through the polygon (visible polygon in this case)
        for (int i = 0; i < poly.size()-2; i++) {
            // //do wraparound for next point
            // int j = i + 1;
            // if (j == poly.size()) {
            //     j = 0;
            // }

            //draw a triangle from input point to the two neighboring points
            //in the vector
            glBegin(GL_TRIANGLES);
            glVertex2f(poly[i].x, poly[i].y);
            glVertex2f(poly[i+1].x, poly[i+1].y);
            glVertex2f(poly[i+2].x, poly[i+2].y);
            glEnd();
        }
        //render last 2 segments between last point and first point
        int last=poly.size()-1;
        glBegin(GL_TRIANGLES);
        glVertex2f(mouse_x, mouse_y);
        glVertex2f(poly[last].x, poly[last].y);
        glVertex2f(poly[0].x, poly[0].y);
        glEnd();

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
        //render last segment between last point and first point
        int last=poly.size()-1;
        glBegin(GL_LINES);
        glVertex2f(poly[last].x, poly[last].y);
        glVertex2f(poly[0].x, poly[0].y);
        glEnd();
    }
}

/* Draws all polygons in polys*/ 
void draw_polygons(){
  if (polys.size() == 0) return;

  for(int i = 0; i < polys.size(); i++){
    //If we try to fill the polygons, need to change Jack's code
    draw_polygon(polys[i], false);

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

  //Draw all polygons
  draw_polygon(polyTemp, false);

  //Draw start and end points if initialized
  if(startInit){
    glColor3fv(blue);
    draw_circle(start.x, start.y);
  }
  if(endInit){
    glColor3fv(blue);
    draw_circle(start.x, start.y);
    draw_circle(endPt.x, endPt.y);
  }

  // point2D point;
  // point.x = mouse_x;
  // point.y = mouse_y;

  //Not drawing
  if (!poly_init_mode) {

      //Draw all visible points here
      calculate_VG();
      printVG();
      makeAdjacency();
      printAdjList();
      runDijkstra();
  }

  draw_polygons();

  if (!poly_init_mode) {
      printShortestPath();

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
bool is_inside(point2D point, vector<point2D> poly) {
    //this is just to allow drawing when in drawing mode (affects only the call
    //to is_inside within the mouse press handler)
    if ((poly_init_mode == 1) || (poly.size() == 0)) {
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
bool simple(vector<point2D> poly){
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

            //draw intersection if they intersect
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
                isSimple = false; //not simple
            }
        }
    }
    return isSimple;
}

//draw visibility graph
void printVG(){
  for(int i=0; i< VG.size(); i++){
      glColor3fv(cyan);
      glBegin(GL_LINES);
      glVertex2f(VG[i].p1.x, VG[i].p1.y);
      glVertex2f(VG[i].p2.x, VG[i].p2.y);
      glEnd();
  }
}

//runs through each point in each polygon and figures out what it can see
void calculate_VG() {
    VG.clear();
    vector<point2D> pt1;
    pt1.push_back(start);
    vector<point2D> pt2;
    pt2.push_back(endPt);
    polys.push_back(pt1);
    polys.push_back(pt2);

    //add points from the bounding box to polys
    // addBoundingBox();

    //for each point in each polygon
    for (int polygon = 0; polygon < polys.size(); polygon++) {
        for (int point = 0; point < polys[polygon].size(); point++) {
            //check the lines that can be drawn from the point to other polygon
            //vertices
            check_all_points(polygon, point);

            //also check if there are any vertices on the specific polygon
            //where we can draw a line between
            check_local_points(polygon, point);
        }
    }

    //remove all four bounding box points and the two (start and end) points from polys
    //if boundingbox is added, for loop should go to 6, not 2. currently boudning box
    //is not added, so we just remove the end and start point
    for(int i = 0; i < 2; i++) {
        polys.pop_back();
    }
}

//add the bounding box to the polys vector
void addBoundingBox() {
    vector<point2D> box;
    point2D p;

    p.x = 0;
    p.y = 0;
    box.push_back(p);
    polys.push_back(box);
    box.clear();

    p.x = 0;
    p.y = WINDOWSIZE;
    box.push_back(p);
    polys.push_back(box);
    box.clear();

    p.x = WINDOWSIZE;
    p.y = 0;
    box.push_back(p);
    polys.push_back(box);
    box.clear();

    p.x = WINDOWSIZE;
    p.y = WINDOWSIZE;
    box.push_back(p);
    polys.push_back(box);
    box.clear();
}

//for each point in the polygon, check all points after that polygon
void check_all_points(int poly, int pt) {
    point2D point = polys[poly][pt];
    //extend a segment to each vertex after the current polygon
    //we can go from poly + 1 to the end because we don't need to
    //go back and re-check things that we've already checked
    for (int i = poly + 1; i < polys.size(); i++) {
        //if (i != polys.size()) {
            for (int j = 0; j < polys[i].size(); j++) {
                //make the segment
                lineSegment2D seg1;
                seg1.p1 = point;
                seg1.p2 = polys[i][j];

                //check if this segment intersects any of the segments of other polygons.
                //if it doesn't add it to the VG
                if(check_seg_intersections(seg1) == false) {
                    VG.push_back(seg1);
                }
            }
        //}
    }
}

//check if a particular segment intersects any other segments in the whole scene
bool check_seg_intersections(lineSegment2D seg1) {
    //go through each polygon
    for (int k1 = 0; k1 < polys.size(); k1++){
        //go through each segment in each polygon
        if (polys[k1].size() != 1) {
            for (int l1 = 0; l1 < polys[k1].size(); l1++) {
                lineSegment2D seg2;
                seg2.p1 = polys[k1][l1];

                //deal with wrap-around
                int l2 = l1 + 1;
                if (l2 == polys[k1].size()) {
                    l2 = 0;
                }
                seg2.p2 = polys[k1][l2];
                //Segments do not share a point
                if (check_seg_points_not_equal(seg1, seg2)){
                    //Segments intersect
                    if (intersect(seg1, seg2)) {
                        return true;
                    }
                }
            }
        }
    }
    //no intersections means add the segment to the VG
    return false;
}

//returns false if two segments ARE equal. returns true if they are NOT equal
bool check_seg_points_not_equal(lineSegment2D seg1, lineSegment2D seg2) {
    if (samePoint(seg1.p1, seg2.p1)) {
        return false;
    }
    if (samePoint(seg1.p1, seg2.p2)) {
        return false;
    }
    if (samePoint(seg1.p2, seg2.p1)) {
        return false;
    }
    if (samePoint(seg1.p2, seg2.p2)) {
        return false;
    }

    return true;
}

/**************************************************/
//check points within the polygon poly that are visible from a point pt
void check_local_points(int poly, int pt) {
    point2D point = polys[poly][pt];

    //Iterate through current polygon
    for (int i = 0; i < polys[poly].size(); i++) {
        //Get next and previous point indices
        int j = pt + 1;
        if (j == polys[poly].size()) {
            j = 0;
        }
        int h = pt - 1;
        if (h == -1) {
            h = polys[poly].size() - 1;
        }
        //Avoid checking adjacent points and current point with itself
        if ((i != pt) && (i != j) && (i != h)) {
            //Create segment from pt to a point
            lineSegment2D seg1;
            seg1.p1 = point;
            seg1.p2 = polys[poly][i];

            //Segment does not intersect with any other segments in the scene
            if (!check_seg_intersections(seg1)) {
                bool intersects = false;
                //Iterate through polygon again and check for adjacent segs
                for (int l1 = 0; l1 < polys[poly].size(); l1++) {

                    //Create segment from adjacent points in polygon
                    lineSegment2D seg2;
                    seg2.p1 = polys[poly][l1];
                    int l2 = l1 + 1;
                    if (l2 == polys[poly].size()) {
                        l2 = 0;
                    }
                    seg2.p2 = polys[poly][l2];

                    //None of the points on seg1 and seg2 are the same
                    if (check_seg_points_not_equal(seg1, seg2)){
                        //Segments don't intersect
                        if (intersect(seg1, seg2)) {
                            intersects = true;
                            break;
                        }
                    }
                }

                //No proper intersections found within polygon
                if (!intersects) {
                    //if the midpoint is inside the polygon, don't
                    //add the segment.
                    //if its not inside the polygon,
                    //add the segment
                    point2D midPoint;
                    midPoint.x = (seg1.p1.x + seg1.p2.x)/2;
                    midPoint.y = (seg1.p1.y + seg1.p2.y)/2;

                    if (!is_inside(midPoint, polys[poly])) {
                        VG.push_back(seg1);
                    }
                }
            }
        }
    }
    //Add adjacent points to the VG
    lineSegment2D seg;
    seg.p1 = polys[poly][pt];
    int next = pt + 1;
    if (next == polys[poly].size()) {
        next = 0;
    }
    seg.p2 = polys[poly][next];

    VG.push_back(seg);
}

//Create adjacency list for each point in VG
void makeAdjacency(){
    //Create two temporary vectors to hold points in adjacency list vectors
    vector<point2D> a;
    vector<point2D> b;

    adjList.clear();

    //Iterate through segments in VG
    for(int i=0; i < VG.size(); i++) {
        //extract the points from the segment
        lineSegment2D seg = VG[i];
        point2D p1 = seg.p1;
        point2D p2 = seg.p2;

        //initialize indices of each point to -1
        //this will be used if we find that the point is already in the adjList
        int index1 = -1;
        int index2 = -1;

        //go through the current adjList, if we find either of the points,
        //add the other point to that point's list
        for (int j = 0; j < adjList.size(); j++) {
            if (samePoint(p1, adjList[j][0])) {
                index1 = j; //indicate that we found it
                adjList[j].push_back(p2);
            }
            if (samePoint(p2, adjList[j][0])) {
                index2 = j;
                adjList[j].push_back(p1);
            }
        }
        //if we didn't find p1 in the list,
        //add it to the end of the list, and add p2 as an adjacent point
        if (index1 == -1) {
            a.push_back(p1);
            a.push_back(p2);
            adjList.push_back(a);
            a.clear();
        }

        //same but for p2
        if (index2 == -1) {
            b.push_back(p2);
            b.push_back(p1);
            adjList.push_back(b);
            b.clear();
        }
    }

    //set indices of all points in the adjList by iterating through the first column
    //and then setting all other instances of that point to the same index i
    //this comes in very very very useful in Dijkstra's
    for (int i = 0; i < adjList.size(); i++) {
        point2D p = adjList[i][0];
        for (int k = 0; k < adjList.size(); k++) {
            for (int j = 0; j < adjList[k].size(); j++) {
                if (samePoint(p, adjList[k][j])) {
                    adjList[k][j].index = i;
                }
            }
        }
    }
}

//Algorithm to run Dijkstra's to find Shortest Path from source to destination
//This algorith was NOT implemented with a priority queue, because we figured
//that the slowness of the visibility graph algorithm was bad enough that we
//didn't have to speed this algorithm up.
void runDijkstra() {
    shortestPath.clear();

    //initialize Q, weights, and previous' vector.
    vector<point2D> Q; //holds all points we can't "see", but will remove a point each iteration
    vector<double> weights; //keeps track of how far away each point is from the source
    vector<int> prevs; //keeps track of the previous point on the shortest path
    //initialize target and source points from globals start and endPt
    point2D source = start;
    point2D target = endPt;

    //build Q, weights, and prevs by taking first point in each adjacency list
    for (int i = 0; i < adjList.size(); i++) {
        //set weight to MAX_INT, unless we're at the source point. in that case set it to 0
        if (samePoint(source, adjList[i][0])) {
            weights.push_back(0);
        } else {
            weights.push_back(INT_MAX);
        }
        //initialize all previous' to -1
        prevs.push_back(-1);
        //add each point to Q
        Q.push_back(adjList[i][0]);
    }

    //make a copy of Q, to use later
    vector<point2D> QOriginal = Q;

    //iterate through Q until it is empty
    while (Q.size() != 0) {
        //find the index of the lowest weight
        int uIndex = get_index_min_distance(Q, weights);
        //save that point
        point2D u = Q[uIndex];
        //remove it from Q
        Q.erase(Q.begin() + uIndex);

        //if you've found the target, save it
        if (samePoint(u, target)) {
            target = u;
        }

        //grab the index that you will use to get into the adj list
        int adjIndex = u.index;

        //interate through u's adjacency list, skipping over u itself by starting at i=1
        for (int i = 1; i < adjList[adjIndex].size(); i++) {
            //calculate the distance by adding the weight of u to the euclidean
            //distance between u and it's adjacent point
            int distance = weights[u.index] + get_distance(u, adjList[adjIndex][i]);
            //if the calculated distance is lower, set is as the weight and set
            //u as the previous of its adjacent point
            if (distance < weights[adjList[adjIndex][i].index]) {
                weights[adjList[adjIndex][i].index] = distance;
                prevs[adjList[adjIndex][i].index] = u.index;
            }
        }
    }

    //now build the shortest path by going from the target and iterating back
    //until we hit the source.
    while (prevs[target.index] != -1) {
        lineSegment2D seg;
        seg.p2 = target;
        seg.p1 = QOriginal[prevs[target.index]];

        shortestPath.insert(shortestPath.begin(), seg);

        target = QOriginal[prevs[target.index]];
    }
}

//draw shortest path with a slightly larger width
void printShortestPath() {
    glLineWidth(3.0);
    for (int i = 0; i < shortestPath.size(); i++) {
        glColor3fv(green);
        glBegin(GL_LINES);
        glVertex2f(shortestPath[i].p1.x, shortestPath[i].p1.y);
        glVertex2f(shortestPath[i].p2.x, shortestPath[i].p2.y);
        glEnd();
    }
    glLineWidth(1.0);
}

//returns euclidean distance between two points
double get_distance(point2D a, point2D b) {
    return sqrt(pow((a.x - b.x),2) + pow((a.y - b.y),2));
}

//find the point in Q with the lowest weight
int get_index_min_distance(vector<point2D> Q, vector<double> weights) {
    int indexOfMin = -1;
    int min = INT_MAX;
    //iterate through all of Q
    for (int i = 0; i < Q.size(); i++) {
        //find the weight for Q by indexing into the weights vector using Q[i].index
        if (weights[Q[i].index] < min) {
            min = weights[Q[i].index];
            indexOfMin = i; //save index
        }
    }

    //return index of the lowest weight point
    return indexOfMin;
}

//print the adjList points and their relevant index
void printAdjList(){
    for(int i=0; i < adjList.size(); i++){
        for(int j=0; j < adjList[i].size(); j++){
            cout << "(" << adjList[i][j].x << ", " << adjList[i][j].y << "); ";
            cout << adjList[i][j].index << " ";
        }
        cout << endl;
    }
}

/* ****************************** */
void keypress(unsigned char key, int x, int y) {
    switch(key) {
        case 'q':
        exit(0);
        break;

    //expected behaviour: press 's', then click on the points you
    //want, and press 'j' when you're done with each polygon and 'e' when you're done with all polygons. the points will be saved
    //in global 'poly'

    case 's':
        //start re-initializeing polygons
        //re-set all vectors
        polys.clear();
        polyTemp.clear();
        VG.clear();
        //re-set the position of the mouse
        mouse_x = mouse_y=0;
        //change to drawing mode
        poly_init_mode = 1;
        startInit = 0;
        endInit = 0;
        glutPostRedisplay();
        break;

    case 'f':
        if(!simple(polyTemp)){
          polyTemp.clear();
          cout << "Polygon is not simple. Please erase, ('r') and draw another one." << endl;
          break;
        }
        //Add polygon to vector of polygons and clear temporary polygon
        polys.push_back(polyTemp);
        polyTemp.clear();
        poly_init_mode=1;
        glutPostRedisplay();
        break;

    //User wants to erase current polygon and start over
    case 'r':
        polyTemp.clear();
        poly_init_mode=1;
        glutPostRedisplay();
        break;

    //Initial beginning point
    case 'b':
        cout << "You are done entering polygons. Now enter start point." << endl;
        polys.push_back(polyTemp);
        startInit = 1;
        endInit = 0;
        poly_init_mode=1;
        glutPostRedisplay();
        break;

    //Initial beginning point
    case 'e':
        cout << "Now enter end point." << endl;
        startInit = 0;
        endInit = 1;
        poly_init_mode=1;
        glutPostRedisplay();
        break;
    //Calculate and show VG
    case 'v':
        cout << "Calculating VG" << endl;
        poly_init_mode = 0;
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
