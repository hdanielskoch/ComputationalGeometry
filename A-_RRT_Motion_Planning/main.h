#ifndef MAIN_H
#define MAIN_H

#include "geom.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <ctime>
#include <climits>
#include <cmath>
#include <fstream>
#include <vector>
#include <queue>
#include "rtimer.h"
#include <assert.h>
#include <cfloat>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif



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

typedef struct _state {
	int x,y, theta;
	double weight;
	double prevWeight;
} state;

typedef struct _stateProps {
	bool explored, occupied;
	state prev;
	double weight;
} stateProps;

typedef struct _robotStruct {
	int shortSide, longSide;
	state currentState;
} robotStruct;

typedef struct _boundingBox {
	int xMin, xMax, yMin, yMax;
} boundingBox;

//adapted from:
//http://stackoverflow.com/questions/16111337/declaring-a-priority-queue-in-c-with-a-custom-comparator
class Compare {
public:
	bool operator() (state s1, state s2) {
		return s1.weight > s2.weight;
	}
};

//Functions
vector<state> aStar();
void setPrevs(state suc, state s);
//Determines if a state is IN the configuration space (doesn't intersect obstacles)
bool isFree(state s);
//Make bounding boxes for each polygon
void makeBoundBoxes();
//Determine if an edge intersects any polygon's bounding box
bool boundBoxIntersect(vector<point2D> robot);
//Generate all successor states
vector<state> generateSuccessors(state s);
//Gets euclidian distance to end point
double heuristicEuc(state s);
//scoes a state
double scoreState(state s, double heurWeight);
//pre processes the grid to hold each the value of each state (free or not)
void resetGrid();
void preProcessGrid();
//allocate space for the grid if not preprocessing, that way the explored can be accessed
void allocateGrid();
//builds a polygon for the robot given a position
vector<point2D> generateRobot(state s);

void display(void);
void initializeRobot();
void initializeCosSin();
bool is_inside(point2D point, vector<point2D> poly);
bool simple(vector<point2D> poly);

//Drawing functions
void mousepress(int button, int state, int x, int y);
void drawMotionPath();
void draw_robot(state s);
void draw_polygon(vector<point2D> poly, bool fill);
void draw_polygons();
void draw_circle(double x, double y);
void keypress(unsigned char key, int x, int y);
void timerfunc();

//Global variables
const int WINDOWSIZE = 500;

//coordinates of last mouse click
double mouse_x=-10, mouse_y=-10;


//when this is 1, then clicking the mouse results in those points being stored in poly
int poly_init_mode = 1;
int astar_mode = 0;
int move_robot = 0;
int firstTime = 1;

//When this is one, the user is done drawing polygones
int startInit = 0;
int endInit = 0;

point2D start;
point2D endPt;

//Queue to store all potential successor states with weights(doubles)
priority_queue<state, vector<state>, Compare> pq;

//vectors to store values for cos and sin to speed up processing
vector<double> cosVect;
vector<double> sinVect;

//Stores all polygons, temporary polygon, bounding box, VG, adjacency list, and shortest path
vector<vector<point2D> > polys;
vector<boundingBox> boundBoxes;
vector<point2D> polyTemp;
vector<state> shortestPath;

int robotPos; // governs the index in the shortest path vector that the robot is currently at
robotStruct robot;

double HEUR_WEIGHT = 1.0;
bool rotations = true; //Determines whether or not we are allowing the robot to rotate
bool preProcess = false;
int dtheta = 1; //Angle increment
int dX = 1; //X resolution
int dY = 1; //Y resolution
vector<vector<vector<stateProps> > > stateGrid; //3D vector contains properties of ALL states

#endif
