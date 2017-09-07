/* main.cpp

Henry and Bo

A* Motion Planning and RRT Motion Planning
*/
#include <omp.h>
#include "main.h"
#include "geom.h"

const static int thread_num = 16;

using namespace std;

/*
    Determine if two states are equal
*/
bool equalStates(state s1, state s2) {
    return ((s1.x == s2.x) && (s1.y == s2.y) && (s1.theta == s2.theta));
}

/*
    Governs the running of the A* algorithm.
    Returns the final path of states
*/
vector<state> aStar() {
    makeBoundBoxes();
    clock_t startT, endT;
    double totalT;
    startT = clock();

    if(firstTime){
	    allocateGrid();
	    if (preProcess)
	        preProcessGrid();
	    endT = clock();
	    totalT = ((double) (endT - startT)) / CLOCKS_PER_SEC;
	    cout << "Grid is ready (took " << totalT << " seconds to build)." << endl;
	} else {
		resetGrid();
	}
    startT = clock();

    state startState;
    startState.x = start.x;
    startState.y = start.y;
    startState.theta = 0;
    startState.weight = 0.0;
    // stateGrid[startState.x][startState.y][startState.theta].prevX = -1;
    // stateGrid[startState.x][startState.y][startState.theta].prevY = -1;
    stateGrid[startState.x][startState.y][startState.theta].weight = 0.0;
    stateGrid[startState.x][startState.y][startState.theta].explored = true;

    pq.push(startState);

    state finalState;
    finalState.x = endPt.x;
    finalState.y = endPt.y;
    finalState.theta = 0;


    bool goalFound = false;
    int numExplored = 0;

    cout <<  "Solving..." << endl;
    while (!goalFound) {
        endT = clock();
        totalT = ((double) (endT - startT)) / CLOCKS_PER_SEC;
        cout << "\r" << totalT;
        state s = pq.top();
        pq.pop();


        stateGrid[s.x][s.y][s.theta].explored = true;
        numExplored++;
        // state temp;
        // temp.x = s.prevX;
        // temp.y = s.prevY;
        // temp.theta = s.prevTheta;
        // setPrevs(s, temp);

        vector<state> successors = generateSuccessors(s);


        // cout << "successor size: " << successors.size() << endl;

        for (int i = 0; i < successors.size(); i++) {
            // cout << successors[i].x << ", " << successors[i].y << ", " << successors[i].theta << endl;
            state suc = successors[i];
            //Sucessor state is the final state
            if (equalStates(suc, finalState)) {
                if (stateGrid[suc.x][suc.y][suc.theta].occupied) {
                    goalFound = true;
                    finalState = suc;
                    setPrevs(suc, s);
                    cout << "\rFound final state." << endl;
                    break; //done searching
                }
            }
            //Successor state has not been explored
            if (!stateGrid[suc.x][suc.y][suc.theta].explored) {
            	//Set state to having been explored
                suc.prevWeight = stateGrid[s.x][s.y][s.theta].weight;
                suc.weight = scoreState(suc, HEUR_WEIGHT);
                if (preProcess) {
                    if (stateGrid[suc.x][suc.y][suc.theta].occupied && suc.weight < stateGrid[suc.x][suc.y][suc.theta].weight) {
                        setPrevs(suc, s);
                        stateGrid[suc.x][suc.y][suc.theta].weight = suc.weight;
                        pq.push(suc);
                    }
                } else {
                    if ((stateGrid[suc.x][suc.y][suc.theta].occupied || isFree(suc)) && suc.weight < stateGrid[suc.x][suc.y][suc.theta].weight) {
                        stateGrid[suc.x][suc.y][suc.theta].occupied = true;
                        setPrevs(suc, s);
                        stateGrid[suc.x][suc.y][suc.theta].weight = suc.weight;
                        pq.push(suc);
                    }
                }
            }
        }
        // exit(0);
    }
    endT = clock();
    totalT = ((double) (endT - startT)) / CLOCKS_PER_SEC;
    cout << "Explored " << numExplored << " states. Total time spent looking: " << totalT << " seconds." << endl;

    vector<state> path;
    while (finalState.x != -1) {
        path.insert(path.begin(), finalState);
        // draw_circle(finalState.x, finalState.y);
        finalState = stateGrid[finalState.x][finalState.y][finalState.theta].prev;
    }

    path.insert(path.begin(), startState);

    return path;
}

void setPrevs(state suc, state s) {
    // stateGrid[suc.x][suc.y][suc.theta].prevX = s.x;
    // stateGrid[suc.x][suc.y][suc.theta].prevY = s.y;
    // stateGrid[suc.x][suc.y][suc.theta].prevTheta = s.theta;
    stateGrid[suc.x][suc.y][suc.theta].prev = s;
}

/*
    Compute a bool that determines if a configuration of a robot is possible
*/
bool isFree(state s) {
    //extract point
    point2D p;
    p.x = s.x;
    p.y = s.y;

    vector<point2D> robotPoly = generateRobot(s);

    //Check if robot first intersects any bounding box edge
    if (!boundBoxIntersect(robotPoly)) {
        return true;
    }

    //Parallelize Check
    #pragma omp parallel num_threads(thread_num)
	{
	#pragma omp for collapse(3) private(r,r2,,i,j,j2)

    for (int r = 0; r < robotPoly.size(); r++) { //go through robot vertices
        lineSegment2D robotEdge;
        int r2 = r + 1;
        if (r2 == robotPoly.size()) {
            r2 = 0;
        }
        //Create edge of robot
        robotEdge.p1 = robotPoly[r];
        robotEdge.p2 = robotPoly[r2];

        //Check robot edge with every polygon edge
        for (int i = 0; i < polys.size(); i++) { //go through all polygons
            for (int j = 0; j < polys[i].size(); j++) { //go through each polygon's edges
                int j2 = j + 1;
                if (j2 == polys[i].size()) {
                    j2 = 0;
                }
                lineSegment2D polyEdge;
                polyEdge.p1 = polys[i][j];
                polyEdge.p2 = polys[i][j2];

                if (intersect(robotEdge, polyEdge)) {
                    return false;
                }
            }
        }
    }
	}

	//Parallelize
	#pragma omp parallel num_threads(thread_num)
	{
	#pragma omp for collapse(1) private(i)

    //check if the point is inside any polygons
    for (int i = 0; i < polys.size(); i++) {
        if (is_inside(p, polys[i])) {
            return false;
        }
    }
	}

    return true;
}

//Generate bounding boxes
void makeBoundBoxes(){

    //Set up box
	int xMin = INT_MAX;
    int xMax = INT_MIN;
    int yMin = INT_MAX;
    int yMax = INT_MIN;

	for (int i = 0; i < polys.size(); i++) { //go through all polygons
	    for (int j = 0; j < polys[i].size(); j++) { //go through each polygon's vertices

	        if(polys[i][j].x < xMin)
	        	xMin = polys[i][j].x;
	        if(polys[i][j].x > xMax)
	        	xMax = polys[i][j].x;
	        if(polys[i][j].y < yMin)
	        	yMin = polys[i][j].y;
	        if(polys[i][j].y > yMax)
	        	yMax = polys[i][j].y;
	    }

        boundingBox box;
        box.xMin = xMin;
        box.xMax = xMax;
        box.yMin = yMin;
        box.yMax = yMax;

	    boundBoxes.push_back(box);
	}
	return;

}

//Determines if edge of robot intersects any bounding box of a polygon
bool boundBoxIntersect(vector<point2D> robot){
    //Set up box
	int xMin = INT_MAX;
    int xMax = INT_MIN;
    int yMin = INT_MAX;
    int yMax = INT_MIN;

    for (int i = 0; i < robot.size(); i++) { //go through each polygon's vertices
        if(robot[i].x < xMin)
            xMin = robot[i].x;
        if(robot[i].x > xMax)
            xMax = robot[i].x;
        if(robot[i].y < yMin)
            yMin = robot[i].y;
        if(robot[i].y > yMax)
            yMax = robot[i].y;
    }

    #pragma omp parallel num_threads(thread_num)
	{
	#pragma omp for collapse(1) private(i)

	//Check all polygons
	for(int i = 0; i < boundBoxes.size(); i++){

        //Check if nothing touches. If anything does, return true (found an intersection
        //with the bouding box)
        if (xMin > boundBoxes[i].xMax || xMax < boundBoxes[i].xMin || yMin > boundBoxes[i].yMax || yMax < boundBoxes[i].yMin) {

        } else {
            return true;
        }
	}
	}
	return false;
}

vector<point2D> generateRobot(state s) {
    vector<point2D> robotPoly;

    int shortSide = robot.shortSide;
    int longSide = robot.longSide;
    int theta = s.theta;

    point2D p;
    p.x = s.x;
    p.y = s.y;
    robotPoly.push_back(p);

    p.x = p.x + shortSide * cosVect[theta];
    p.y = p.y + shortSide * sinVect[theta];
    robotPoly.push_back(p);

    p.x = p.x + longSide * sinVect[theta];
    p.y = p.y - longSide * cosVect[theta];
    robotPoly.push_back(p);

    p.x = p.x - shortSide * cosVect[theta];
    p.y = p.y - shortSide * sinVect[theta];
    robotPoly.push_back(p);

    return robotPoly;
}

/*
    Generate all successors for a state s, even if we've already been there and no
    matter if it isn't free
*/
vector<state> generateSuccessors(state s) {
    //Generate vector of successors
    vector<state> successors;

    //Iterate through ALL 4 (maybe 5) translational options (maybe rotations too)
    for(int i= -dX; i <=dX; i+=dX){
        for(int j= -dY; j<=dY; j+=dY){

            //Skip Diagonals
            if(((i==-dX)&&(j==-dY)) || ((i==-dX)&&(j==dY))
                || ((i==dX)&&(j==-dY)) || ((i==dX)&&(j==dY))){
                continue;
            }

            //Robot can rotate
            if(rotations){
                //Rotations
                if (i == 0 || j == 0) {
                    for(int k= dtheta; k>=-dtheta; k-=dtheta){
                        //Skip state with no rotations or translation
                        if((i==0)&&(j==0)&&(k==0)) {
                            continue;
                        } else if (s.x + i >= WINDOWSIZE || s.y + j >= WINDOWSIZE /*|| s.theta + k >= 360*/
                            || s.x + i <= -1 || s.y + j <= -1 /*|| s.theta + k < 0*/) {
                            continue;
                        } else{
                            //Generate 3D state
                            state newState;
                            newState.x = s.x + i;
                            newState.y = s.y + j;
                            if (s.theta + k >= 360) {
                                newState.theta = 0;
                            } else if (s.theta + k < 0) {
                                newState.theta = 359;
                            } else {
                                newState.theta = s.theta + k;
                            }
                            successors.push_back(newState);
                        }
                    }
                } else {
                    //Skip state with no traslation
                    if((i==0)&&(j==0)) {
                            continue;
                    } else if (s.x + i >= WINDOWSIZE || s.y + j >= WINDOWSIZE
                        || s.x + i <= -1 || s.y + j <= -1) {
                        continue;
                    } else {
                        //Generate 3D state
                        state newState;
                        newState.x = s.x + i;
                        newState.y = s.y + j;
                        newState.theta = 0;
                        successors.push_back(newState);
                    }
                }
            //Robot CAN'T rotate
            }else{
                //Skip state with no traslation
                if (s.x + i >= WINDOWSIZE || s.y + j >= WINDOWSIZE
                    || s.x + i <= -1 || s.y + j <= -1) {
                    continue;
                } else {
                    //Generate 3D state
                    state newState;
                    newState.x = s.x + i;
                    newState.y = s.y + j;
                    newState.theta = 0;
                    successors.push_back(newState);
                }
            }
        }
    }
    return successors;
}

//Returns the Euclidian distance from state point to end point
double heuristicEuc(state s){
    int theta;
    if (s.theta > 180) {
        theta = 360 - s.theta;
    } else {
        theta = s.theta;
    }
    return sqrt( (s.x - endPt.x)*(s.x - endPt.x) + (s.y - endPt.y)*(s.y - endPt.y) + (theta * theta));
}

/*
    Returns the cost of a state by adding the previous cost from the start to the
    state, and then adding the weighted heuristic cost to the end point
*/
double scoreState(state s, double heurWeight) {
    double weight = s.prevWeight + 1 + heurWeight * heuristicEuc(s);
    return weight;
}

/*
	Reset parts of the grid. Only gets called if the isFree() part has already
	been computed.
*/
void resetGrid() {
	for (int x = 0; x < WINDOWSIZE; x++) {
        //vector<vector<stateProps> > vect1;
        for (int y = 0; y < WINDOWSIZE; y++) {
            //vector<stateProps> vect2;
            int angles;
            if (!rotations) {
                angles = 1;
            } else {
                angles = 360;
            }
            for (int theta = 0; theta < angles; theta += dtheta) {
                state g;
                g.x = -1;
                g.y = -1;
                g.theta = -1;

                stateGrid[x][y][theta].prev = g;
                stateGrid[x][y][theta].weight = INT_MAX;
                stateGrid[x][y][theta].explored = false;
            }
        }
    }
}

/*
    Allocate the grid if not pre processing that way the "explored" property
    can be used without a segfault
*/
void allocateGrid() {
    stateGrid.clear();

    cout << "Allocating grid..." << endl;
    for (int x = 0; x < WINDOWSIZE; x++) {
        vector<vector<stateProps> > vect1;
        for (int y = 0; y < WINDOWSIZE; y++) {
            vector<stateProps> vect2;
            int angles;
            if (!rotations) {
                angles = 1;
            } else {
                angles = 360;
            }
            for (int theta = 0; theta < angles; theta += dtheta) {
                stateProps prop;
                state g;
                g.x = -1;
                g.y = -1;
                g.theta = -1;

                prop.prev = g;
                prop.weight = INT_MAX;
                prop.explored = false;
                prop.occupied = false;
                vect2.push_back(prop);
            }
            vect1.push_back(vect2);
        }
        stateGrid.push_back(vect1);
    }
}

/*
    Pre process the grid in all three dimensions to determine if a state is
    free. This does not have to happen, it can happen on the fly if you want.
*/
void preProcessGrid() {
    //stateGrid.clear();

    #pragma omp parallel num_threads(thread_num)
	{
	#pragma omp for collapse(3) private(x,y,theta,angles)

    cout << "Pre processing grid..." << endl;
    for (int x = 0; x < WINDOWSIZE; x++) {
        //vector<vector<stateProps> > vect1;
        for (int y = 0; y < WINDOWSIZE; y++) {
            //vector<stateProps> vect2;
            int angles;
            if (!rotations) {
                angles = 1;
            } else {
                angles = 360;
            }
            for(int theta = 0; theta < angles; theta += dtheta){
                stateProps prop;
                state s;
                s.x = x;
                s.y = y;
                s.theta = theta;
                if (isFree(s)) {
                    prop.occupied = true;
                } else {
                    prop.occupied = false;
                }

                state g;
                g.x = -1;
                g.y = -1;
                g.theta = -1;

                prop.prev = g;
                prop.weight = INT_MAX;
                prop.explored = false;
                stateGrid[x][y][theta] = prop;
                //vect2.push_back(prop);
            }
            //vect1.push_back(vect2);
        }
        //stateGrid.push_back(vect1);
    }
	}
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
void mousepress(int button, int status, int x, int y) {

    if(status == GLUT_DOWN) {
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
            state s;
            s.x = mouse_x;
            s.y = mouse_y;
            s.theta = 0;
            if (isFree(s)) {
                astar_mode = 0;
                move_robot = 0;
                endPt.x = mouse_x;
                endPt.y = mouse_y;
                cout << "End: " << endPt.x << " " << endPt.y << endl;
            } else {
                cout << "Point (" << mouse_x << ", " << mouse_y << ") is not in free space. Please try again." << endl;
            }

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
        state s;
        s.x = start.x;
        s.y = start.y;
        s.theta = 0;
        draw_robot(s);
    }
    if(endInit){
    	//Why do we draw the robot at the end/at 45degrees?
        glColor3fv(blue);
        draw_circle(start.x, start.y);
        draw_circle(endPt.x, endPt.y);
        if (!move_robot) {
            state s;
            s.x = start.x;
            s.y = start.y;
            s.theta = 0;
            draw_robot(s);
        }
    }

    // point2D point;
    // point.x = mouse_x;
    // point.y = mouse_y;

    //Not drawing
    if (!poly_init_mode && astar_mode) {
        //Draw all visible points here
        drawMotionPath();
    }

    if (move_robot) {
        // cout << robotPos << endl;
        draw_robot(shortestPath[robotPos]);
    }

    draw_polygons();

    /* execute the drawing commands */
    glFlush();
}

void initializeRobot() {
    robot.shortSide = 20;
    robot.longSide = 40;
}

/* Draws the path that the robot will follow */
void drawMotionPath() {
    for (int i = 0; i < shortestPath.size() - 1; i++) {
        int j = i + 1;

        glColor3fv(red);
        glBegin(GL_LINES);
        glVertex2f(shortestPath[i].x, shortestPath[i].y);
        glVertex2f(shortestPath[j].x, shortestPath[j].y);
        glEnd();
    }
}

/* Draws the robot */
void draw_robot(state s) {
    vector<point2D> robotPoly = generateRobot(s);

    glColor3fv(white);
    glBegin(GL_POLYGON);
    for (int i = 0; i < robotPoly.size(); i++) {
        glVertex2f(robotPoly[i].x, robotPoly[i].y);
    }
    glEnd();
}

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

void initializeCosSin() {

    for (int i = 0; i < 360; i++) {
        cosVect.push_back(cos(i*2*M_PI/360.0));
        sinVect.push_back(sin(i*2*M_PI/360.0));
    }
}

int main(int argc, char** argv) {

	#pragma omp parallel
  	{


    //int thread_Num = omp_get_num_threads();
    cout << "Num threads " << thread_num << endl;
	}


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
    glutIdleFunc(timerfunc);

    initializeCosSin();
    initializeRobot();
    start.x = WINDOWSIZE/2;
    start.y = WINDOWSIZE/2;
    endPt.x = 0;
    endPt.y = 0;

    /* init GL */
    /* set background color black*/
    glClearColor(0.1, 0.1, 0.1, 0.0);
    /* here we can enable depth testing and double buffering and so
       on */

    /* give control to event handler */

    glutMainLoop();
    return 0;
}


//Using the ray-casting algorithm, determine if a point is inside the polygon
//This technique "casts" a ray from the point, to the right, and then counts
//the intersections of that ray with the polygon. If the number of intersections
//is even, the point is outside the polygon. if its odd, the point is inside.
bool is_inside(point2D point, vector<point2D> poly) {
    //this is just to allow drawing when in drawing mode (affects only the call
    //to is_inside within the mouse press handler)
    if (poly_init_mode && !(startInit || endInit)) {
        return true;
    }
    if (poly.size() == 0) {
        return false;
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
bool simple(vector<point2D> poly) {
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
        //re-set the position of the mouse
        mouse_x = mouse_y=0;
        firstTime = 1;
        //change to drawing mode
        poly_init_mode = 1;
        startInit = 0;
        endInit = 0;
        astar_mode = 0;
        move_robot = 0;
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
        astar_mode = 0;
        move_robot = 0;
        glutPostRedisplay();
        break;

    //User wants to erase current polygon and start over
    case 'r':
        polyTemp.clear();
        poly_init_mode=1;
        astar_mode = 0;
        move_robot = 0;
        glutPostRedisplay();
        break;

    //Initial beginning point
    case 'b':
        cout << "You are done entering polygons. Now enter start point." << endl;
        polys.push_back(polyTemp);
        startInit = 1;
        endInit = 0;
        poly_init_mode=1;
        astar_mode = 0;
        move_robot = 0;
        glutPostRedisplay();
        break;

    //Initial beginning point
    case 'e':
        cout << "Now enter end point." << endl;
        startInit = 0;
        endInit = 1;
        poly_init_mode=1;
        astar_mode = 0;
        move_robot = 0;
        glutPostRedisplay();
        break;
    //Calculate and show VG
    case 'v':
        poly_init_mode = 0;
        astar_mode = 1;
        move_robot = 0;
        //aStar();
        shortestPath = aStar();
        firstTime = 0;
        glutPostRedisplay();
        break;
    case 'a':
        move_robot = !move_robot;
        robotPos = 0;
        glutPostRedisplay();
    }
}

void timerfunc() {
    if (move_robot && astar_mode) {
        if (robotPos + 1 < shortestPath.size()) {
            robotPos++;
        }

        glutPostRedisplay();
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
