#include "kdtree.h"

/* Libraries */
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <vector>

/* OpenGL */
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

void display();
void keypress(unsigned char key, int x, int y);
void initialize_points_random();
void initialize_points_diag();
void initialize_points_grid();
void switch_cases();
void renderTree(TreeNode* node, int xMin, int xMax, int yMin, int yMax);
void draw_points();
void draw_point(point2D p);
void reshape(GLsizei width, GLsizei height);

const int WINDOWSIZE = 500;

int n;
vector<point2D> p;
Kdtree* k;
GLfloat *colors[5] = {red, yellow, blue, white, white};
int init_case = 0;
const int NB_TEST_CASES = 3;

int main(int argc, char** argv) {
    //seed random generator
    srand(time(NULL));

    //parse command line arguments - make sure user enters 2 args
    if (argc == 2) {
        n = atoi(argv[1]);
    } else {
        cout << "Incorrect usage. Please enter the number of points you wish to include" << endl;
        cout << "e.g.: ./main 50 for 50 points" << endl;
        exit(0);
    }
    cout << "You entered n = " << n << endl;
    // initialize_points_random();
    initialize_points_diag();
    
    // initialize_points_grid();

    // Make a new Kdtree with n points, p as vector of points, and first-cut = vert;
    // in this program 1 is vertical cut, and 0 is horizontal cut
    k = new Kdtree(p, n, 1);

    //Print Tree
    k->printTree(k->root);

    //Get height of tree and number of nodes
    int treeHeight = k->getHeight(k->root);
    cout << "Height of KD Tree is: " << treeHeight-1 << endl; //Height is actually 1 less

    //The number of nodes counted includes
    cout << "Number of nodes in KD Tree is: " << k->numNodes << endl;


    /*********************************************/
    /* The following is copied from Laura's code */
    /* initialize GLUT  */
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
	glutInitWindowPosition(100,100);
	glutCreateWindow(argv[0]);

	/* register callback functions */
	glutDisplayFunc(display);
	glutKeyboardFunc(keypress);

	/* init GL */
	/* set background color white*/
	glClearColor(1, 1, 1, 1);


	/* give control to event handler */
	glutMainLoop();
	return 0;
}

// Simple point initialization made by Bo and Henry
void initialize_points_random() {
    p.clear();

    point2D a;
    for (int i = 0; i < n; i++) {
        a.x = (int)(.3*WINDOWSIZE)/2 + rand() % ((int)(.7*WINDOWSIZE));
        a.y =  (int)(.3*WINDOWSIZE)/2 + rand() % ((int)(.7*WINDOWSIZE));
        p.push_back(a);
    }
}

// Bo Bleckel & Henry Daniels Koch
// Initialize the points in a downward sloping diagonal
// works best with n <= 100
void initialize_points_diag() {
    p.clear();

    point2D point;
    for (int i = 0; i < n; i++) {
        double x = i * (WINDOWSIZE/n);
        double y = (n-i) * (WINDOWSIZE/n);
        point.x = x;
        point.y = y;

        p.push_back(point);
    }
}

// Jack Ward
// Initialize the `points` vector to a grid.
void initialize_points_grid() {

    // clear the vector
    p.clear();

    double window = (double) WINDOWSIZE;
    double padding = window / 8;
    double width = window - 2 * padding;
    int side = sqrt(n);
    double spacing = width / (double) side;

    for (int row = 0; row < side; row++) {
        for (int col = 0; col < side; col++) {
                     // col * spacing + padding
          double x = fma(col, spacing, padding);
                     // row * spacing + padding
          double y = fma(row, spacing, padding);

          point2D point;
          point.x = x;
          point.y = y;
          p.push_back(point);
        }
    }
}

// Render the Kdtree. This function is recursive, and goes through the whole tree
void renderTree(TreeNode* node, int xMin, int xMax, int yMin, int yMax) {
    // null leaf node
    // should never go here, but just in case
    if(node == NULL) {
        return;
    }

    // leaf node
    if(node->type == 'l') {
        // make a background rectangle that is black (to use as frame)
        // and fill the node's window with a colored rectangle
        // on top of the background

        //offset for background rect
        int width = 4;
        glColor3fv(black);
        glBegin(GL_POLYGON);
        glVertex2f(xMin-width, yMin-width);
        glVertex2f(xMin-width, yMax+width);
        glVertex2f(xMax+width, yMax+width);
        glVertex2f(xMax+width, yMin-width);
        glEnd();

        //choose a random color from the random color list
        //and make a rectangle of that color in the window-frame
        int randNum = rand()%5;
        GLfloat *randColor = colors[randNum];
        glColor3fv(randColor);
        glBegin(GL_POLYGON);
        glVertex2f(xMin, yMin);
        glVertex2f(xMin, yMax);
        glVertex2f(xMax, yMax);
        glVertex2f(xMax, yMin);
        glEnd();

        return;
    }

    //vertical
    if(node->type == 'v') {
        //identify the endpoints p1 and p2 of the line segment that you
        //need to draw
        int xVal = node->p.x;

        //make recursive call to left and right
        renderTree(node->left, xMin, xVal, yMin, yMax);
        renderTree(node->right, xVal, xMax, yMin, yMax);
    }

    //horizontal
    if(node->type == 'h') {
        //identify the endpoints p1 and p2 of the line segment that you
        //need to draw
        int yVal = node->p.y;

        // make recursive call to left and right
        renderTree(node->left, xMin, xMax, yMin, yVal);
        renderTree(node->right, xMin, xMax, yVal, yMax);
    }
}

/* ****************************** */
/* Copied from Laura's code */
void display() {

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

    //render the tree
    renderTree(k->root, 0, WINDOWSIZE, 0, WINDOWSIZE);

    /** uncomment below if you want to see points **/
    // draw_points();

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
            // go to next case of points
            switch_cases();
            // re-build Kdtree with new points
            k = new Kdtree(p, n, 1);

            //re-render display
			glutPostRedisplay();
			break;
        case '1':
            initialize_points_random();
            // re-build Kdtree with new points
            k = new Kdtree(p, n, 1);

            //re-render display
			glutPostRedisplay();
			break;
        case '2':
            initialize_points_grid();
            // re-build Kdtree with new points
            k = new Kdtree(p, n, 1);

            //re-render display
            glutPostRedisplay();
            break;
        case '3':
            initialize_points_diag();
            // re-build Kdtree with new points
            k = new Kdtree(p, n, 1);

            //re-render display
			glutPostRedisplay();
			break;
	}
}

/** Adapted directly from Laura's code on Orthogonal Line intersection **/
/* Add more test cases as needed here and update NB_TEST_CASES at the */
/* top of the file */
void switch_cases() {
    switch (init_case)  {
		case 0:
			initialize_points_random();
			break;
		case 1:
			initialize_points_diag();
			break;
		case 2:
			initialize_points_grid();
			break;
		default:
			initialize_points_random();
	}

    init_case = (init_case+1) % NB_TEST_CASES;
	return;
}

/* To draw all points */
/** Adapted from Laura's code **/
void draw_points() {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    //set color
    glColor3fv(cyan);
    for (int i = 0; i < n; i++) {
        //draw a small square centered at (points[i].x, points[i].y)
        draw_point(p[i]);
    }
}

/* To draw one point */
void draw_point(point2D p) {

    const int R = 2;
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    //set color
    glColor3fv(green);

    //draw a small square centered at (points[i].x, points[i].y)
    glBegin(GL_POLYGON);
    glVertex2f(p.x - R, p.y - R);
    glVertex2f(p.x + R, p.y - R);
    glVertex2f(p.x + R, p.y + R);
    glVertex2f(p.x - R, p.y + R);
    glEnd();

}

/** Taken directly from Laura's code **/
/* Handler for window re-size event. Called back when the window first appears and
 whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer

	// Set the viewport to cover the new window
	glViewport(0, 0, width, height);

	glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
	glLoadIdentity();             // Reset
	gluOrtho2D(0.0, (GLdouble) width, 0.0, (GLdouble) height);
}
