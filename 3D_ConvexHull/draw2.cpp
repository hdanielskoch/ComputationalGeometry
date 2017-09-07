/*

  Drawing exercise. 

  Allows to translate/rotate when user presses l/r/u/d/x/X,y/Y,z/Z).
  
  OpenGL 1.x
  Laura Toma
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
//this allows this code to compile both on apple and linux platforms
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#include <vector>

using namespace std; 



/* global variables */


const int WINDOWSIZE = 500; 

//we predefine some colors for convenience 
GLfloat red[3] = {1.0, 0.0, 0.0};
GLfloat green[3] = {0.0, 1.0, 0.0};
GLfloat blue[3] = {0.0, 0.0, 1.0};
GLfloat black[3] = {0.0, 0.0, 0.0};
GLfloat white[3] = {1.0, 1.0, 1.0};
GLfloat gray[3] = {0.5, 0.5, 0.5};
GLfloat yellow[3] = {1.0, 1.0, 0.0};
GLfloat magenta[3] = {1.0, 0.0, 1.0};
GLfloat cyan[3] = {0.0, 1.0, 1.0};


//keep track of global translation and rotation 
GLfloat pos[3] = {0,0,0};
GLfloat theta[3] = {0,0,0};


/* forward declarations of functions */
void display(void);
void keypress(unsigned char key, int x, int y);


void cube(GLfloat side, int mode); 
void draw_sphere(float rad);
void draw_my_scene(); 
void draw_circle(float rad); 
void draw_flower(); 


int main(int argc, char** argv) {

  

    /* open a window and initialize GLUT stuff */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
  glutInitWindowPosition(100,100);
  glutCreateWindow(argv[0]);

  /* register callback functions */
  glutDisplayFunc(display); 
  glutKeyboardFunc(keypress);
  
  /* OpenGL init */
  /* set background color black*/
  glClearColor(0, 0, 0, 0);  
    glEnable(GL_DEPTH_TEST); 

  /* setup the camera (i.e. the projection transformation) */ 
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60, 1 /* aspect */, 1, 10.0); /* the frustrum is from z=-1 to z=-10 */
  /* camera is at (0,0,0) looking along negative y axis */
  
  //the points are in z=[-1,1].  initialize the translation to bring
  //the points in the view frustrum which is [-1, -10]
  pos[2] = -3;

  //move it down, look at it from above 
  //pos[1] = -1.3;

  /* start the event handler */
  glutMainLoop();

  return 0;
}




/* this function is called whenever the window needs to be rendered */
void display(void) {

  //clear the screen
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  //clear all modeling transformations 
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity();

  /* The default GL window is x=[-1,1], y= [-1,1], z=[-1,1] with the origin in
     the center.  The view frustrum was set up from z=-1 to z=-10. The
     camera is at (0,0,0) looking along negative z axis.
  */ 
  

 /* First we translate/rotate our local reference system with the
    current translation and rotation that have been accumulated in the
    globals pos[] and theta[]
 */
  glTranslatef(pos[0], pos[1], pos[2]);  
  glRotatef(theta[0], 1,0,0); //rotate theta[0] around x-axis, etc 
  glRotatef(theta[1], 0,1,0);
  glRotatef(theta[2], 0,0,1);
  
  /* Now we draw the scene. Our local reference system is still
     [-1,1]x [-1,1] x [-1,1] */
  
    draw_my_scene(); 
  
  glFlush();
}




void draw_my_scene() {

  //draw a cube of side 1, centered at (0,0,0)
  glColor3fv(green); 
  cube(1, 0); 

  //glColor3fv(red); 
  //  cube(.5); 

  //draw_sphere(.5); 

  //glColor3fv(yellow); 
  //draw_circle(.5); 

  glRotatef(45, 0, 0, 1);

  draw_flower();
  glRotatef(-45, 0, 0, 1);
  draw_sphere(.6);

}


void draw_circle(float rad) {

  int N = 100; 
  float angle = 2* M_PI/N; 
  float x, y; 
  int i;   

  glBegin(GL_POLYGON); 
  for (i=0; i<N; i++) {
    x = rad * cos(angle *i); 
    y = rad * sin (angle * i); 
    glVertex3f(x, y, 0);
  }
  glEnd(); 
}



void draw_flower() {


  float  r1  = .2, r2 = .1; 
  int  PETALS = 10; 

  
  //draw  a stem 
  glColor3fv(green); 
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glBegin(GL_POLYGON); 
  glVertex3f(-.01, -.5, 0);
  glVertex3f(-.01, .5, 0);
  glVertex3f(-.01, .5, 0);
  glVertex3f(.01, -.5, 0);
  glEnd(); 


  //draw middle 
  glTranslatef(0, .5+r1, 0);
  draw_circle(r1);


  glColor3fv(magenta); 
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  float   angle  = 2* M_PI/PETALS; 
  float x, y; 
  int i; 
  for (i=0; i< PETALS; i++) {
    //draw petal i
    x = (r1+r2)*cos (i* angle);
    y = (r1+r2)*sin (i* angle);

    glTranslatef(x,y,0);
    draw_circle(r2);
    glTranslatef(-x, -y, 0); 
  }
  glTranslatef(0, -.5-r1, 0);
  

}


//draw points sampled on a sphere of radius rad
void draw_sphere(float rad) {

  glColor3fv(blue);
  int N = 10; 
  float u = 2*M_PI/N; 
  float v = M_PI/N; 

  int i, j; 
  float x, y, z; 
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      x = rad * cos (i*u) * sin(j*v); 
      y = rad * sin (i*u) * sin(j*v); 
      z = rad * cos(j*v); 
     
      //draw a cube centered at that point (x,y,z) 
      glTranslatef(x,y,z);
      cube(.01, 1); 
      glTranslatef(-x,-y,-z);

    }
  }
}



//draw a  cube centered at origin  [-side,side]^3
void cube(GLfloat s, int mode) {
  
  //the polygon mode 
  if(mode) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  } else {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  } 

  /* back, front faces */
  glBegin(GL_POLYGON); 
  glVertex3f(-s, -s, -s);
  glVertex3f(-s, s, -s);
  glVertex3f(s, s, -s);
  glVertex3f(s, -s, -s); 
  glEnd(); 

  glBegin(GL_POLYGON); 
  glVertex3f(-s, -s, s);
  glVertex3f(-s, s, s);
  glVertex3f(s, s, s);
  glVertex3f(s, -s, s); 
  glEnd(); 

  
  /* left, right faces*/
  glBegin(GL_POLYGON); 
  glVertex3f(-s, -s, -s);
  glVertex3f(-s, -s, s);
  glVertex3f(-s, s, s);
  glVertex3f(-s, s, -s); 
  glEnd(); 

  glBegin(GL_POLYGON); 
  glVertex3f(s, -s, -s);
  glVertex3f(s, -s, s);
  glVertex3f(s, s, s);
  glVertex3f(s, s, -s); 
  glEnd(); 

  
  /* up, down  faces  */
  glBegin(GL_POLYGON); 
  glVertex3f(-s, -s, -s);
  glVertex3f(-s, -s, s);
  glVertex3f(s, -s, s);
  glVertex3f(s, -s, -s); 
  glEnd(); 

  glBegin(GL_POLYGON); 
  glVertex3f(-s, s, -s);
  glVertex3f(-s, s, s);
  glVertex3f(s, s, s);
  glVertex3f(s, s, -s); 
  glEnd(); 
}





/* this function is called whenever  key is pressed */
void keypress(unsigned char key, int x, int y) {

  switch(key) {

 
    
    //ROTATIONS 
  case 'x':
    theta[0] += 5.0; 
    glutPostRedisplay();
    break;
  case 'y':
    theta[1] += 5.0;
    glutPostRedisplay();
    break;
  case 'z':
    theta[2] += 5.0;
    glutPostRedisplay();
    break;
  case 'X':
    theta[0] -= 5.0; 
    glutPostRedisplay();
    break;
  case 'Y':
    theta[1] -= 5.0; 
    glutPostRedisplay();
    break;
  case 'Z':
    theta[2] -= 5.0; 
    glutPostRedisplay();
    break;
    
    //TRANSLATIONS 
    //backward (zoom out)
  case 'b':
    pos[2] -= 0.1; 
    glutPostRedisplay();
    break;
    //forward (zoom in)
  case 'f':
    pos[2] += 0.1; 
    //glTranslatef(0,0, 0.5);
    glutPostRedisplay();
    break;
    //down 
  case 'd': 
     pos[1] -= 0.1; 
    //glTranslatef(0,0.5,0);
    glutPostRedisplay();
    break;
    //up
  case 'u': 
    pos[1] += 0.1; 
    //glTranslatef(0,-0.5,0);
    glutPostRedisplay();
    break;
    //left 
  case 'l':
    pos[0] -= 0.1; 
    glutPostRedisplay();
    break;
    //right
  case 'r':
    pos[0] += 0.1; 
    glutPostRedisplay();
    break;

    case 'q':
    exit(0);
    break;
  } 
}//keypress

