Convex Hull in 3D

By Bo Bleckel and Henry Daniels Koch
Computational Geometry
Bowdoin College
February, 2017

This code implements a brute-force algorithm for finding the convex hull of a set of points in 3D. All code written in geom.cpp/.h is attributed to Bo Bleckel and Henry Daniels Koch. Code written in hull3d.cpp is attributed to Laura Toma, Bo Bleckel, Henry Daniels Koch, and members of the CS3250 class.

The brute force algorithm checks all triplets of points to determine if the triplet makes up an extreme face of the 3d hull. To determine an extreme face, the algorithm checks if all points are to the left of the plane using the 4 points in question’s signed volume. 

We note that the convex hull we find does not triangulate the faces. All possible point triplets that are coplanar are returned in the hull. A possible solution is to allow faces to made up of more than 3 points, run a 2D convex hull algorithm on the set of coplanar points and then return the points in the 2D hull. To return triangulated faces of coplanar points, one could keep track of faces and only add a point to 2 existing points on a face if the two faces don’t overlap. 

To execute this program, compile using the given Makefile (simply type make from the command line) and then execute the hull3d executable that will be made. The program requires that the user enter the number of points that they would like to make as a second command line argument, after the executable. Additionally, to test out different configurations of points, enter values 0-4 and ‘i’ to try premade test cases. We note that for the test case of points in a plane, we made all the coplanar faces yellow to produce one face. In the rest of the test cases, we use random colors and make the colors transparent. 