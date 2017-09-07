#ifndef KDTREE_H
#define KDTREE_H

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <ctime>
#include <climits>
#include <cmath>
#include <fstream>

using namespace std;

typedef struct _point2d {
	double x,y;
} point2D;

class TreeNode {
    public:
        TreeNode(point2D p, char type);
        ~TreeNode();
        TreeNode *left, *right;
        char type;
        point2D p;
        int printNode();
};

class Kdtree {
    public:
        Kdtree(vector<point2D> p, int n, int cuttype);
        ~Kdtree();
        TreeNode* root;
        void printTree(TreeNode *node);
        int getHeight(TreeNode* node);
        int numNodes;
        TreeNode* buildKdtree(vector<point2D> sortedbyx, vector<point2D> sortedbyy, int n, int cuttype);
};

#endif
