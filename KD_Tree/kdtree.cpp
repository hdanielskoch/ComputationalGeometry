#include "kdtree.h"
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>

using namespace std;

/*******************************/
//Comparison functions
bool compareX(point2D i, point2D j);
bool compareY(point2D i, point2D j);

/*******************************/
/* KD Tree Constructor */
Kdtree::Kdtree(vector<point2D> p, int n, int cuttype) {

    numNodes = 0;

    vector<point2D> sortedByX = p;
    vector<point2D> sortedByY = p;

    //Sort 2 arrays, one by x and one by y
    sort(sortedByX.begin(), sortedByX.end(), compareX);
    sort(sortedByY.begin(), sortedByY.end(), compareY);

    root = buildKdtree(sortedByX, sortedByY, n, cuttype);
}
/*******************************/
/* KD Tree Destructor */
Kdtree::~Kdtree() {
}

/*******************************/
/* TreeNode Constructor */
TreeNode::TreeNode(point2D p, char type) {
    this->p = p;
    this->type = type;
}
/*******************************/
/* TreeNode Destructors */
TreeNode::~TreeNode() {
}
/*******************************/
/* Print information about node including type and location*/
int TreeNode::printNode(){
    if(type == 'h'){
        cout << "Horizontal line at y coord: " << p.y << endl;
    } else if (type == 'v'){
        cout << "Vertical line at x coord: " << p.x << endl;
    } else{
        cout << "Point at: (" << p.x << ", " << p.y << ")" << endl;
    }
    return 0;
}
/*******************************/

/* Compare functions for sorting vectors of point2D objects by X & Y */
/* Takes two points as input to compare */
bool compareX(point2D i, point2D j) {
    //x values are different
    if (i.x != j.x)
        return (i.x < j.x);
    //x values are the same, compare by y
    else
        return (i.y < j.y);
}
bool compareY(point2D i, point2D j) {
    //y values are different
    if (i.y != j.y)
        return (i.y < j.y);
    //y values are the same, compare by x
    else
        return (i.x < j.x);
}

/*******************************/
/* Kdtree Class Functions */

/* **************************************** */
//Print Each node of tree
void Kdtree::printTree(TreeNode *node){
    //Base Case
    if((node->left == NULL) && (node->right == NULL)){
        node->printNode();
        return;
    //Left node is empty
    }else if(node->left == NULL){
        node->printNode();
        printTree(node->right);
    //Right node is empty
    }else if(node->right == NULL){
        node->printNode();
        printTree(node->left);
    //Both nodes are full
    }else{
        node->printNode();
        printTree(node->left);
        printTree(node->right);
    }
}
/* **************************************** */
// Compute height of KD Tree
int Kdtree::getHeight(TreeNode* node){
    if(node == NULL){
        return 0;
    } else{
        int left = getHeight(node->left);
        int right = getHeight(node->right);
        int maxHeight = max(left,right) + 1;
        return maxHeight;
    }
}

/* Builds KD Tree recursively until leaves with one point are reached */
/* Takes inputs of arrays sorted by x, by y, the size, 
   and whether the current node has a verytical or horizontal cut */
TreeNode* Kdtree::buildKdtree(vector<point2D> sortedByX, vector<point2D> sortedByY, int n, int cuttype) {
    //Base Case: Node does not contain any points
    if (n <= 0) {
        return NULL;
    }
    //Base Case: Only one point. Create leaf node
    if (n == 1) {
        TreeNode* x = new TreeNode(sortedByX[0], 'l');
        x->left = NULL;
        x->right = NULL;
        numNodes++;
        return x;
    }
    //Horizontal cut
    if (cuttype == 0) {
        //Create child node
        TreeNode* child = new TreeNode(sortedByY[sortedByY.size()/2], 'h');

        //cut sortedByY vector in half
        vector<point2D>::iterator mid = sortedByY.begin();
        advance(mid, sortedByY.size()/2);
        vector<point2D> P1(sortedByY.begin(), mid); //first half
        vector<point2D> P2(mid, sortedByY.end()); //second half

        //sort P1 by X & Y
        vector<point2D> P1X = P1;
        vector<point2D> P1Y = P1;
        sort(P1X.begin(), P1X.end(), compareX);
        sort(P1Y.begin(), P1Y.end(), compareY);

        //sort P2 by X & Y
        vector<point2D> P2X = P2;
        vector<point2D> P2Y = P2;
        sort(P2X.begin(), P2X.end(), compareX);
        sort(P2Y.begin(), P2Y.end(), compareY);

        //make recursive call to left
        child->left = buildKdtree(P1X, P1Y, P1.size(), 1);
        //make recursive call to right
        child->right = buildKdtree(P2X, P2Y, P2.size(), 1);

        //Add two nodes to number of nodes
        numNodes ++;

        return child;

    //Vertical cut
    } else {
        //Create child node
        TreeNode* child = new TreeNode(sortedByX[sortedByX.size()/2], 'v');

        //cut sortedByY vector in half
        vector<point2D>::iterator mid = sortedByX.begin();
        advance(mid, sortedByX.size()/2);
        vector<point2D> P1(sortedByX.begin(), mid); //first half
        vector<point2D> P2(mid, sortedByX.end()); //second half

        //sort P1 by X & Y
        vector<point2D> P1X = P1;
        vector<point2D> P1Y = P1;
        sort(P1X.begin(), P1X.end(), compareX);
        sort(P1Y.begin(), P1Y.end(), compareY);

        //sort P2 by X & Y
        vector<point2D> P2X = P2;
        vector<point2D> P2Y = P2;
        sort(P2X.begin(), P2X.end(), compareX);
        sort(P2Y.begin(), P2Y.end(), compareY);

        //make recursive call to left
        child->left = buildKdtree(P1X, P1Y, P1.size(), 0);
        //make recursive call to right
        child->right = buildKdtree(P2X, P2Y, P2.size(), 0);

        //Add two nodes to number of nodes
        numNodes ++;

        return child;
    }
}
/*******************************/

/*******************************/
/* TreeNode Class Functions */

/*******************************/
