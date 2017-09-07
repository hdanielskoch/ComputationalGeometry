#include "pointStack.h" 
#include "geom.h"

#include <stdlib.h>
#include <stdio.h>

/*

A stack class with methods. All code here was adapted from
my data structures implementation of a stack

*/

pointStack *createPointStack() {
	pointStack *stack = (pointStack *) malloc(sizeof(pointStack));

	if (stack == NULL) {
		printf("malloc error in creating stack\n");

		return NULL;
	}

	stack->firstPoint = NULL;
	
	stack->size = 0;

	return stack;
}

// "Push" pushes new pointNode nodes onto the stack.
// It has no return type, and takes the stack and
// new node as input parameters.
void push(pointStack *stack, pointNode *node) {
	node->next = stack->firstPoint;

	stack->firstPoint = node;
	
	++stack->size; //keep track of your stack's size
}

// "Pop" pops off the top node from the stack
// and resets the stack to point to the new top.
// the function returns the pointNode move node that
// is being popped, so that it can be used in other
// sections of this program if needed. The function takes the
// stack as its input parameter.
pointNode *pop(pointStack *stack) {
	//if there is nothing to pop, return NULL
	if (stack->size == 0)
		return NULL;

	//keep track of what you're popping off
	pointNode *node = stack->firstPoint;

	//reset the stack to point to the new head
	stack->firstPoint = stack->firstPoint->next;

	--stack->size;

	//if the stack is now empty, make it point to nothing
	if (stack->size == 0)
		stack->firstPoint = NULL;

	//make sure that what you're popping off is no longer
	//linked to the stack's linked list
	node->next = NULL;

	return node;
}

// "Peek" lets you see what is on top without
// altering it in any way. It returns a pointNode,
// and takes the stack as its input parameter.
pointNode *peek(pointStack *stack) {
	if (stack->size == 0)
		return NULL;

	return stack->firstPoint;
}

// "Size" lets the program see how big the stack is.
// It returns the size as an int number, and takes
// the stack as its only input parameter.
int size(pointStack *stack) {
	return stack->size;
}

// "isEmpty" is very self explanatory. It tells you
// if the stack is empty. Boolean return type, stack
// input.
int isEmpty(pointStack *stack) {
	if (stack->size == 0)
		return 1;
	return 0;
}
