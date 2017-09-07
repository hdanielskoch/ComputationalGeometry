#ifndef __pointStack_h
#define __pointStack_h

#include "geom.h"

typedef struct Stack {
	pointNode *firstPoint;
	int size;
} pointStack;

pointStack *createPointStack();
void push(pointStack *stack, pointNode *node);
pointNode *pop(pointStack *stack);
pointNode *peek(pointStack *stack);
int size(pointStack *stack);
int isEmpty(pointStack *stack);

#endif
