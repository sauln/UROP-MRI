//this will be the point class

#ifndef POINT_H
#define POINT_H


#include <iostream>
#include <stdlib.h> 
#include <vector>
#include <math.h>
#include <algorithm>

class Point{
		
	float x;
	float y;
	float z;
	
	public:
		
		Point();
		Point(float, float, float);
		void setValues(float, float, float);
		void printCoord();
		float normEuclid();   //returns the euclidean norm of our point
		void add(Point next); //adds another to our point
		void divide(float f); // divides our point by f pointwise
		
		friend Point differences(Point f, Point g); //returns a point that is the difference between the two points.
};

void Point::printCoord(){
	cout << "x: "<< x << "   y:  " << y<< "   z:  " << z<< endl;
}
void Point::setValues(float xi, float yi, float zi){
	x = xi;
	y = yi;
	z = zi;
}



Point::Point(float xi, float yi, float zi){
	setValues(xi,yi,zi);
}
float Point::normEuclid(){
	return sqrt(pow(x,2) + pow(y,2) + pow(z,2));
}




void Point::divide(float f)
{
	x = x/f;
	y = y/f;
	z = z/f;
}
void Point::add(Point next)
{
	x += next.x;
	y += next.y;
	z += next.z;
	
}


Point::Point(){
	setValues(0.0,0.0,0.0);
}

Point differences(Point f, Point g){
	
	Point diff;
	
	diff.x = f.x - g.x;
	diff.y = f.y - g.y;
	diff.z = f.z - g.z;
	
	return diff;
}




#endif 
