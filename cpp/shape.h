#ifndef SHAPE_H
#define SHAPE_H






#include <iostream>
#include <stdlib.h> 
#include <vector>
#include <math.h>
#include <algorithm>


#include "point.h"



//this will be the shape class
class Shape{
	//make an array of points
	public:
		int initT; //initiated time
	
		vector<Point> points;
		Shape();
		Shape(int);
		Shape(vector<Shape>);
		Shape(float var);
		void printShape();
		friend Shape difference(Shape, Shape);
};


Shape::Shape(){}
Shape::Shape(int r){
	if (r ==0){
		//make a triangle
		//points.push_back(Point(box_muller(0,0.1), box_muller(0,0.1), box_muller(0,0.1)));
		//points.push_back(Point(box_muller(0,0.1), box_muller(1,0.1), box_muller(1,0.1)));
		//points.push_back(Point(box_muller(0,0.1), box_muller(0,0.1), box_muller(1,0.1)));
		points.push_back(Point(0.0, 0.0, 0.0));
		points.push_back(Point(0.0, 1.0, 1.0));
		points.push_back(Point(0.0, 0.0, 1.0));
	}
	
	//for(vector<Point>::iterator it = points.begin(); it != points.end(); it++)
	//	cout<< "another point added" << endl;

	
}

Shape::Shape(float var){

	points.push_back(Point(box_muller(0,var), box_muller(0,var), box_muller(0,var)));
	points.push_back(Point(box_muller(0,var), box_muller(0,var), box_muller(1,var)));
	points.push_back(Point(box_muller(0,var), box_muller(1,var), box_muller(1,var)));
	points.push_back(Point(box_muller(0,var), box_muller(1,var), box_muller(0,var)));

	points.push_back(Point(box_muller(1,var), box_muller(0,var), box_muller(0,var)));
	points.push_back(Point(box_muller(1,var), box_muller(0,var), box_muller(1,var)));
	points.push_back(Point(box_muller(1,var), box_muller(1,var), box_muller(1,var)));
	points.push_back(Point(box_muller(1,var), box_muller(1,var), box_muller(0,var)));
		
}
/* */ 

Shape::Shape(vector<Shape> stack){
	
	//constructor that builds a new average shape out of  a stack of shapes
	
	//This will assume that the individual shape vectors are ordered (rotated)
	// also assume they all have the same number of points.
	
	//this is a vector of points we're adding together for averaging - summing each point
	vector<Point> tt;
	int s = (*stack.begin()).points.size() ;
	for(int i = 0; i < s; i++)
		tt.push_back(Point());
	cout << "SIZE OF : " << (*stack.begin()).points.size() << endl;;
	
	vector<Point>::iterator ttit = tt.begin();
	
	
	for(vector<Shape>::iterator it = stack.begin(); it != stack.end(); it++){
		ttit = tt.begin();
		for(vector<Point>::iterator pit = (*it).points.begin(); pit != (*it).points.end(); pit++){
			(*ttit).add(*pit);
			ttit++;	
		}
	}
	
	//cout<< "trying to add all the points together" <<endl;
	for(vector<Point>::iterator it = tt.begin(); it != tt.end(); it++)
		(*it).divide(float(stack.size()));
	
	points = tt;
	
}
void Shape::printShape(){
	for(vector<Point>::iterator it = points.begin(); it != points.end(); it++)
		(*it).printCoord();
	
}


Shape difference(Shape a, Shape b){
	//for each point in shape a and b, make a new point that is the difference
	
	Shape n;
	vector<Point>::iterator bit = b.points.begin();
	for(vector<Point>::iterator ait = a.points.begin(); ait != a.points.end(); ait++)
	{
		n.points.push_back(differences( (*ait), (*bit) ) );
		bit++;
	}
	cout << "made it here" << endl;
	
	n.printShape();
	
	
}



#endif
