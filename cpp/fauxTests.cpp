//To compile this,  type 'make genData' into the cygwin terminal


#include <iostream>
#include <stdlib.h> 
#include <vector>
#include <math.h>
#include <algorithm>
#include <typeinfo>
#include <string>


#include <gnuplot>


#include "histogram.h"
#include "statFuncts.h"
#include "genData.h"
#include "shape.h"
#include "point.h"




using namespace std;


int main() {
	
	srand(time(NULL));
	//float mean, variance;
	//mean = 2;
	float var = 0.1;
	vector<Shape> patients;
	
	
	///make 20 shapes
	for(int i = 0; i<20; i++){
		patients.push_back(Shape(var));
	}
	
	///print those 20 shapes out
	//for(vector<Shape>::iterator it = patients.begin();it != patients.end(); it++){
	//	(*it).printShape();
	//}
	
	Shape aveShape = Shape(patients);
	aveShape.printShape();
	
	
	//generate a L=100 length random vector
	//supposed to see if it is actually normal
	//vector<float> r = fauxData(mean, variance);
	//histogram(r);
	//normality(r,float(mean));

	
	
	//look at 2 different shapes and see
	//if the differences are normal
	
	//vector<Shape> stack;
	//stack.push_back(Shape(0));
	//stack.push_back(Shape(0));
	
	
	//cout << typeid(stack).name() << endl;
	//cout << typeid(stack[0]).name() << endl;
	
	//stack[0].printShape();
	//stack[1].printShape();
	
	//Shape n = difference(stack[0], stack[1]);
	
	//
	//cout << "make it here 2.0" << endl;
	
	//Point diff = differences(stack[0].points.begin(), stack[1].points.begin());

	
	
	//cout << " the difference vector between 2 points"<< endl;
	//diff.printCoord();
	
	//float n = diff.normEuclid();
	//cout << "norm of this difference vector: " << diff.normEuclid()<< endl;
	
	//Shape avShape = Shape(stack);
	
	//cout << "woohoo this is what the average shape looks like! " <<endl;
	//avShape.printShape();
	//now take 2 triangles and find the average shape. 
	

    return 0;
}




