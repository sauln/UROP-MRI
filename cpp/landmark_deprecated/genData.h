
#include <iostream>
#include <stdlib.h> 
#include <vector>
#include <math.h>
#include <algorithm>

#include "boxmuller.h"


using namespace std;



/*

pseduocode of the average shape calculation

build a vector of points that is long enough.

for each shape
	add each point to the vector of points



*/



vector<float> dim1(float mean, float variance);



vector<float> dim1(float mean, float variance){
	
	int length = 100;
	vector<float> r;
	for( int i =0; i < length; i++)
	{
		r.push_back(box_muller(mean, sqrt(variance)));
	}
	
	return r;
	
}




