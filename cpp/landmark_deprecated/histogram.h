#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <iostream>
#include <stdlib.h> 
#include <vector>
#include <math.h>
#include <algorithm>


using namespace std;


/* makes an ugly histogram of a vector . 

this can be built up considerably, 
but this is only used to visually check that
some data looks how you'd expect it to.

*/




class box{ 
	public: 
		int c; 
		double r; 
		box(); 
		box(int x, double y);
	};



void histogram(vector <float> r);

	
	
box::box(){}	
box::box(int x, double y){c =x; r=y;}

void histogram(vector <float> r)
{  
 
	double start, end;
	double interval = 1;
	
	double test = *min_element(r.begin(),r.end());
	
	cout << "min element!!: " << test << endl;
	start = floor( *min_element(r.begin(), r.end()));
	end = ceil( *max_element(r.begin(), r.end()));

	vector<box> histArray;

	//very inefficient histogram maker.
	for(double s = start; s <= end; s+=interval){

		
		int count = 0;
		for (vector<float>::iterator it = r.begin(); it != r.end(); it++){
			if (s < (*it) && (*it) < s+interval){
				count++;
			}	
		}
		box n = box(count,s);
		histArray.push_back(n);
	}
	
	for(vector<box>::iterator h = histArray.begin(); h != histArray.end(); h++)
	{
		cout << (*h).c << " of " << (*h).r  << endl;	
	}
	
	
}

#endif