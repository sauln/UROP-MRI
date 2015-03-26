//To compile this,  type 'make genData' into the cygwin terminal





#include <iostream>
#include <stdlib.h> 
#include <vector>
#include <math.h>
#include <algorithm>

#include "histogram.h"

using namespace std;


float sampleMean(vector<float> r);
float varianceIsh(vector<float> r, int pow);
void jarqueBeraTest(vector<float> r);

float box_muller(float m, float s);



float box_muller(float m, float s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * float(rand())/RAND_MAX - 1.0;
			x2 = 2.0 * float(rand())/RAND_MAX - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}

class Point{
	
		float x;
		float y;
		float z;
		
	public:
		Point();
		Point(float, float, float);
		void setValues(float, float, float);
		void printCoord();
};
Point::Point(){}
Point::Point(float xi, float yi, float zi){
	setValues(xi,yi,zi);
}
void Point::printCoord(){
	cout << x << y<< z<< endl;
}
void Point::setValues(float xi, float yi, float zi){
	x = xi;
	y = yi;
	z = zi;
}

class Shape{
	//make an array of points
	public:
		int initT; //initiated time
	
		vector<Point> first;
	
};

void jarqueBeraTest(vector<float> r){

	float mThree = varianceIsh(r,3);
	float mFour = varianceIsh(r,4);
	
	float eThree = pow(varianceIsh(r,2), 1.5);
	float eFour = pow( varianceIsh(r,2),2);
	
	float S = mThree / eThree;
	float K = mFour  / eFour;
	float jb = (r.size()/6)*(pow(S,2) + pow( .25*(K-3), 2));
	
	cout << "Jarque-Bera Test gives: " << jb << endl;
	
}

float varianceIsh(vector<float> r, int p){
	float tot = 0;
	float s = sampleMean(r);
	
	for(vector<float>::iterator it = r.begin(); it != r.end(); it++){
		tot += pow((*it) - s,p);
	}
	
	return tot/r.size();
	
}

float sampleMean(vector<float> r){
	
	float sum =0;
	for (vector<float>::iterator it = r.begin(); it != r.end(); it++)
		sum += *it;
	
	return sum/r.size();
	
}

float sampleVariance(vector<float> r){
	float sum = 0.0;
	float sumsqr = 0.0;
	for (vector<float>::iterator it = r.begin(); it != r.end(); it++){
		sumsqr += (*it)*(*it);
		sum += *it;
	}
	
	return (sumsqr - (sum*sum)/r.size())/(r.size()-1);
}

/*float standardError(vector<float> r){
	for (vector<float>::iterator it = r.begin(); it != r.end(); it++){
		sum += *it;
		sumsqr += (*it)*(*it);
	}
	float standardError = float(sd / sqrt(r.size()));
}*/

float confidenceOfMean(vector<float> r, float actual, float sd){
		

	float meanHat = sampleMean(r);
	float varHat = sampleVariance(r);
	
	
	float standardError = float(sd / sqrt(r.size()));
	
	float sampleError = float(varHat / sqrt(r.size()));
	cout << "standard error: " << standardError <<endl ;
	cout << "sample error: " << sampleError << endl;
	
	float lowerBound = meanHat - 2*standardError;
	float upperBound = meanHat + 2*standardError;
	
	cout << "95% sure that mean is between " << lowerBound << " and " << upperBound << endl; 
	
	
	/*if(lowerBound < actual && upperBound > actual)
		{cout << "mean: " << actual << " is in bounds" << endl;}
	else
		{cout <<"mean: " << actual << " is out of bounds"  << endl;}
	*/
	return meanHat;
}

void normality(vector<float> r, float actual){
	//take a bunch of random numbers and see if they fall on a normal curve.
	
	//first find the mean,
	
	/*  Asuume known variance!  variance = 1 */

	float variance = 1; // assume known
	float sd = variance*variance;
	
	float mean = confidenceOfMean(r, actual, sd);
	
	jarqueBeraTest(r);
	/*
	variance = (sumsqr - (sum*sum)/r.size())/(r.size()-1);
	
	int count =0;
	for (vector<float>::iterator it = r.begin(); it != r.end(); it++){
		if( abs(*it) >= (variance*variance)   ){
			count++;
		}	
	}
	
	
	cout << "proportion of items out of 2sc: " << float(count)/total << endl;
	//then find the variance..
	
	//then fit to a normal curve with these mean and variance
	
	cout << "success" << endl;
	*/
	
}

vector<float> fauxData(float mean, float variance){
	
	int length = 100;
	vector<float> r;
	for( int i =0; i < length; i++)
	{
		r.push_back(box_muller(mean, sqrt(variance)));
	}
	
	return r;
	
}

int main() {
	srand(time(NULL));
	Shape sh;
	sh.initT = 0;
	sh.first.push_back(Point(1,1,0));
	sh.first.push_back(Point(1,2,0));
	sh.first.push_back(Point(2,1,0));
	sh.first.push_back(Point(2,2,0));
		
	float mean, variance;
	mean = 2;
	variance = 1;
	
	vector<float> r = fauxData(mean, variance);
	
	histogram(r);

	normality(r,float(mean));

	
    return 0;
}




