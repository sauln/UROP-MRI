#ifndef STATFUNCTS_H
#define STATFUNCTS_H


#include <iostream>
#include <stdlib.h> 
#include <vector>
#include <math.h>
#include <algorithm>

using namespace std;


void jarqueBeraTest(vector<float> r);
void normality(vector<float> r, float actual);
float varianceIsh(vector<float> r, int pow);
float sampleMean(vector<float> r);
float sampleVariance(vector<float> r);
float confidenceOfMean(vector<float> r, float actual, float sd);





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



#endif