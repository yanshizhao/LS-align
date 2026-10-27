#ifndef POINT3D_H
#define POINT3D_H

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

class Point3D{
private:
	double x;
	double y;
	double z;
public:
	Point3D(double x, double y, double z){
		this->x = x;
		this->y = y;
		this->z = z;
	}

	Point3D(const Point3D& point){
		this->x = point.x;
		this->y = point.y;
		this->z = point.z;
	}	

	Point3D(): x(0.0), y(0.0), z(0.0){}
		
	~Point3D(){}
		
	Point3D* operator+ (const Point3D& point){
		double ans_x = x + point.x;
		double ans_y = y + point.y;
		double ans_z = z + point.z;
		
		return new Point3D(ans_x, ans_y, ans_z);
	}
	
	Point3D* operator- (const Point3D& point){
		double ans_x = x - point.x;
		double ans_y = y - point.y;
		double ans_z = z - point.z;
		
		return new Point3D(ans_x, ans_y, ans_z);
	}

	Point3D* operator/ (int n){
		double ans_x = x / n;
		double ans_y = y / n;
		double ans_z = z / n;
		
		return new Point3D(ans_x, ans_y, ans_z);
	}
	
	/***************************************
	 * @param ind: The ind value must in {0, 1, 2}
	 * @return ans[0] = x; ans[1] = y; ans[2] = z;
	 ***************************************/
	double 	operator[](int ind){
		double ans;
		if (0 == ind)
			ans = x;
		else if (1 == ind)
			ans = y;
		else if (2 == ind)
			ans = z;
		return ans;
	}
	
	double getX(){
		return x;
	}
	
	double getY(){
		return y;
	}
	
	double getZ(){
		return z;
	}
	
	double get(int ind){
		double ans;
		if (0 == ind)
			ans = x;
		else if (1 == ind)
			ans = y;
		else if (2 == ind)
			ans = z;
		return ans;
	}
	
	/*************************************************
	 * calculate the distance between pt and this Point3D
	 *************************************************/
	double dist(const Point3D& point){
		double distPow2 = (x - point.x)*(x - point.x) + (y - point.y)*(y - point.y) + (z - point.z)*(z - point.z);
		return sqrt(distPow2);
	}

	/*************************************************
	 * calculate the distance pow 2 between pt and this Point3D
	 *************************************************/
	double distPow2(const Point3D& point){
		double distPow2 = (x - point.x)*(x - point.x) + (y - point.y)*(y - point.y) + (z - point.z)*(z - point.z);
		return distPow2;
	}
	
	/************************************************************
	 * calculate the nearest distance pow 2 between pt and each 
	 * Point3D in pts
	 ************************************************************/
	double distPow2(const vector<Point3D>& pts){
		double minDistPow2 = 999999999999.0;
		for (int j = 0; j < pts.size(); j++){
			double tmp = this->distPow2(pts[j]); 
			if (tmp < minDistPow2){
				minDistPow2 = tmp;
			}
		}
		return minDistPow2;
	}

	string toString(){
		char buf[100];
		sprintf(buf, "%.6f\t%.6f\t%.6f", x, y, z);
		
		string ans(buf);
		return ans;
	}
};

#endif // POINT3D_H
