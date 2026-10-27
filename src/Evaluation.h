#ifndef EVALUATION_H
#define EVALUATION_H

#include <vector>
#include "Point3D.h"
#include <cmath>
using namespace std;

/*********************************************************
 * Notice: the size of aPts must equal to that of bPts
 *********************************************************/
double RMSD(vector<Point3D> aPts, vector<Point3D> bPts){
	double ans = 0.0;
	for (int i = 0; i < aPts.size(); i++){
		ans += aPts[i].distPow2(bPts[i]);
	}
	ans = sqrt(ans/aPts.size());
	
	return ans;
}

/*********************************************************
 * the RMSD upper bound value, which 
 * is the normal RMSD.
 * Notice: the size of aPts must equal to that of bPts
 *********************************************************/
double RMSD_ub(vector<Point3D> aPts, vector<Point3D> bPts){
	double ans = 0.0;
	for (int i = 0; i < aPts.size(); i++){
		ans += aPts[i].distPow2(bPts[i]);
	}
	ans = sqrt(ans/aPts.size());
	
	return ans;
}

/*********************************************************
 * the RMSD lower bound value, whose definition can be found
 * in the following paper: 
 * 	"Software News and Update AutoDock Vina: Improving the Speed 
 * 	and Accuracy of Docking with a New Scoring Function, Efficient 
 * 	Optimization, and Multithreading".
 *********************************************************/
double RMSD_lb(vector<Point3D> aPts, vector<Point3D> bPts){
	int i = 0, j = 0;
	
	double** p_disPow2Mtx = new double*[aPts.size()];
	for (i = 0; i < aPts.size(); i++){
		p_disPow2Mtx[i] = new double[bPts.size()];
		for (j = 0; j < bPts.size(); j++){
			p_disPow2Mtx[i][j] = aPts[i].distPow2(bPts[j]);
		}
	}
	
	double ansA_B = 0.0;
	for (i = 0; i < aPts.size(); i++){
		double minDisPow2 = 99999999999999999.0;
		for (j = 0; j < bPts.size(); j++){
			if (minDisPow2 > p_disPow2Mtx[i][j]){
				minDisPow2 = p_disPow2Mtx[i][j];
			}
		}
		ansA_B += minDisPow2;
	}
	ansA_B = sqrt(ansA_B / aPts.size());

	double ansB_A = 0.0;
	for (j = 0; j < bPts.size(); j++){
		double minDisPow2 = 99999999999999999.0;
		for (i = 0; i < aPts.size(); i++){
			if (minDisPow2 > p_disPow2Mtx[i][j]){
				minDisPow2 = p_disPow2Mtx[i][j];
			}
		}
		ansB_A += minDisPow2;
	}
	ansB_A = sqrt(ansB_A / bPts.size());
	
	// release temp memory
	release2DArr(aPts.size(), p_disPow2Mtx);
	return ansA_B > ansB_A ? ansA_B : ansB_A;
}

/*********************************************************
 * the RMSD lower bound value, whose definition can be found
 * in the following paper: 
 * 	"Software News and Update AutoDock Vina: Improving the Speed 
 * 	and Accuracy of Docking with a New Scoring Function, Efficient 
 * 	Optimization, and Multithreading".
 * Here, we chose the minimum value.
 *********************************************************/
double RMSD_lb_min(vector<Point3D> aPts, vector<Point3D> bPts){
	int i = 0, j = 0;
	
	double** p_disPow2Mtx = new double*[aPts.size()];
	for (i = 0; i < aPts.size(); i++){
		p_disPow2Mtx[i] = new double[bPts.size()];
		for (j = 0; j < bPts.size(); j++){
			p_disPow2Mtx[i][j] = aPts[i].distPow2(bPts[j]);
		}
	}
	
	double ansA_B = 0.0;
	for (i = 0; i < aPts.size(); i++){
		double minDisPow2 = 99999999999999999.0;
		for (j = 0; j < bPts.size(); j++){
			if (minDisPow2 > p_disPow2Mtx[i][j]){
				minDisPow2 = p_disPow2Mtx[i][j];
			}
		}
		ansA_B += minDisPow2;
	}
	ansA_B = sqrt(ansA_B / aPts.size());

	double ansB_A = 0.0;
	for (j = 0; j < bPts.size(); j++){
		double minDisPow2 = 99999999999999999.0;
		for (i = 0; i < aPts.size(); i++){
			if (minDisPow2 > p_disPow2Mtx[i][j]){
				minDisPow2 = p_disPow2Mtx[i][j];
			}
		}
		ansB_A += minDisPow2;
	}
	ansB_A = sqrt(ansB_A / bPts.size());
	
	// release temp memory
	release2DArr(aPts.size(), p_disPow2Mtx);
	return ansA_B > ansB_A ? ansB_A : ansA_B;
}

#endif
