#ifndef BASIC_FUNC
#define BASIC_FUNC

#include "Point3D.h"
#include <stdio.h>
#include <cmath>
using namespace std;

// calculate distance
double distance(double *x, double *y){
	return sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]));
}

double dist2pow(double *x, double *y){
	return (x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]);
}

double dot(double *a, double *b)
{
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

Point3D rotateAndTrans(Point3D pt, double trans[3], double rot[3][3]){
	double x = pt[0];
	double y = pt[1];
	double z = pt[2];
	
	double xt = trans[0] + rot[0][0]*x + rot[0][1]*y + rot[0][2]*z;
	double yt = trans[1] + rot[1][0]*x + rot[1][1]*y + rot[1][2]*z;
	double zt = trans[2] + rot[2][0]*x + rot[2][1]*y + rot[2][2]*z;
	
	Point3D ans(xt, yt, zt);	
	return ans;
}

/****************************************
 * u = rot; t = trans
 ****************************************/
vector<Point3D> rotateAndTrans(vector<Point3D> pts, double trans[3], double rot[3][3]){
	vector<Point3D> ans;
	for (int i = 0; i < pts.size(); i++){
			double x = pts[i][0];
			double y = pts[i][1];
			double z = pts[i][2];
			
			double xt = trans[0] + rot[0][0]*x + rot[0][1]*y + rot[0][2]*z;
			double yt = trans[1] + rot[1][0]*x + rot[1][1]*y + rot[1][2]*z;
			double zt = trans[2] + rot[2][0]*x + rot[2][1]*y + rot[2][2]*z;
			
			Point3D pt(xt, yt, zt);
			ans.push_back(pt);
	}
	
	return ans;
}

/****************************************
 * u = rot;
 ****************************************/
vector<Point3D> rotate(vector<Point3D> pts, Point3D rotation_point, double rot[3][3]){
	vector<Point3D> ans;
	
	double rp_x = rotation_point[0];
	double rp_y = rotation_point[1];
	double rp_z = rotation_point[2];
	
	for (int i = 0; i < pts.size(); i++){
			double x = pts[i][0];
			double y = pts[i][1];
			double z = pts[i][2];
			
			double xt = rp_x + rot[0][0]*(x-rp_x) + rot[0][1]*(y-rp_y) + rot[0][2]*(z-rp_z);
			double yt = rp_y + rot[1][0]*(x-rp_x) + rot[1][1]*(y-rp_y) + rot[1][2]*(z-rp_z);
			double zt = rp_z + rot[2][0]*(x-rp_x) + rot[2][1]*(y-rp_y) + rot[2][2]*(z-rp_z);
			
			Point3D pt(xt, yt, zt);
			ans.push_back(pt);
	}
	
	return ans;
}

void transform(double t[3], double u[3][3], double *x, double *x1)
{
    x1[0]=t[0]+dot(&u[0][0], x);
    x1[1]=t[1]+dot(&u[1][0], x);
    x1[2]=t[2]+dot(&u[2][0], x);
}

void do_rotation(double **x, double **x1, int len, double t[3], double u[3][3])
{
    for(int i=0; i<len; i++)
    {
        transform(t, u, &x[i][0], &x1[i][0]);
    }
}

void transform(double* t, double** u, double *x, double *x1)
{
    x1[0]=t[0]+dot(&u[0][0], x);
    x1[1]=t[1]+dot(&u[1][0], x);
    x1[2]=t[2]+dot(&u[2][0], x);
}

void do_rotation(double **x, double **x1, int len, double *t, double** u)
{
    for(int i=0; i<len; i++)
    {
        transform(t, u, x[i], x1[i]);
    }    
}

// define a 2-D arry
double** new2DArr(int row, int col){
	double **ans=new double*[row];
	for(int i=0;i<row;i++){
		ans[i]=new double[col];
		for(int j=0; j<col; j++)
			ans[i][j] = 0.0;
	}

	return ans;
}

void release2DArr(int n, double ** Arr){
	for(int i = 0; i < n; i++){
		delete [] Arr[i];
	}
	delete [] Arr;
	
	Arr = NULL;
}

char** new2DCharArr(int row, int col){
	char **ans=new char*[row];
	for(int i=0;i<row;i++){
		ans[i]=new char[col];
		for(int j=0; j<col; j++)
			ans[i][j] = '\0';
	}
	
	return ans;
}

void release2DCharArr(int n, char ** Arr){
	for(int i = 0; i < n; i++){
		delete [] Arr[i];
	}
	delete [] Arr;
	
	Arr = NULL;
}

int** new2DIntArr(int row, int col){
	int **ans=new int*[row];
	for(int i=0;i<row;i++){
		ans[i]=new int[col];
		for(int j=0; j<col; j++)
			ans[i][j] = 0;
	}
	
	return ans;
}

int*** new3DIntArr(int row, int col, int thd){
	int ***ans=new int**[row];
	for(int i=0;i<row;i++){
		ans[i]=new int*[col];
		for (int j=0; j<col;j++){
			ans[i][j]=new int[thd];
			for (int k=0; k<thd; k++){
				ans[i][j][k] = 0;
			}
		}
	}
	
	return ans;
}

void release3DIntArr(int row, int col, int*** Arr){
	for(int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			delete[] Arr[i][j];
		}
		delete[] Arr[i];
	}
	delete[] Arr;
	Arr = NULL;
}

void release2DIntArr(int n, int ** Arr){
	for(int i = 0; i < n; i++){
		delete [] Arr[i];
	}
	delete [] Arr;
	Arr = NULL;
}

string** new2DStringArr(int row, int col){
	string **ans=new string*[row];
	for(int i=0;i<row;i++){
		ans[i]=new string[col];
	}

	return ans;
}

void release2DStringArr(int n, string ** Arr){
	for(int i = 0; i < n; i++){
		delete [] Arr[i];
	}
	delete [] Arr;
	Arr = NULL;
}

/*
double abs(double value){
	return value >= 0 ? value : -value;
}
*/

double average(vector<double> vd){
	double ans = vd[0];
	for (int i = 1; i < vd.size(); i++)
		ans += vd[i];
	ans /= vd.size();
	
	return ans;
}

bool isContain(int* arr, int arrlen, int value){
	bool ans = false;
	for (int i = 0; i < arrlen; i++){
		if (value == arr[i]){
			ans = true;
			break;
		}
	}
	
	return ans;
}

bool isContain(vector<int> vec, int value){
	bool ans = false;
	for (int i = 0; i < vec.size(); i++){
		if (value == vec[i]){
			ans = true;
			break;
		}
	}
	
	return ans;
}

bool isSame(vector<int> vecA, vector<int> vecB){
	bool ans = true;
	if (vecA.size() == vecB.size()){
		int i = 0, j = 0;

		bool* isMatchB = new bool[vecB.size()];
		for (i = 0; i < vecA.size(); i++){
			isMatchB[i] = false;
		}
		
		for (i = 0; i < vecA.size(); i++){
			bool isMatch = false;

			for (j = 0; j < vecB.size(); j++){
				if (isMatchB[j])	// just use once time
					continue;

				if (vecA[i] == vecB[j]){
					isMatch = true;
					isMatchB[j] = true;
					break;
				}
			}

			if (!isMatch){
				ans = false;
				break;
			}
		}

		delete[] isMatchB;
		isMatchB = NULL;
	}else{
		return false;
	}

	return ans;
}

int indexOf(vector<int> vec, int value){
	int ans = -1;
	for (int i = 0; i < vec.size(); i++){
		if (value == vec[i]){
			ans = i;
			break;
		}
	}
	
	return ans;
}


/****************************************************************
 * For some applications, it is helpful to be able to make a rotation with a given axis. 
 * Given a unit vector u = (ux, uy, uz), where ux*ux + uy*uy + uz*uz = 1, the matrix for a rotation
 * by an angle of ¦È about an axis in the direction of u is following:
 * 		R = [[cos¦È+ux*ux*(1-cos¦È)		ux*uy*(1-cos¦È)-uz*sin¦È		ux*uz*(1-cos¦È)+uy*sin¦È]
 * 			 [uy*ux*(1-cos¦È)+uz*sin¦È	cos¦È+uy*uy*(1-cos¦È)			uy*uz*(1-cos¦È)-ux*sin¦È]
 * 			 [uz*ux*(1-cos¦È)-uy*sin¦È	uz*uy*(1-cos¦È)+ux*sin¦È		cos¦È+uz*uz*(1-cos¦È)	  ]]
 * 
 * 
 * The detail can be found in https://en.wikipedia.org/wiki/Rotation_matrix
 * @param axis : the O(0, 0, 0) to axis(awx, awy, awz) means the axis or u
 * @param angle : the rotated angle , ¦È
 * @return the rotation matrix
 ***************************************************************/
double** generateRotatedCoordinate(Point3D axis, double ang){		// the memory need to release
	double awx = axis.getX();
	double awy = axis.getY();
	double awz = axis.getZ();
	
	Point3D Opt(0.0, 0.0, 0.0);
	double axisLen = axis.dist(Opt);
	awx /= axisLen;
	awy /= axisLen;
	awz /= axisLen;
			
	double acos=cos(ang);
	double asin=sin(ang);
	double a11=awx*awx*(1-acos)+acos;
	double a12=awx*awy*(1-acos)-awz*asin;
	double a13=awx*awz*(1-acos)+awy*asin;
	double a21=awx*awy*(1-acos)+awz*asin;
	double a22=awy*awy*(1-acos)+acos;
	double a23=awy*awz*(1-acos)-awx*asin;
	double a31=awx*awz*(1-acos)-awy*asin;
	double a32=awy*awz*(1-acos)+awx*asin;
	double a33=awz*awz*(1-acos)+acos;
	
	double** rotmtx = new2DArr(3, 3);
	rotmtx[0][0] = a11;
	rotmtx[0][1] = a12;
	rotmtx[0][2] = a13;
	
	rotmtx[1][0] = a21;
	rotmtx[1][1] = a22;
	rotmtx[1][2] = a23;
	
	rotmtx[2][0] = a31;
	rotmtx[2][1] = a32;
	rotmtx[2][2] = a33;
	
	return rotmtx;
}

double** calculateDisPow2Mtx(vector<Point3D> pts){
	int i = 0, j = 0;
	double** p_disPow2Mtx = new double*[pts.size()];
	for (i = 0; i < pts.size(); i++){
		p_disPow2Mtx[i] = new double[pts.size()];
		for (j = 0; j < pts.size(); j++){
			double dis = 0.0;
			if (i == j) dis = 0.0;
			else if (i > j) dis = p_disPow2Mtx[j][i];
			else dis = pts[i].distPow2(pts[j]);
				p_disPow2Mtx[i][j] = dis;
		}
	}
	return p_disPow2Mtx;
}

/***********************************************************************************************
 * @param pts: the search point map
 * @param p_disPowMtx: the search map data, calculate from pts
 * @return	disPow2_AB: the reference far distance , the distance pow2 between A and B
 *			disPow2_AC: the reference distance, the distance pow2 between A and C
 *			disPow2_BC: the reference distance, the distance pow2 between B and C
 ***********************************************************************************************/
double* calculateThreeSpecicalDistPow2(int lig_size, double** p_disPow2Mtx){
	int i = 0, j = 0;
	int AInd = -1;
	int BInd = -1;
	double disPow2_AB = 0.0;
	for (i = 0; i < lig_size; i++){
		for (j = i+1; j < lig_size; j++){
			if (p_disPow2Mtx[i][j] > disPow2_AB){
				disPow2_AB = p_disPow2Mtx[i][j];
				AInd = i;
				BInd = j;
			}
		}
	}
	
	double disPow2_ACB = 0.0;
	double disPow2_AC = 0.0;
	double disPow2_BC = 0.0;
	for (i = 0; i < lig_size; i++){
		if (i != AInd && i != BInd){
			double disACB = p_disPow2Mtx[i][AInd] + p_disPow2Mtx[i][BInd];
			if (disACB > disPow2_ACB){
				disPow2_ACB = disACB;
				disPow2_AC = p_disPow2Mtx[i][AInd];
				disPow2_BC = p_disPow2Mtx[i][BInd];
			}
		}
	}
	double* ans = new double[3];
	ans[0] = disPow2_AB;
	ans[1] = disPow2_AC;
	ans[2] = disPow2_BC;
	return ans;
}

/***********************************************************************************************
 * @param pts: the search point map
 * @param disPowMtx: the search map data, calculate from pts
 * @param ref_disPow2_AB: the reference far distance , the distance pow2 between A and B
 * @param ref_disPow2_AC: the reference distance, the distance pow2 between A and C
 * @param ref_disPow2_BC: the reference distance, the distance pow2 between B and C
 * @return the similar score, the lower the better
 * @description the time complex is O(n^2)
 ***********************************************************************************************/
double scoreWithThreeSpecialPts(int lig_size, 
	double** disPow2Mtx, double ref_disPow2_far, double ref_disPow2_dis1, double ref_disPow2_dis2){
        int i = 0, j = 0/*, k = 0*/;
	
	double sco = 999999999999999999999.0;
	for (i = 0; i < lig_size; i++){
	        /*int Aind = i; */
		int Bind = -1;
		int Cind = -1;
		double ABDiff = 9999999999999999999.0;
		double ACDiff = 9999999999999999999.0;
		for (j = i+1; j < lig_size; j++){
			double _ABDiff = abs(disPow2Mtx[i][j] - ref_disPow2_far);
			if (_ABDiff < ABDiff){
				ABDiff = _ABDiff;
				Bind = j;
			}
			double _ACDiff = abs(disPow2Mtx[i][j] - ref_disPow2_dis1);
			if (_ACDiff < ACDiff){
				ACDiff = _ACDiff;
				Cind = j;
			}
		}
		
		if (-1 == Bind || -1 == Cind){
			continue;
		}
		double BCDiff = abs(disPow2Mtx[Bind][Cind] - ref_disPow2_dis2);
		double _sco = (ABDiff + ACDiff + BCDiff);
		if (_sco < sco){
			sco = _sco;
		}
	}
	return sco;
}

#endif 
