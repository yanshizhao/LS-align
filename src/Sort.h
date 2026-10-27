#ifndef SORT_H
#define SORT_H

#include <vector>
#include <algorithm>
#include <iterator>
using namespace std;

bool descentFactor(double a, double b){
	if (a > b)
		return true;
	return false;
}

void descentSort(vector<double>& vec){
	sort(vec.begin(), vec.end(), descentFactor);
}

void ascentSort(vector<double>& vec){
	sort(vec.begin(), vec.end());
}

int* descentInsertSortIndex(const vector<double>& vec){
	int i;

	int* index = new int[vec.size()];
	for (i = 0; i < vec.size(); i++)
		index[i] = i;

	for (i = 1; i < vec.size(); i++){
		for (int j = 0; j < i; j++){
			if (vec[index[i]] > vec[index[j]]){
				int tmp = index[i];
				index[i] = index[j];
				index[j] = tmp;
			}
		}
	}

	return index;
}

int* descentInsertSortIndex(const vector<double>& vec, int sn){
	int i;

	if (sn > vec.size()){
		sn = vec.size();
	}

	int* index = new int[vec.size()];
	for (i = 0; i < vec.size(); i++)
		index[i] = i;

	for (i = 1; i < vec.size(); i++){
		for (int j = 0; j < i && j < sn; j++){
			if (vec[index[i]] > vec[index[j]]){
				int tmp = index[i];
				index[i] = index[j];
				index[j] = tmp;
			}
		}
	}

	int* ans = new int[sn];
	for (i = 0; i < sn; i++)
		ans[i] = index[i];

	delete[] index;
	index = NULL;
	return ans;
}

int* ascentInsertSortIndex(const vector<double>& vec){
	int i;

	int* index = new int[vec.size()];
	for (i = 0; i < vec.size(); i++)
		index[i] = i;

	for (i = 1; i < vec.size(); i++){
		for (int j = 0; j < i; j++){
			if (vec[index[i]] < vec[index[j]]){
				int tmp = index[i];
				index[i] = index[j];
				index[j] = tmp;
			}
		}
	}

	return index;
}

int* ascentInsertSortIndex(const vector<double>& vec, int sn){
	int i;

	if (sn > vec.size()){
		sn = vec.size();
	}

	int* index = new int[vec.size()];
	for (i = 0; i < vec.size(); i++)
		index[i] = i;

	for (i = 1; i < vec.size(); i++){
		for (int j = 0; j < i && j < sn; j++){
			if (vec[index[i]] < vec[index[j]]){
				int tmp = index[i];
				index[i] = index[j];
				index[j] = tmp;
			}
		}
	}

	int* ans = new int[sn];
	for (i = 0; i < sn; i++)
		ans[i] = index[i];

	delete[] index;
	index = NULL;
	return ans;
}

int* ascentInsertSortIndex(const vector<int>& vec){
	int i;

	int* index = new int[vec.size()];
	for (i = 0; i < vec.size(); i++)
		index[i] = i;

	for (i = 1; i < vec.size(); i++){
		for (int j = 0; j < i; j++){
			if (vec[index[i]] < vec[index[j]]){
				int tmp = index[i];
				index[i] = index[j];
				index[j] = tmp;
			}
		}
	}

	return index;
}

#endif
