#ifndef SHORTESTROAD_H
#define SHORTESTROAD_H

#include "basic_fun.h"
#include <climits>

using namespace std;

class ShortestRoad{
private:
	double** adjmap;
	int* p_node_num;
	int*** allRoad8Floyd;
	int** allRoadLen8Floyd;

	vector<int> longestShortestRoad8Floyd;
	bool m_isValid;
public:
	// Upper bound on the atom count of one ligand: 300^3 ints ~= 108 MB.
	static const int MAX_FLOYD_NODES = 300;

	ShortestRoad(double** adjmap, int node_num){
		// Reject before allocating: n^3 ints is ~108 MB at n = 300.
		// The n > INT_MAX/n/n test keeps the n^3 product from overflowing.
		if (node_num <= 0 || node_num > MAX_FLOYD_NODES
			|| node_num > INT_MAX / node_num / node_num){
			this->adjmap = NULL;
			p_node_num = NULL;
			allRoad8Floyd = NULL;
			allRoadLen8Floyd = NULL;
			m_isValid = false;
			return;
		}
		m_isValid = true;

		this->adjmap = new2DArr(node_num, node_num);
		for (int i = 0; i < node_num; i++){
			for (int j = 0; j < node_num; j++){
				this->adjmap[i][j] = adjmap[i][j];
			}
		}
		
		this->p_node_num = new int[1];
		p_node_num[0] = node_num;
		
		allRoad8Floyd = new3DIntArr(node_num, node_num, node_num);
		allRoadLen8Floyd = new2DIntArr(node_num, node_num);
		
		floyd4LongestLengthRoadWithShortestDis();
	}		
	
	int getNodeNum(){
		return (NULL != p_node_num) ? *p_node_num : 0;
	}
	
	bool isValid() const {
		return m_isValid;
	}
	
	vector<int> getLongestShortestRoad8Floyd(){
		return longestShortestRoad8Floyd;
	}

	// Shortest hop count from i to j; -1 when unreachable.
	// The stored path includes the start atom, so a direct edge is len 2.
	int getHopCount(int i, int j) const {
		if (!m_isValid || i < 0 || j < 0
			|| i >= *p_node_num || j >= *p_node_num){
			return -1;
		}
		if (i == j) return 0;
		int len = 0;
		while (len < *p_node_num && allRoad8Floyd[i][j][len] >= 0){
			len++;
		}
		return (len <= 1) ? -1 : len - 1;
	}

	~ShortestRoad(){
		if (NULL != adjmap){
			release2DArr(*p_node_num, adjmap);
		}
		
		if (NULL != allRoad8Floyd){
			 release3DIntArr(*p_node_num, *p_node_num, allRoad8Floyd);
		}
		if (NULL != allRoadLen8Floyd){
			release2DIntArr(*p_node_num, allRoadLen8Floyd);
		}
		
		delete[] p_node_num;
	}

	int*** getAllRoad8Floyd(){
		return this->allRoad8Floyd;
	}

private: 	
	void floyd4LongestLengthRoadWithShortestDis(){
		int i = 0;

		double** disMtx = new2DArr(*p_node_num, *p_node_num);
		int** spot = new2DIntArr(*p_node_num, *p_node_num);
		
		for (i = 0; i < *p_node_num; i++) {
			for (int j = 0; j < *p_node_num; j++) {
				spot[i][j] = -1;
				disMtx[i][j] = adjmap[i][j];
				
				if (disMtx[i][j] <= 0)
					disMtx[i][j] = 999999999999999;
				if (i == j)
					disMtx[i][j] = 0;
			}
		}

		for (int k = 0; k < *p_node_num; k++){
			for (int i = 0; i < *p_node_num; i++){
				for (int j = 0; j < *p_node_num; j++){
					if (disMtx[i][j] > disMtx[i][k] + disMtx[k][j]) {
						disMtx[i][j] = disMtx[i][k] + disMtx[k][j];
						spot[i][j] = k;
					}
				}
			}
		}
	
		// search the longest length road with the longest shortest distance
		int maxRoadLen = 0;
		double corrRoadDis = 0;
		int* maxLenRoadInfo = NULL; // road info is "maxLenRoadInfo[0] -> ... -> maxLenRoadInfo[i] -> maxLenRoadInfo[i+1] -> ... -> maxLenRoadInfo[N]"
		for (i = 0; i < *p_node_num; i++){
			for (int j = 0; j < *p_node_num; j++){
				int* tmpPath = new int[*p_node_num];		// do not delete it
				for (int k = 0; k < *p_node_num; k++){
					tmpPath[k] = -1;
				}
				int* roadLen = new int[1];
				roadLen[0] = 0;
				
				tmpPath[roadLen[0]++] = i;
				recursionSearchPath8FloydAlgorithm(spot, i, j, tmpPath, roadLen);
				
				double dis = disMtx[i][j];
				if (roadLen[0] > maxRoadLen){
					maxRoadLen = roadLen[0];
					maxLenRoadInfo = tmpPath;
					
					corrRoadDis = dis;
				}else if (roadLen[0] == maxRoadLen){
					if (dis > corrRoadDis){
						maxRoadLen = roadLen[0];
						maxLenRoadInfo = tmpPath;
						
						corrRoadDis = dis;
					}
				}
				
				delete[] allRoad8Floyd[i][j];
				allRoad8Floyd[i][j] = tmpPath;
				allRoadLen8Floyd[i][j] = roadLen[0];
				
				delete[] roadLen;
				roadLen = NULL;
			}
		}
	
		for (i = 0; i < *p_node_num; i++){
			if (-1 == maxLenRoadInfo[i])
				break;
			longestShortestRoad8Floyd.push_back(maxLenRoadInfo[i]);
		}

		release2DIntArr(*p_node_num, spot);
		release2DArr(*p_node_num, disMtx);
	}
	
	/****************************************************************************
	 * @param spot : the i to j throw node index, floyd algorithm spot matrix
	 * @param i : the 1th position
	 * @param j : the 2nd position
	 * @param onePath : record the shortest road from i, j
	 * @param roadLen : roadLen[0] record the road length
	 ****************************************************************************/
	void recursionSearchPath8FloydAlgorithm(int** spot, int i, int j, int* path, int* roadLen) {
		if (i == j)return;
		if (spot[i][j] == -1){
			// spot == -1 covers both "direct edge" and "no path".  Only a
			// real edge (adjmap entry > 0) appends j; an unreachable pair
			// keeps the path at just the start atom instead of faking [i, j].
			if (adjmap[i][j] > 0.0)
				path[roadLen[0]++] = j;
		}
		else {
			recursionSearchPath8FloydAlgorithm(spot, i, spot[i][j], path, roadLen);
			recursionSearchPath8FloydAlgorithm(spot, spot[i][j], j, path, roadLen);
		}
	}
	
	ShortestRoad(const ShortestRoad&);
	ShortestRoad& operator=(const ShortestRoad&);
};

#endif
