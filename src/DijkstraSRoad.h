#ifndef DIJKSTRA_S_ROAD_H
#define DIJKSTRA_S_ROAD_H

#include "basic_fun.h"
#include <vector>
using namespace std;

class DijkstraSRoad{
private:
	double** _adjmap;
	int node_num;

public:
	DijkstraSRoad(double** adjmap, int node_num){	// not mem copy
		this->_adjmap = adjmap;
		this->node_num = node_num;
	}

	/****************************************************************************
	 * @param start_pos : the road start index (index start from 0)
	 * @param end_pos : the road end index (index start from 0)
	 * @param ignoreDirectRoad : whether ignore the directly road, i.e., {start_pos}-{end_pos}
	 ****************************************************************************/
	vector<int> getRoad(int start_pos, int end_pos, bool ignoreDirectRoad){
		int i = 0, j = 0;
		
		//-----------------------------------------------------------------------//
		//--------  Using dijkstra algorthm to search the shortest road----------//
		//-----------------------------------------------------------------------//
		double* dist = new double[node_num];
		int* path = new int[node_num];
		for (i = 0; i < node_num; i++){
			dist[i] = -1.0;
			path[i] = -1;
		}

		double** adjmap = new double*[node_num];
		for (i = 0; i < node_num; i++){
			adjmap[i] = new double[node_num];

			for (j = 0; j < node_num; j++){
				if (_adjmap[i][j] <= 0){
					adjmap[i][j] = -1;
				}
				
				if (i == j)
					adjmap[i][j] = 0;
				else adjmap[i][j] = _adjmap[i][j];

				if (ignoreDirectRoad){
					if ((i == start_pos && j == end_pos)
						|| (i == end_pos && j == start_pos)){
						adjmap[i][j] = -1;
					}
				}
			}
		}
		
		dijkstra(adjmap, node_num, dist, start_pos, path);

		//-----------------------------------------------------------------------//
		//--------  Extract the shortest road                          ----------//
		//-----------------------------------------------------------------------//
		vector<int> ans;

		ans.push_back(end_pos);
		int previous_pos = path[end_pos];
		bool isWalk2Start = false;
		while (true){
			if (previous_pos == start_pos){
				ans.push_back(start_pos);
				isWalk2Start = true;
				break;
			}
			
			if (previous_pos < 0){
				break;
			}
			
			ans.push_back(previous_pos);
			
			previous_pos = path[previous_pos];
		}
		
		if (!isWalk2Start){
			ans.erase(ans.begin(), ans.end());
		}

		delete[] dist;
		dist = NULL;
		delete[] path;
		path = NULL;
		release2DArr(node_num, adjmap);

		return ans;
	}

private:
	void dijkstra(double** adjmap, int node_num, double* dist, int start_position, int* path) {
		int i = 0; 

		int nNode = node_num; // the node number

		bool* flag = new bool[nNode];	// tag array, judge is processed.
		for (i = 0; i < nNode; i++)
			flag[i] = false;
		dist[start_position] = 0;		// the self length of start position is zero
		while (true) {
			int v = -1; // initial unknow
			for (i = 0; i != nNode; i++)
				if (!flag[i] && dist[i] >= 0) // search the unprocessed node
												// with the shortest distance
					if (v < 0 || dist[i] < dist[v])
						v = i;
			if (v < 0)
				return; // if all nodes is processed, return
			flag[v] = true; // tag
			for (i = 0; i != nNode; ++i)
				if (adjmap[v][i] >= 0)
					if (dist[i] < 0 || dist[v] + adjmap[v][i] < dist[i]) {
						dist[i] = dist[v] + adjmap[v][i];
						path[i] = v;
					}
		}

		delete[] flag;
		flag = NULL;
	}

};

#endif