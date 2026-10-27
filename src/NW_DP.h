#ifndef NW_DP_H
#define NW_DP_H

int* runNWDP(double** scores, int row_num, int col_num, double gap_cost){
	//NW dynamic programming for alignment
	//not a standard implementation of NW algorithm
	//Input: vectors x, y, rotation matrix t, u, scale factor d02, and gap_cost
	//Output: j2i[1:col_num] \in {1:row_num} U {-1}
	//path[0:row_num, 0:col_num]=1,2,3, from diagonal, horizontal, vertical
	int i, j;
	double h, v, d;

	int* ans = new int[row_num+1];
	double** val = new double*[row_num+1];
	bool** path = new bool*[row_num+1];
	for (i = 0; i < row_num+1; i++){
		val[i] = new double[col_num+1];
		path[i] = new bool[col_num+1];
		for (j = 0; j < col_num+1; j++){
			val[i][j] = 0.0;
			path[i][j] = false;
		}
	}			

	//initialization
	val[0][0]=0;
	for(i=0; i<=row_num; i++){
		val[i][0]=0;
		path[i][0]=false; //not from diagonal
		ans[i]=-1;	//all are not aligned
	}

	for(j=0; j<=col_num; j++){
		val[0][j]=0;
		path[0][j]=false; //not from diagonal
	}      
	double xx[3], dij;

	//decide matrix and path
	for(i=1; i<=row_num; i++){	
		for(j=1; j<=col_num; j++){
			d=val[i-1][j-1] +  scores[i-1][j-1];

			//symbol insertion in horizontal (= a gap in vertical)
			h=val[i-1][j];
			if(path[i-1][j]) //aligned in last position
				h += gap_cost;				

			//symbol insertion in vertical
			v=val[i][j-1];
			if(path[i][j-1]) //aligned in last position
				v += gap_cost;

			if(d>=h && d>=v){
				path[i][j]=true; //from diagonal
				val[i][j]=d;
			}else{
				path[i][j]=false; //from horizontal
				if(v>=h)
					val[i][j]=v;
				else					
					val[i][j]=h;
			}
		} //for i
	} //for j

	//trace back to extract the alignment
	i=row_num;
	j=col_num;
	while(i>0 && j>0){
		if(path[i][j]){ //from diagonal
			ans[i-1] = j - 1;		
			i--;
			j--;
		}else{
			h=val[i-1][j];
			if(path[i-1][j]) h +=gap_cost;
			v=val[i][j-1];
			if(path[i][j-1]) v +=gap_cost;

			if(v>=h)
				j--;
			else
				i--;
		}
	}

	// release memory
	delete[] val;
	delete[] path;
	
	return ans;	
}

#endif
