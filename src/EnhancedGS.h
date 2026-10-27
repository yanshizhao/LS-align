#ifndef ENHACEDGS_H
#define ENHACEDGS_H

#include <vector>
#include "Sort.h"
using namespace std;

/*************************************************************
 * define a struct to save each aligned atom pair indexes
 *************************************************************/
typedef struct{
	int rowInd;		// query index
	int colInd;		// templ index
}AlignedPair;

class AlignInfo{
private:
	vector<AlignedPair> aliPairs;	// temp
	int* alivec; 	// the index is the row index, value means the col index, alivec[i] = -1 means the i-th element have no alignment
public :
	double sco;	
	int len;
	
	AlignInfo(double sco, const int* alivec, int len){
		this->sco = sco;
		this->len = len;
		this->alivec = new int[len];
		for (int i = 0; i < len; i++)
			this->alivec[i] = alivec[i];
	}
	
	AlignInfo(const AlignInfo& aliinfo){
		int i = 0;
		this->sco = aliinfo.sco;
		this->len = aliinfo.len;
		this->alivec = new int[len];
		for (i = 0; i < len; i++)
			this->alivec[i] = aliinfo.alivec[i];
		
		for (i = 0; i < aliinfo.aliPairs.size(); i++)
			this->aliPairs.push_back(aliinfo.aliPairs[i]);
	}
	
	int getIthAliVec(int ind){
		return alivec[ind];
	}
	
	int* getAliVec(){
		return alivec;
	}
	
	~AlignInfo(){
		if (NULL != alivec){
			delete[] alivec;
			alivec = NULL;
			
			len = 0;
		}
	}
	
	vector<AlignedPair> getAliPairs(){	// the aliInfo and alivec can exchange each other.
		if (0 == aliPairs.size()){
			for (int i = 0; i < len; i++){
				if (-1 != alivec[i]){
					AlignedPair alignedPair;
					alignedPair.rowInd = i;
					alignedPair.colInd = alivec[i];
	
					aliPairs.push_back(alignedPair);
				}
			}	
		}
		
		return aliPairs;
	}
};

/****************************************************
 * The Enhanced Greedy Search alogrithm can be found in the following paper
 * 	Brown, Peter, et al. "Fast and accurate non-sequential protein 
 *	structure alignment using a new asymmetric linear sum assignment
 *	 heuristic." Bioinformatics (2015): btv580. 
 ****************************************************/
class EnhancedGS{
private:
	vector<vector<double> >	scoMtx;
	vector<int*> rowSortedInds;		// temp save the sort row content, the vector size is equal to the row of scoMtx
	AlignInfo* pAli;			// the align information in scoMtx
public:
	EnhancedGS(vector<vector<double> > scoMtx){
		this->scoMtx = scoMtx;
		
		pAli = enhanceGreedySeach();
//		pAli = LoopEnhanceGreedySearch();
	}

	EnhancedGS(vector<vector<double> > scoMtx, int type){
		this->scoMtx = scoMtx;
		
		if (0 == type){
			pAli = greedySearch();			
		}else if (1 == type){
			pAli = enhanceGreedySeach();
		}else{
			pAli = LoopEnhanceGreedySearch();
		}
		
	}

	~EnhancedGS(){
		if (0 != rowSortedInds.size()){
			for (int i = 0; i < rowSortedInds.size(); i++){
				delete(rowSortedInds[i]);
				rowSortedInds[i] = NULL;
			}
			rowSortedInds.clear();
			vector<int*>().swap(rowSortedInds);
		}
		
		if (0 != scoMtx.size()){
			for (int i = 0; i < scoMtx.size(); i++){
				scoMtx[i].clear();
				vector<double>().swap(scoMtx[i]);
			}
			scoMtx.clear();
			vector<vector<double> >().swap(scoMtx);
		}

		delete pAli;
		pAli = NULL;
	}

	AlignInfo* getPAliInfo(){
		return pAli;
	}

private :

	AlignInfo* LoopEnhanceGreedySearch(){
		int i = 0;

		// the original greedy search procedure
		AlignInfo* gsPAli = greedySearch();

		int rownum = scoMtx.size();
		int colnum = scoMtx[0].size();
		int* alivec = new int[rownum];	// new memory
		for (i = 0; i < rownum; i++)
			alivec[i] = gsPAli->getIthAliVec(i);
		
		double oScore = gsPAli->sco;

		if (isTheHighestValuesOfAllRowsSelected(alivec, rownum, colnum)){
			return gsPAli;
		}

		// initial
		double newScore = oScore;
		int* newAlivec = new int[rownum];
		for (i = 0; i < rownum; i++)
			newAlivec[i] = alivec[i];

		int max_iter = rownum>colnum ? colnum : rownum;
		max_iter = max_iter>1 ? max_iter : 1;
		int iter = 0;
		while (iter++ < max_iter){
//			cout << endl << iter << "sco " << newScore << endl;
//			for (i = 0; i < rownum; i++)
//				cout << iter << " : " << i << " : " << alivec[i] << endl;			

			// then we want to enhance the gs alignment
			// sorted the row indexes based on the conditional max values from each rows
			vector<double> rowedCondMaxVals;
			for (i = 0; i < rownum; i++){
				if (-1 != newAlivec[i]){
					rowedCondMaxVals.push_back(scoMtx[i][newAlivec[i]]);
				}else{
					rowedCondMaxVals.push_back(0.0);
				}
			}
			int* descendInds = descentInsertSortIndex(rowedCondMaxVals);	// new memory

			for (i = 0; i < rownum; i++){
				int row_ind = descendInds[i];
				int col_ind = newAlivec[row_ind];

				double maxScoreChange = 0.0;
				int exchange_row_ind = -1;
				int exchange_col_ind = -1;
				for (int j = 0; j < colnum; j++){
					int new_col_ind = rowSortedInds[row_ind][j];

					// calculate the increase value 
					double incVal = scoMtx[row_ind][new_col_ind];
					if (-1 != col_ind) 
						incVal = scoMtx[row_ind][new_col_ind] - scoMtx[row_ind][col_ind];
				
					// calculate the decrease value
					int rela_another_row_ind = -1;
					for (int k = 0; k < rownum; k++){
						if (newAlivec[k] == new_col_ind){
							rela_another_row_ind = k;
						}
					}

					double decVal = 0.0;
					if (-1 != rela_another_row_ind)
						decVal = scoMtx[rela_another_row_ind][new_col_ind];
					if (-1 != col_ind && -1 != rela_another_row_ind)
						decVal = scoMtx[rela_another_row_ind][new_col_ind] - scoMtx[rela_another_row_ind][col_ind];

					if (incVal - decVal > maxScoreChange){
						maxScoreChange = incVal - decVal;
						exchange_row_ind = rela_another_row_ind;
						exchange_col_ind = new_col_ind;
					}
				}

				if (maxScoreChange > 0.0){
					newScore += maxScoreChange;
					newAlivec[row_ind] = exchange_col_ind;
					if (-1 != exchange_row_ind)
						newAlivec[exchange_row_ind] = col_ind;
				}
			}
			delete[] descendInds;
			descendInds = NULL;		

			if (newScore > oScore){
				oScore = newScore;
				for (i = 0; i < rownum; i++)
					alivec[i] = newAlivec[i];
			}else if (newScore < oScore){
				cout << newScore << "\t" << oScore << endl;  
				cout << "You have a wrong in testLoopEnhanceGreedySearch" << endl;
			}else {
				break;
			}
		}
	
		AlignInfo* pAns = new AlignInfo(oScore, alivec, rownum);
		
		delete gsPAli;
		gsPAli = NULL;
		delete[] alivec;
		alivec = NULL;
		delete[] newAlivec;
		newAlivec = NULL;
		return pAns;
	}	

	AlignInfo* enhanceGreedySeach(){
		int i = 0;

		// the original greedy search procedure
		AlignInfo* gsPAli = greedySearch();

		int rownum = scoMtx.size();
		int colnum = scoMtx[0].size();	
		int* alivec = new int[rownum];	// new memory
		for (i = 0; i < rownum; i++)
			alivec[i] = gsPAli->getIthAliVec(i);
		
		double oScore = gsPAli->sco;

		if (isTheHighestValuesOfAllRowsSelected(alivec, rownum, colnum)){
			return gsPAli;
		}

		
		// then we want to enhance the gs alignment
		// sorted the row indexes based on the conditional max values from each rows
		vector<double> rowedCondMaxVals;
		for (i = 0; i < rownum; i++){
			if (-1 != alivec[i]){
				rowedCondMaxVals.push_back(scoMtx[i][alivec[i]]);
			}else{
				rowedCondMaxVals.push_back(0.0);
			}
		}
		
		int* descendInds = descentInsertSortIndex(rowedCondMaxVals);	// new memory

		for (i = 0; i < rownum; i++){
			int row_ind = descendInds[i];
			int col_ind = alivec[row_ind];

			double maxScoreChange = 0.0;
			int exchange_row_ind = -1;
			int exchange_col_ind = -1;
			for (int j = 0; j < colnum; j++){
				int new_col_ind = rowSortedInds[row_ind][j];
				if (new_col_ind == col_ind){
					break;
				}

				// calculate the increase value 
				double incVal = scoMtx[row_ind][new_col_ind];
				if (-1 != col_ind) 
					incVal = scoMtx[row_ind][new_col_ind] - scoMtx[row_ind][col_ind];
				if (incVal < 0){
					cout << "scoMtx[row_ind][new_col_ind] - scoMtx[row_ind][col_ind]; has a problem" << endl;
				}
				
				// calculate the decrease value
				int rela_another_row_ind = -1;
				for (int k = 0; k < rownum; k++){
					if (alivec[k] == new_col_ind){
						rela_another_row_ind = k;
					}
				}

				double decVal = 0.0;
				if (-1 != rela_another_row_ind)
					decVal = scoMtx[rela_another_row_ind][new_col_ind];
				if (-1 != col_ind && -1 != rela_another_row_ind)
					decVal = scoMtx[rela_another_row_ind][new_col_ind] - scoMtx[rela_another_row_ind][col_ind];

				if (incVal - decVal > maxScoreChange){
					maxScoreChange = incVal - decVal;
					exchange_row_ind = rela_another_row_ind;
					exchange_col_ind = new_col_ind;
				}
			}
			
			if (maxScoreChange > 0.0 && exchange_row_ind < 0){
				cout << " ---------------88-------------------- " << endl;
			}

			if (maxScoreChange > 0.0){
				oScore += maxScoreChange;
				alivec[row_ind] = exchange_col_ind;
				alivec[exchange_row_ind] = col_ind;
			}
		}
		delete[] descendInds;
		descendInds = NULL;
	
		AlignInfo* pAns = new AlignInfo(oScore, alivec, rownum);
		
		delete gsPAli;
		gsPAli = NULL;
		delete[] alivec;
		alivec = NULL;
		return pAns;
	}

	AlignInfo* greedySearch(){
		int i = 0; 
		int rownum = scoMtx.size();
		int colnum = scoMtx[0].size();

		// sort each row content index
		for (i = 0; i < rownum; i++){
			int* tmp = descentInsertSortIndex(scoMtx[i]);
			rowSortedInds.push_back(tmp); 
		}		

		int* alivec = new int[rownum];
		for (i = 0; i < rownum; i++)
			alivec[i] = -1;
		bool* isColVisited = new bool[colnum];
		for (i = 0; i < colnum; i++)
			isColVisited[i] = false;
		double score = 0.0;
		for (i = 0; i < rownum; i++){	// i means the ith number, not row
			double maxValue = 0.0;
			int corrRowInd = -1;
			int corrColInd = -1;
			for (int rr = 0; rr < rownum; rr++){	// rr means the rr-th row in rowSortedInds
				if (-1 == alivec[rr]){ // the rr-th is not aligned
					int colInd = -1;
					double value = 0.0;
					for (int j = 0; j < colnum; j++){ // j means the j-th number not column
						if (false == isColVisited[rowSortedInds[rr][j]]){	// rowSortedInds[rr][j] is 
							colInd = rowSortedInds[rr][j];
							value = scoMtx[rr][colInd];
							
							break;
						}
					}
					
					if (maxValue < value){
						maxValue = value;
						corrRowInd = rr;
						corrColInd = colInd;
					}
				}
			}
			
			if (-1 == corrRowInd || -1 == corrColInd){
				continue;
			}
			
			score += maxValue;
			alivec[corrRowInd] = corrColInd;
			isColVisited[corrColInd] = true;

		}

		AlignInfo* pAns = new AlignInfo(score, alivec, rownum);
		
		delete[] isColVisited;
		delete[] alivec;
		alivec = NULL;
		return pAns;
	}

	bool isTheHighestValuesOfAllRowsSelected(int* alivec, int rownum, int colnum){
		bool ans = true;
		bool isBreak = false;
		int minusOneNum = 0;
		for (int i = 0; i < rownum; i++){
			int r = i;
			int c = alivec[i];
			
			if (c != -1 && c != rowSortedInds[r][0]){
				ans = false;
				isBreak = true;
				break;
			}
			
			if (c == -1){
				minusOneNum++;
			}
		}
		
		if (!isBreak
			&& rownum > colnum 
			&& rownum - colnum == minusOneNum){
			ans = true;
		}
		
		return ans;
	}
};

#endif
