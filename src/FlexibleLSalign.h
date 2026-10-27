#ifndef FLEXIBLELSALIGN_H
#define FLEXIBLELSALIGN_H

#include "RigidLSalign.h"
#include "LigRotamerDetector.h"
#include "Ligands.h"
#include "Point3D.h"
using namespace std;

class FlexibleLSalign{
private :
	Ligand query;
	Ligand templ;
	
	double d0;		// default is 2.0

	bool isNeedMoreCorrectAlign;

	RigidLSalign* p_initalRLSalign; 
	vector<Ligand> queryRotamers;
public : 
	FlexibleLSalign(Ligand query, Ligand templ){
		this->query = query;
		this->templ = templ;
		this->d0 = 2.0;

		p_initalRLSalign = new RigidLSalign(query, templ, d0);

		LigRotamerDetector ligRD(query);
		queryRotamers = ligRD.generateRotamers();

		isNeedMoreCorrectAlign = false;
	}

	FlexibleLSalign(Ligand query, Ligand templ, double d0){
		this->query = query;
		this->templ = templ;
		this->d0 = d0;
		p_initalRLSalign = new RigidLSalign(query, templ, d0);
		LigRotamerDetector ligRD(query);
		queryRotamers = ligRD.generateRotamers();

		isNeedMoreCorrectAlign = false;
	}

	FlexibleLSalign(Ligand query, Ligand templ, double d0, 
		bool isNeedMoreSmallAngleStep, bool isWantedSaveTimeUsingRandomSelect,
		bool isNeedMoreCorrectAlign){
		this->query = query;
		this->templ = templ;
		this->d0 = d0;
		p_initalRLSalign = new RigidLSalign(query, templ, d0);
		LigRotamerDetector ligRD(query, isNeedMoreSmallAngleStep, isWantedSaveTimeUsingRandomSelect);
//		queryRotamers = ligRD.generateRotamers();
		queryRotamers = ligRD.generateRotamers(p_initalRLSalign);

		//////////////////////////////////////////////////////////////
		// Temp codes
//		string saveRotamersFolder = "/nfs/amino-home/hujunum/LSalign/0.NewWorkFolder/rotamers/DUDE/";
//		for (int i = 0; i < queryRotamers.size() && i < 100; i++){
//			char buf[500];
//			sprintf(buf, "%s_%d.pdb", query.getLigname().c_str(), i+1);
//			string savepath = saveRotamersFolder + buf;
//			queryRotamers[i].save8PDBFormat(savepath);
//		}
		//////////////////////////////////////////////////////////////


		this->isNeedMoreCorrectAlign = isNeedMoreCorrectAlign;
	}

	FlexibleLSalign(Ligand query, Ligand templ, RigidLSalign initalRLSalign){
		this->query = query;
		this->templ = templ;

		this->d0 = initalRLSalign.get_d0();

		this->p_initalRLSalign = new RigidLSalign(initalRLSalign);
		
		LigRotamerDetector ligRD(query);
		queryRotamers = ligRD.generateRotamers();

		isNeedMoreCorrectAlign = false;
	}

	~FlexibleLSalign(){
		if (NULL != p_initalRLSalign){
			delete p_initalRLSalign;
			p_initalRLSalign = NULL;
		}

		if (0 != queryRotamers.size()){
			vector<Ligand>().swap(queryRotamers);
		}
	}

	RigidLSalign* getInitialRigidLSalign(){
		return p_initalRLSalign;
	}

	RigidLSalign* align(int select_rotamer_num){
		int i = 0/*, j = 0*/;
		vector<Ligand> select_rotamers = selectRotamers(select_rotamer_num);
		
		RigidLSalign* p_ans = new RigidLSalign(*p_initalRLSalign);
		for (i = 0; i < select_rotamers.size(); i++){
			RigidLSalign* p_ithLSalign = new RigidLSalign(select_rotamers[i], templ, d0);
			if (p_ithLSalign->getLSscore() > p_ans->getLSscore()){
				if (NULL != p_ans){
					delete p_ans;
					p_ans = NULL;
				}

				p_ans = p_ithLSalign;
				p_ithLSalign = NULL;
			}
		}

		return p_ans;
	}
private :
	/*************************************************
	 * @param sn : the select rotamer number
     *************************************************/
	vector<Ligand> selectRotamers(int sn){
		int i = 0, j = 0;
		
		vector<Ligand> ans;
		if (queryRotamers.size() == 0){
			return ans;
		}
		
		if (queryRotamers.size() <= sn){
			for (i = 0; i < queryRotamers.size(); i++){
				ans.push_back(queryRotamers[i]);
			}

			return ans;
		}
		
		// across the init RigidLSalign information 
		vector<AlignedPair> initAlign = p_initalRLSalign->getFinalAlignedPair();
		vector<Point3D> templInitAlignPts;
		for (i = 0; i < initAlign.size(); i++){
			int tInd = initAlign[i].colInd;
			templInitAlignPts.push_back(templ.getIthPoints(tInd));
		}

		// 使用旋转后的query rotamer坐标计算原有align位置的旋转平移矩阵，然后计算RMSD_lb
		vector<double> rmsd_lb_arr;
		for (i = 0; i < queryRotamers.size(); i++){
			vector<Point3D> qrotamerInitAlignPts;
			for (j = 0; j < initAlign.size(); j++){
				int qInd = initAlign[j].rowInd;
				qrotamerInitAlignPts.push_back(queryRotamers[i].getIthPoints(qInd));
			}

			// using Kabsch's rotation matrix to superpose query to templ
			double u[3][3];
			double *t = new double[3];
			Kabsch(qrotamerInitAlignPts, templInitAlignPts, t, u);
			
			vector<Point3D> qRotamerPts_T = rotateAndTrans(queryRotamers[i].getPoints(), t, u);
			double rmsd_lb = RMSD_lb(qRotamerPts_T, templ.getPoints());
			
			rmsd_lb_arr.push_back(rmsd_lb);
			delete[] t;
			t = NULL;
		}

		bool* isSelected = new bool[queryRotamers.size()];
		for (i = 0; i < queryRotamers.size(); i++)
			isSelected[i] = false;
		int* index1 = ascentInsertSortIndex(rmsd_lb_arr, sn/2);
		for (i = 0; i < sn/2; i++){
			ans.push_back(queryRotamers[index1[i]]);
			isSelected[index1[i]] = true;
		}

		// 首先在query和templ原子数目少的ligand上计算距离两个距离最远的点和一个距离最远平均点
		// 再在原子数目多的ligand上计算最佳匹配的三个点，计算相似
		vector<double> dot_Rela_arr;

		double** p_tDisPow2Mtx = calculateDisPow2Mtx(templ.getPoints());
		if (query.size() > templ.size()){
			double* tSpecThreeDisPow2 = calculateThreeSpecicalDistPow2(templ.size(), p_tDisPow2Mtx);	
			for (i = 0; i < queryRotamers.size(); i++){
				if (isSelected[i]){
					dot_Rela_arr.push_back(9999999999999999.0);
					continue;
				}


				double** p_ithQDisPow2Mtx = calculateDisPow2Mtx(queryRotamers[i].getPoints());

				double sco = scoreWithThreeSpecialPts(queryRotamers[i].size(), p_ithQDisPow2Mtx,
					tSpecThreeDisPow2[0], tSpecThreeDisPow2[1], tSpecThreeDisPow2[2]);		// the sco is the lower the better

				dot_Rela_arr.push_back(sco);

				// release memory
				release2DArr(queryRotamers[i].size(), p_ithQDisPow2Mtx);
			}

			delete[] tSpecThreeDisPow2;
			tSpecThreeDisPow2 = NULL;
		}else{
			for (i = 0; i < queryRotamers.size(); i++){
				if (isSelected[i]){
					dot_Rela_arr.push_back(9999999999999999.0);
					continue;
				}

				double** p_ithQDisPow2Mtx = calculateDisPow2Mtx(queryRotamers[i].getPoints());
				double* qSpecThreeDisPow2 = calculateThreeSpecicalDistPow2(queryRotamers[i].size(), p_ithQDisPow2Mtx);

				double sco = scoreWithThreeSpecialPts(templ.size(), p_tDisPow2Mtx,
					qSpecThreeDisPow2[0], qSpecThreeDisPow2[1], qSpecThreeDisPow2[2]);

				dot_Rela_arr.push_back(sco);

				// release memory
				release2DArr(queryRotamers[i].size(), p_ithQDisPow2Mtx);
				delete[] qSpecThreeDisPow2;
				qSpecThreeDisPow2 = NULL;
			}

		}
		
		int* index2 = ascentInsertSortIndex(dot_Rela_arr, sn/2);
		for (i = 0, j = sn/2; j < sn; i++, j++){
			ans.push_back(queryRotamers[index2[i]]);
			isSelected[index2[i]] = true;
		}
		
		if (isNeedMoreCorrectAlign){
			//--------------------------------------------------------------------------//
			//-----------------     Time consuming using Raw RigidLSalign      ---------//
			//----------------   使用粗糙的RigidAlign去选择2个最佳的rotamers   ---------//
			//--------------------------------------------------------------------------//
			vector<double> rawRigidLSscoreVec;
			for (i = 0; i < queryRotamers.size(); i++){
				if (isSelected[i]){
					rawRigidLSscoreVec.push_back(-999999.9);
					continue;
				}
				
				int rough_sliding_step = 8;
				int detail_sliding_step = 4;
				RigidLSalign* p_rawRLSali = new RigidLSalign(queryRotamers[i], templ, d0, 
					rough_sliding_step, detail_sliding_step);

				double LS_score = p_rawRLSali->getLSscore();
				rawRigidLSscoreVec.push_back(LS_score);

				delete p_rawRLSali;
				p_rawRLSali = NULL;
			}
				
			int* index3 = descentInsertSortIndex(rawRigidLSscoreVec, 2);
			for (i = 0; i < 2; i++){
				ans.push_back(queryRotamers[index3[i]]);
			}

			delete[] index3;
			index3 = NULL;
		}
		
		// release memory 
		delete[] index2;
		index2 = NULL;
		release2DArr(templ.size(), p_tDisPow2Mtx);

		delete[] isSelected;
		isSelected = NULL;
		delete[] index1;
		index1 = NULL;

		return ans;
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
};

#endif 
