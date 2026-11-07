#ifndef LIG_ROTAMER_DETECTOR_H
#define LIG_ROTAMER_DETECTOR_H

#include "ShortestRoad.h"
#include "basic_fun.h"
#include "Ligands.h"
#include <vector>
using namespace std;

#define DOUBLE_PI 6.283185307179586	// 2 * 3.141592653589793

#define NRS_ROTATABLE_BOND_NUM 3
#define NRS_MIN_ROTATABLE_BOND_NUM 3
#define NRS_LEGACY_UNPARAM_BOND_NUM 5	// differs from the paper (3); see R5
#define ROTATABLE_BOND_MIN_SIDE_ATOMS 4

class LigRotamerDetector{
private:
	Ligand lig;
	ShortestRoad* p_sroad;
	
	// tmp
	double** origLigDisMtx;		
	bool isWantedSaveTimeUsingRandomSelect;	// true: cap rotamer sets > 5000 by deterministic decimation (X1)

	// const final 
	vector<double> rotAngles;
	
public:
	LigRotamerDetector(Ligand lig){
		int i = 0, j = 0;
		this->lig = lig;	                    
		p_sroad = shortestRoadBetweenTwoAtoms();
			
		origLigDisMtx = new double*[lig.size()];
		for (i = 0; i < lig.size(); i++){
			origLigDisMtx[i] = new double[lig.size()];
			for (j = 0; j < lig.size(); j++){
				double dis = 0.0;
				if (i == j) dis = 0.0;
				else if (j > i) dis = lig.getIthPoints(i).dist(lig.getIthPoints(j));
				else dis = origLigDisMtx[j][i];
				
				origLigDisMtx[i][j] = dis;
			}
		}

		rotAngles.push_back(0.0);
		rotAngles.push_back(-120.0/360.0 * DOUBLE_PI);
		rotAngles.push_back(+120.0/360.0 * DOUBLE_PI);
		rotAngles.push_back(-60.0/360.0 * DOUBLE_PI);
		rotAngles.push_back(+60.0/360.0 * DOUBLE_PI);

		isWantedSaveTimeUsingRandomSelect = false;
	}	
	
	LigRotamerDetector(Ligand lig, bool isNeedMoreSmallAngleStep){
		int i = 0, j = 0;
		this->lig = lig;	                    
		p_sroad = shortestRoadBetweenTwoAtoms();
			
		origLigDisMtx = new double*[lig.size()];
		for (i = 0; i < lig.size(); i++){
			origLigDisMtx[i] = new double[lig.size()];
			for (j = 0; j < lig.size(); j++){
				double dis = 0.0;
				if (i == j) dis = 0.0;
				else if (j > i) dis = lig.getIthPoints(i).dist(lig.getIthPoints(j));
				else dis = origLigDisMtx[j][i];
				
				origLigDisMtx[i][j] = dis;
			}
		}

		rotAngles.push_back(0.0);
		rotAngles.push_back(-120.0/360.0 * DOUBLE_PI);
		rotAngles.push_back(+120.0/360.0 * DOUBLE_PI);
		rotAngles.push_back(-60.0/360.0 * DOUBLE_PI);
		rotAngles.push_back(+60.0/360.0 * DOUBLE_PI);

		if (isNeedMoreSmallAngleStep){
			rotAngles.push_back(-180.0/360.0 * DOUBLE_PI);
			rotAngles.push_back(+180.0/360.0 * DOUBLE_PI);
			rotAngles.push_back(-90.0/360.0 * DOUBLE_PI);
			rotAngles.push_back(+90.0/360.0 * DOUBLE_PI);
			rotAngles.push_back(-30.0/360.0 * DOUBLE_PI);
			rotAngles.push_back(+30.0/360.0 * DOUBLE_PI);
		}

		isWantedSaveTimeUsingRandomSelect = false;
	}

	LigRotamerDetector(Ligand lig, bool isNeedMoreSmallAngleStep, bool isWantedSaveTimeUsingRandomSelect){
		int i = 0, j = 0;
		this->lig = lig;	                    
		p_sroad = shortestRoadBetweenTwoAtoms();
			
		origLigDisMtx = new double*[lig.size()];
		for (i = 0; i < lig.size(); i++){
			origLigDisMtx[i] = new double[lig.size()];
			for (j = 0; j < lig.size(); j++){
				double dis = 0.0;
				if (i == j) dis = 0.0;
				else if (j > i) dis = lig.getIthPoints(i).dist(lig.getIthPoints(j));
				else dis = origLigDisMtx[j][i];
				
				origLigDisMtx[i][j] = dis;
			}
		}

		rotAngles.push_back(0.0);

		rotAngles.push_back(-120.0/360.0 * DOUBLE_PI);
		rotAngles.push_back(+120.0/360.0 * DOUBLE_PI);
		rotAngles.push_back(-60.0/360.0 * DOUBLE_PI);
		rotAngles.push_back(+60.0/360.0 * DOUBLE_PI);

		if (isNeedMoreSmallAngleStep){
			rotAngles.push_back(-180.0/360.0 * DOUBLE_PI);
			rotAngles.push_back(+180.0/360.0 * DOUBLE_PI);
			rotAngles.push_back(-90.0/360.0 * DOUBLE_PI);
			rotAngles.push_back(+90.0/360.0 * DOUBLE_PI);
			rotAngles.push_back(-30.0/360.0 * DOUBLE_PI);
			rotAngles.push_back(+30.0/360.0 * DOUBLE_PI);
		}

		this->isWantedSaveTimeUsingRandomSelect = isWantedSaveTimeUsingRandomSelect;
	}

	~LigRotamerDetector(){
		int i;

		if (NULL != origLigDisMtx){
			for (i = 0; i < lig.size(); i++){
				delete[] origLigDisMtx[i];
			}
			delete[] origLigDisMtx;
			origLigDisMtx = NULL;
		}

		if (NULL != p_sroad){
			delete p_sroad;
			p_sroad = NULL;
		}
		
		vector<double>().swap(rotAngles);
	}
	
	double** getOrigLigDisMtx(){
		return origLigDisMtx;
	}

	/****************************************************
	 * rotate the rotatable single bond to generate the 
	 * available rotamers.
	 ****************************************************/
	vector<Ligand> generateRotamers(){
		vector<Ligand> ans;
		vector<BondInfo> singleBonds = findAllSingleBond();
		if (0 == singleBonds.size()){
			return ans;
		}

		vector<vector<vector<int> > > twoPartsIndss = getUsefulBondTwoPartIndsInfo(&singleBonds);
		if (0 == twoPartsIndss.size()){
			return ans;
		}

		// select 5 avaliable single bond
		vector<int> selectNIndsOfSignleBonds = select_N_RotableSingleBondOnUniformDistribution(twoPartsIndss, singleBonds, NRS_LEGACY_UNPARAM_BOND_NUM);
		
		int all_rotamer_num = (int)pow(rotAngles.size(), selectNIndsOfSignleBonds.size());

		int max_iter_num = all_rotamer_num;
		if (isWantedSaveTimeUsingRandomSelect && all_rotamer_num > 5000){
			max_iter_num = 5000;
		}

		int* mod_ref_arr = new int[selectNIndsOfSignleBonds.size()];
		for (int i = 0; i < selectNIndsOfSignleBonds.size(); i++){
			mod_ref_arr[i] = (int)pow(rotAngles.size(), (int)selectNIndsOfSignleBonds.size()-1-i);
		}

		for (int iter = 1; iter < max_iter_num; iter++){			// the rotamer set does not include itself.
			int i = 0;
			
			int recordNum = iter;
			if (isWantedSaveTimeUsingRandomSelect && all_rotamer_num > 5000){
				// X1: deterministic decimation (evenly spaced, reproducible)
				recordNum = (int)((double)iter * all_rotamer_num / max_iter_num);
				if (recordNum >= all_rotamer_num) recordNum = all_rotamer_num - 1;
			}

			double* angles = new double[selectNIndsOfSignleBonds.size()];
			for (i = 0; i < selectNIndsOfSignleBonds.size(); i++){                                                                         
				int mod_ref = mod_ref_arr[i];
				int iInd = recordNum / mod_ref;                                                                                        
				angles[i] = rotAngles[iInd];                                                                                           
				
				recordNum -= iInd*mod_ref;                                                                                             
			}
			
			// rotated selected single bond with different angles;
			Ligand* p_decoy = new Ligand(lig);
			for (i = 0; i < selectNIndsOfSignleBonds.size(); i++){
				BondInfo bi = singleBonds[selectNIndsOfSignleBonds[i]];
				vector<vector<int> > twoPartInds = twoPartsIndss[selectNIndsOfSignleBonds[i]];
				
				double angle = angles[i];                                                                                              
				if (angle == 0.0)
					continue;

				rotatedLigandOnePartOfOneBond(p_decoy, bi, twoPartInds, angle);
			}
			
			if (!judgeIsOverlapAtomsInLig(p_decoy)){
				ans.push_back(*p_decoy);
			}
			delete p_decoy;
			p_decoy = NULL;

			delete[] angles;
			angles = NULL;
		}		

		delete[] mod_ref_arr;
		mod_ref_arr = NULL;
		return ans;
	}

	/****************************************************
	 * rotate the rotatable single bond (which is not aliged
	 * in pRigidLSalign) to generate the available rotamers.
	 ****************************************************/
	vector<Ligand> generateRotamers(RigidLSalign* pRigidLSalign){
		vector<Ligand> ans;
		vector<BondInfo> _singleBonds = findAllSingleBond();
		if (0 == _singleBonds.size()){
			return ans;
		}
		
		// Remove these BondInfoes which are aligned in RigidAlign
		vector<BondInfo> singleBonds;
		vector<AlignedPair> rigidAlignInfo = pRigidLSalign->getFinalAlignedPair();
		vector<double> rigidAlignDiss = pRigidLSalign->getFinalAlignedPairDist();
		for (int dan = 0; dan < _singleBonds.size(); dan++){
			int aInd = _singleBonds[dan].aInd;
			int bInd = _singleBonds[dan].bInd;
			bool isAIndAligned = false;
			bool isBIndAligned = false;
			double aveDis = 0.0;
			for (int jun = 0; jun < rigidAlignInfo.size(); jun++){
				AlignedPair ai = rigidAlignInfo[jun];
				
				if (!isAIndAligned && ai.rowInd == aInd){
					isAIndAligned = true;
					aveDis += rigidAlignDiss[jun];
				}

				if (!isBIndAligned && ai.rowInd == bInd){
					isAIndAligned = true;
					aveDis += rigidAlignDiss[jun];
				}

				if (isAIndAligned && isBIndAligned){
					break;
				}
			}
			
			aveDis /= 2.0;
			if (!(isAIndAligned && isBIndAligned 
				&& aveDis < 1.0)){
				singleBonds.push_back(_singleBonds[dan]);
			}
		}

		if (0 == singleBonds.size()){
			return ans;
		}

		vector<vector<vector<int> > > twoPartsIndss = getUsefulBondTwoPartIndsInfo(&singleBonds);
		if (0 == twoPartsIndss.size()){
			return ans;
		}

		// select 3 avaliable single bond
		vector<int> selectNIndsOfSignleBonds = select_N_RotableSingleBondOnUniformDistribution(twoPartsIndss, singleBonds, NRS_ROTATABLE_BOND_NUM);
		
		int all_rotamer_num = (int)pow(rotAngles.size(), selectNIndsOfSignleBonds.size());

		int max_iter_num = all_rotamer_num;
		if (isWantedSaveTimeUsingRandomSelect && all_rotamer_num > 5000){
			max_iter_num = 5000;
		}

		int* mod_ref_arr = new int[selectNIndsOfSignleBonds.size()];
		for (int i = 0; i < selectNIndsOfSignleBonds.size(); i++){
			mod_ref_arr[i] = (int)pow(rotAngles.size(), (int)selectNIndsOfSignleBonds.size()-1-i);
		}

		for (int iter = 1; iter < max_iter_num; iter++){			// the rotamer set does not include itself.
			int i = 0;
			
			int recordNum = iter;
			if (isWantedSaveTimeUsingRandomSelect && all_rotamer_num > 5000){
				// X1: deterministic decimation (evenly spaced, reproducible)
				recordNum = (int)((double)iter * all_rotamer_num / max_iter_num);
				if (recordNum >= all_rotamer_num) recordNum = all_rotamer_num - 1;
			}

			double* angles = new double[selectNIndsOfSignleBonds.size()];
			for (i = 0; i < selectNIndsOfSignleBonds.size(); i++){                                                                         
				int mod_ref = mod_ref_arr[i];
				int iInd = recordNum / mod_ref;                                                                                        
				angles[i] = rotAngles[iInd];                                                                                           
				
				recordNum -= iInd*mod_ref;                                                                                             
			}
			
			// rotated selected single bond with different angles;
			Ligand* p_decoy = new Ligand(lig);
			for (i = 0; i < selectNIndsOfSignleBonds.size(); i++){
				BondInfo bi = singleBonds[selectNIndsOfSignleBonds[i]];
				vector<vector<int> > twoPartInds = twoPartsIndss[selectNIndsOfSignleBonds[i]];
				
				double angle = angles[i];                                                                                              
				if (angle == 0.0)
					continue;

				rotatedLigandOnePartOfOneBond(p_decoy, bi, twoPartInds, angle);
			}
			
			if (!judgeIsOverlapAtomsInLig(p_decoy)){
				ans.push_back(*p_decoy);
			}
			delete p_decoy;
			p_decoy = NULL;

			delete[] angles;
			angles = NULL;
		}		

		delete[] mod_ref_arr;
		mod_ref_arr = NULL;
		return ans;
	}
	
private:
	
	// return the select index in singleBonds
	vector<int> select_N_RotableSingleBondOnUniformDistribution(vector<vector<vector<int> > > twoPartsIndss, 
		vector<BondInfo> singleBonds, int N){
		if (N < NRS_MIN_ROTATABLE_BOND_NUM){
			N = NRS_MIN_ROTATABLE_BOND_NUM;
		}

		int i = 0, j = 0;
		vector<int> ans;
		if (N >= singleBonds.size()){
			for (i = 0; i < singleBonds.size(); i++)
				ans.push_back(i);
			return ans;
		}

		vector<int> samllPartAtomNums; 
		vector<int> smallPartAtomInds;	// content just 0 / 1
		for (i = 0; i < singleBonds.size(); i++){
			if (i >= (int)twoPartsIndss.size() || twoPartsIndss[i].size() < 2){
				cerr << "ERROR: rotamer sampling: bond " << i
					<< " lacks two atom sides; aborting" << endl;
				return vector<int>();
			}
			int aPartAtomNum = twoPartsIndss[i][0].size();
			int bPartAtomNum = twoPartsIndss[i][1].size();
			
			int smallAtomNum = aPartAtomNum;
			int smallPartAtomInd = 0;
			if (aPartAtomNum > bPartAtomNum){
				smallAtomNum = bPartAtomNum;
				smallPartAtomInd = 1;
			}

			samllPartAtomNums.push_back(smallAtomNum);
			smallPartAtomInds.push_back(smallPartAtomInd);
		}

		int* index = ascentInsertSortIndex(samllPartAtomNums);
		if (samllPartAtomNums.size() < 2){
			cerr << "ERROR: rotamer sampling: fewer than 2 usable bonds; aborting" << endl;
			delete[] index;
			return vector<int>();
		}

		vector<int> sortedSIInds;
		sortedSIInds.push_back(index[0]);
		sortedSIInds.push_back(index[1]);

		for (i = 2; i < samllPartAtomNums.size(); i++){
			vector<int> ithSmallPartAtomInds = twoPartsIndss[index[i]][smallPartAtomInds[index[i]]];
			for (j = 0; j < i; j++){
				vector<int> jthSmallPartAtomInds = twoPartsIndss[index[j]][smallPartAtomInds[index[j]]];

				bool isOverlap = false;
				for (int k = 0; k < ithSmallPartAtomInds.size(); k++){
					for (int l = 0; l < jthSmallPartAtomInds.size(); l++){
						if (ithSmallPartAtomInds[k] == jthSmallPartAtomInds[l]){
							isOverlap = true;
							break;
						}
					}
					if (true == isOverlap){
						break;
					}
				}
				if (false == isOverlap)
					break;
			}
			
			sortedSIInds.insert(sortedSIInds.begin()+j, *(index+i));
		}

		int step = singleBonds.size() / N;
		int mod = singleBonds.size() % N;
		for (i = 0; i < N; i++){
			int pick = (i >= mod) ? (i*step + mod) : (i*step + i);
			ans.push_back(sortedSIInds[pick]);
		}

		delete[] index;
		index = NULL;

		return ans;
	}


	bool judgeIsOverlapAtomsInLig(Ligand* p_decoy){
		AtomVarDerWaalRadius avdwr;
		BondType** bondTypeMtx = p_decoy->getBondTypeMtx();
		
		bool isOverlap = false;
		for (int i = 0; i < p_decoy->size(); i++){
			Point3D iPos = p_decoy->getIthPoints(i);
			string iAtomtype = p_decoy->getAtomType(i);
			double iAtomRadii = avdwr[iAtomtype];
			
			for (int j = i+1; j < p_decoy->size(); j++){
				if (NOTCONNECTED != bondTypeMtx[i][j])
					continue;
							
				Point3D jPos = p_decoy->getIthPoints(j);
				string jAtomtype = p_decoy->getAtomType(j);
				double jAtomRadii = avdwr[jAtomtype];
					
				double dis = iPos.dist(jPos);
				if (dis < origLigDisMtx[i][j] 
					&& dis/((iAtomRadii+jAtomRadii > origLigDisMtx[i][j] ? origLigDisMtx[i][j] : iAtomRadii+jAtomRadii)) < 0.8){		// TODO here is a constant
					isOverlap = true;
					break;
				}
			}
			
			if (isOverlap){
				break;
			}
		}
		
		return isOverlap;
	}
	
	/*********************************************************************************************************
	 * @param p_lig: the content will be changed.
	 *********************************************************************************************************/
	void rotatedLigandOnePartOfOneBond(Ligand* p_lig, BondInfo bi, vector<vector<int> > twoPartInds, double angle){
		int i = 0;

		vector<int> bigPartInds;
		vector<int> smallPartInds;
		if (twoPartInds[0].size() > twoPartInds[1].size()){
			bigPartInds = twoPartInds[0];
			smallPartInds = twoPartInds[1];
		}else{
			bigPartInds = twoPartInds[1];                                                                                        
			smallPartInds = twoPartInds[0];                                                                                      
		}
		
		int bigPartAtomIndInBond = bi.aInd;                                                                           
		int smallPartAtomIndInBond = bi.bInd;                                                                              
		if (isContain(smallPartInds, bi.aInd)){                                                               
			smallPartAtomIndInBond =bi.aInd;                                                                                
			bigPartAtomIndInBond = bi.bInd;            
		}
		
		Point3D rot_axis = p_lig->getIthPoints(smallPartAtomIndInBond) - p_lig->getIthPoints(bigPartAtomIndInBond);
		if (rot_axis.dist(Point3D(0.0, 0.0, 0.0)) < 1e-12){
			cerr << "LSalign: WARNING: zero-length rotation axis (atoms "
				<< smallPartAtomIndInBond << " and " << bigPartAtomIndInBond
				<< " coincide); rotation skipped for this angle." << endl;
			return;
		}
		double** rotMtx = generateRotatedCoordinate(rot_axis, angle);			// need release
		
		// around the bigPartAtomIndInBond point to rotated the small part
		Point3D rotation_point = p_lig->getIthPoints(bigPartAtomIndInBond);

		vector<Point3D> rottedSmallPartPts;
		for (i = 0; i < smallPartInds.size(); i++){
			double cx = rotation_point[0];
			double cy = rotation_point[1];
			double cz = rotation_point[2];

			Point3D ithPt = p_lig->getIthPoints(smallPartInds[i]);
			double x = ithPt[0];
			double y = ithPt[1];
			double z = ithPt[2];
			
			double xt = cx + rotMtx[0][0]*(x-cx) + rotMtx[0][1]*(y-cy) + rotMtx[0][2]*(z-cz);
			double yt = cy + rotMtx[1][0]*(x-cx) + rotMtx[1][1]*(y-cy) + rotMtx[1][2]*(z-cz);
			double zt = cz + rotMtx[2][0]*(x-cx) + rotMtx[2][1]*(y-cy) + rotMtx[2][2]*(z-cz);
			
			Point3D pt(xt, yt, zt);
			rottedSmallPartPts.push_back(pt);
		}
		
		// load the ligand points information
		bool* isSet = new bool[p_lig->size()];
		for (i = 0; i < p_lig->size(); i++)
			isSet[i] = false;
		Point3D* newPoss = new Point3D[p_lig->size()];
		for (i = 0; i < bigPartInds.size(); i++){
			int ind = bigPartInds[i];
			newPoss[ind] = p_lig->getIthPoints(ind);
			isSet[ind] = true;
		}
		for (i = 0; i < smallPartInds.size(); i++){
			int ind = smallPartInds[i];                                                                                          
			newPoss[ind] = rottedSmallPartPts[i];                                                                                    
			isSet[ind] = true;                     
		}
		
		vector<Point3D> _poss;
		for (i = 0; i < p_lig->size(); i++){                                                                                  
			if (false == isSet[i]) {                                                                                                 
				cerr << "ERROR: 201702161138" << endl;
				// write to log
			} else _poss.push_back(newPoss[i]);                                                                                            
		}
		
		p_lig->setPoints(_poss);
		
		delete[] isSet;
		isSet = NULL;
		release2DArr(3, rotMtx);
		delete[] newPoss;
		newPoss = NULL;
	}
	
	ShortestRoad* shortestRoadBetweenTwoAtoms(){
		int node_num = lig.size();
		
		// construct the atoms adjacency matrix                                                                                    
		double** adjacencyMtx = new2DArr(node_num, node_num); // 0.0 means no adjacency                                          
		vector<BondInfo> bondAtomInfo = lig.getBondInfoes();                                                                       
		for (int i = 0; i < bondAtomInfo.size(); i++){                                                                              
			BondInfo ithBond = bondAtomInfo[i];                                                                               
			double ithBondDis = ithBond.bondDis;                                                                                          
			                                                                                                                         
			adjacencyMtx[ithBond.aInd][ithBond.bInd] = ithBondDis;                                                         
			adjacencyMtx[ithBond.bInd][ithBond.aInd] = ithBondDis;                                                         
		}                                                                                                                          
		                                                                                                                           
		ShortestRoad* p_sroad = new ShortestRoad(adjacencyMtx, node_num); 
		
		release2DArr(node_num, adjacencyMtx);
		// B31: report the rejection on stderr; callers see isValid() == false
		// and skip every bond through leftRightPartOnBond().
		if (NULL != p_sroad && !p_sroad->isValid()){
			int maxNodes = ShortestRoad::MAX_FLOYD_NODES;
			cerr << "ERROR: ligand '" << lig.getLigname() << "' has "
				<< node_num << " atoms, above the supported maximum of "
				<< maxNodes << "; skipped" << endl;
		}
		return p_sroad;
	}	
	
	vector<vector<vector<int> > > getUsefulBondTwoPartIndsInfo(vector<BondInfo>* p_singleBondInfos){
		int i = 0;

		vector<vector<vector<int> > > ans;
		vector<BondInfo> singleBIs;
		for (i = 0; i < p_singleBondInfos->size(); i++){
			BondInfo oneSB = (*p_singleBondInfos)[i];
			vector<vector<int> > twoPartInds = leftRightPartOnBond(oneSB);
			
			
			if (0 == twoPartInds.size()
					|| twoPartInds[0].size() < ROTATABLE_BOND_MIN_SIDE_ATOMS
					|| twoPartInds[1].size() < ROTATABLE_BOND_MIN_SIDE_ATOMS
					/*|| 1.0*twoPartInds[0].size()/ligsize <= 0.25                                                                       
					|| 1.0*twoPartInds[1].size()/ligsize <= 0.25*/){
				continue;
			}    
			
			ans.push_back(twoPartInds);
			singleBIs.push_back(oneSB);
		}
		
		p_singleBondInfos->clear();
		for (i = 0; i < singleBIs.size(); i++){
			p_singleBondInfos->push_back(singleBIs[i]);
		}
		
		return ans;
	}
	
	vector<vector<int> > leftRightPartOnBond(BondInfo onebond){ 
		int i = 0, j = 0, endInd = 0;

		// B31: an oversized ligand leaves ShortestRoad unusable; skip the bond.
		if (NULL == p_sroad || !p_sroad->isValid()){
			return vector<vector<int> >();
		}

		int node_num = lig.size();
		                                                  
		int*** allRoad = p_sroad->getAllRoad8Floyd();
		                                                                          
		int firstBondAtomInd = onebond.aInd;
		int secondBondAtomInd = onebond.bInd;
		                                                                                                                           
		vector<int> firstPartAtomInds;                                                                 
		bool* isFirstPart = new bool[node_num]; 
		for (i = 0; i < node_num; i++) isFirstPart[i] = false;                                                                          
		for (endInd = 0; endInd < node_num; endInd++){                                                                          
			if (isFirstPart[endInd])                                                                                                
				continue;                                                                                                              
			                                                                                                                         
			int* road = allRoad[firstBondAtomInd][endInd];
			if (isContain(road, node_num, secondBondAtomInd)){                                                                         
				continue;                                                                                                              
			}                                                                                                                        
			                                                                                                                         
			for (int j = 0; j < node_num; j++){                                                                                   
				if (road[j] < 0)                                                                                                       
					break;                                                                                                               
				                                                                                                                       
				if (!isFirstPart[road[j]]){                                                                                            
					firstPartAtomInds.push_back(road[j]);                                                                                      
					isFirstPart[road[j]] = true;                                                                                         
				}                                                                                                                      
			}  
			
			road = NULL; // code added by Jun at 20180915
		}
		                                                                                                                         
		vector<int> secondPartAtomInds;                                                                
		bool* isSecondPart = new bool[node_num];
		for (i = 0; i < node_num; i++) isSecondPart[i] = false;  
		for (endInd = 0; endInd < node_num; endInd++){                                                                          
			if (isSecondPart[endInd])
				continue;                                                                                                              
			                                                                                                                         
			int* road = allRoad[secondBondAtomInd][endInd];                                                                         
			if (isContain(road, node_num, firstBondAtomInd)){                                                                          
				continue;                                                                                                              
			}                                                                                                                        
			                                                                                                                         
			for (int j = 0; j < node_num; j++){                                                                                   
				if (road[j] < 0)                                                                                                       
					break;                                                                                                               
				                                                                                                                       
				if (!isSecondPart[road[j]]){                                                                                           
					secondPartAtomInds.push_back(road[j]);                                                                                     
					isSecondPart[road[j]] = true;                                                                                        
				}                                                                                                                      
			}   
			
			road = NULL; // code added by Jun at 20180915
		}
		                                                                                                                           
		// here we could consider whether the single bond is rotatable bond or not                                                 
		bool isRotatable = true;                                                                                                
		for (i = 0; i < node_num; i++){                                                                                         
			if (isFirstPart[i] && isSecondPart[i]){   
				isRotatable = false;                                                                                                   
				break;                                                                                                                 
			}                                                                                                             
		}
		
		// Code added by Jun at 20180915 START
		if (isRotatable){
			for (i = 0; i < node_num; i++){
				for (j = i+1; j < node_num; j++){
					if (i == firstBondAtomInd || j == firstBondAtomInd
						|| i == secondBondAtomInd || j == secondBondAtomInd){
						continue;
					}

					if (   (isFirstPart[i] != isFirstPart[j])
						|| (isSecondPart[j] != isSecondPart[i])){
						int* road = allRoad[i][j];
						
						if (!isContain(road, node_num, firstBondAtomInd)
							&& !isContain(road, node_num, secondBondAtomInd)){                                                                          
							isRotatable = false;  
							road = NULL;
							break;                                                                                                           
						}   
					}
				}
				if (!isRotatable){
					break;
				}
			}
		}
		// Code added by Jun at 20180915 END
		
		// B33: a partition that misses atoms (e.g. a disconnected ligand) must
		// never be rotated -- rotating it would drop those atoms downstream.
		// Only publish the split when both sides together cover every atom.
		int numCovered = 0;
		vector<bool> isCovered(node_num, false);
		for (i = 0; i < (int)firstPartAtomInds.size(); i++){
			if (!isCovered[firstPartAtomInds[i]]){
				isCovered[firstPartAtomInds[i]] = true;
				numCovered++;
			}
		}
		for (i = 0; i < (int)secondPartAtomInds.size(); i++){
			if (!isCovered[secondPartAtomInds[i]]){
				isCovered[secondPartAtomInds[i]] = true;
				numCovered++;
			}
		}

		vector<vector<int> > ans;           
		if (isRotatable && numCovered == node_num){                                                                                    
			ans.push_back(firstPartAtomInds);
			ans.push_back(secondPartAtomInds);  
		}  
		
		delete[] isFirstPart;
		isFirstPart = NULL;
		delete[] isSecondPart;
		isSecondPart = NULL;
		allRoad = NULL;  // Code added at 20180915 by Jun
		return ans;                                                                                                                
	}  
	
	/**                                                                                                                          
	 * @return the single bond (without 'H' atom) int[0] and int[1] means the two atoms, respectively, of the corresponding bond.
	 */                                                                                                                          
	vector<BondInfo> findAllSingleBond(){
		vector<BondInfo> singleBond;                                                                            
                                                                                
		vector<BondInfo> bindInfo = lig.getBondInfoes();                                                                        
		for (int i = 0; i < bindInfo.size(); i++){ // disconsider the first and last bond                                       
			BondInfo ithBI = bindInfo[i];
			
			string aAtomType = lig.getAtomType(ithBI.aInd);
			toUpperString(aAtomType); 
			string bAtomType = lig.getAtomType(ithBI.bInd);
			toUpperString(bAtomType);
			
			if (SINGLE == ithBI.bt                                                                           
					&& 0 != strcmp(aAtomType.c_str(), "H")                                                              
					&& 0 != strcmp(bAtomType.c_str(), "H")){                                                                            
				singleBond.push_back(bindInfo[i]);                                                                                    
			}                                                                                                                        
		}                                                                                                                          
		                                                                                                                           
		return singleBond;                                                                                                         
	} 
	
	LigRotamerDetector(const LigRotamerDetector&);
	LigRotamerDetector& operator=(const LigRotamerDetector&);
};

#endif
