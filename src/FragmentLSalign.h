#ifndef FRAGMENT_LSALIGN_H
#define FRAGMENT_LSALIGN_H

#include "Ligands.h"
#include "ShortestRoad.h"
#include "RigidLSalign.h"
#include "basic_fun.h"
#include "NW_DP.h"
using namespace std;

typedef struct{	
	vector<int> atominds;		// the atom index vector in the original ligand; index start from 0
	vector<Point3D> atompts; 	// the atom point coordinate vector
	vector<string> atomtypes; 	// the atom type vector 
	vector<BondInfo> bondInfoes; 	// the bond information of the atompts (and atomtypes)
	
	int link_ind;	// the interface atom index in the original ligand
					// if the value is -1, it means the last fragment in the .
					// so the fragment is the whole ligand information
}Fragment;

class FragmentLSalign{
private:
	Ligand query;			// [in] the input query ligand
	Ligand templ;			// [in] the input templ ligand
	double d0;			// [in] the input d0 for LS-score calculation
	
	Ligand* p_query_sup;			// [out] the changed query ligand after superposition	
	vector<AlignedPair> finalAlignedPair;	// [out] the final aligned information
 	
	vector<Fragment>	query_fragments;		// record the original fragment information
	vector<Fragment>	templ_fragments;		// record the original fragment information

	double alisco;			// the ali score. LS-score = alisco / (query_fragments.size()>templ_fragments.size() ? query_fragments.size() : templ_fragments.size())

public:
	FragmentLSalign(Ligand orig_query, Ligand templ, double d0){
		RigidLSalign *p_origi_rigid_LSalign = new RigidLSalign(orig_query, templ);
		this->query = p_origi_rigid_LSalign->getRottedQueryLigand();
			
		this->templ = templ;
		this->d0 = d0;
		p_query_sup = NULL;

		align();

		delete p_origi_rigid_LSalign;
		p_origi_rigid_LSalign = NULL;
	}

	~FragmentLSalign(){
		if (NULL != p_query_sup){
			delete p_query_sup;
			p_query_sup = NULL;
		}
	}

	double getLSscore(){
		return alisco / (query.size()>templ.size() ? query.size() : templ.size());
	}

	double getLSscore8QuerySize(){
		return alisco / query.size();
	}

	double getLSscore8TemplSize(){
		return alisco / templ.size();
	}

	vector<Fragment> getQueryFragmets(){
		return query_fragments;
	}

	vector<Fragment> getTemplFragmets(){
		return templ_fragments;
	}

	vector<AlignedPair> getFinalAliPair(){
		return finalAlignedPair;
	}
	
	double getJaccardRatio(){
		double alinum = finalAlignedPair.size();
		return (1.0*alinum) / ((int)query.size() + (int)templ.size() - alinum);
	}

	Ligand* getQuerySupLigand(){
		return p_query_sup;
	}

	vector<string> formatAliInfo(){
		vector<string> ans;

		string queryAtomInd;
		string queryAtomType;
		string aliTag;	// distance < d0/2.0 using :, d0/2.0 <= distance < d0 using .
		string aliDis;
		string templAtomType;
		string templAtomInd;
		
		for (int i = 0; i < finalAlignedPair.size(); i++){
			AlignedPair ithPair = finalAlignedPair[i];

			int qAI = ithPair.rowInd + 1;		// index start from 1
			string qAT = query.getAtomType(ithPair.rowInd);
			double dis = p_query_sup->getIthPoints(ithPair.rowInd).dist(templ.getIthPoints(ithPair.colInd));

			char tag[2];
			tag[0] = '.';
			if (dis < d0/2.0)
				tag[0] = ':';
			tag[1] = '\0';
			
			string tAT = templ.getAtomType(ithPair.colInd);
			int tAI = ithPair.colInd;

			char buf[20];
			sprintf(buf, "%4d ", qAI);
			queryAtomInd += buf;

			sprintf(buf, "%4s ", qAT.c_str());
			queryAtomType += buf;

			sprintf(buf, "%4s ", tag);
			aliTag += buf;
			
			sprintf(buf, "%4.2f ", dis);
			aliDis += buf;

			sprintf(buf, "%4s ", tAT.c_str());
			templAtomType += buf;

			sprintf(buf, "%4d ", tAI);
			templAtomInd += buf;

			if (queryAtomInd.size() >= 120){
				ans.push_back("Query Atom Index:\t" + queryAtomInd);
				queryAtomInd.erase(queryAtomInd.begin(), queryAtomInd.end());
				
				ans.push_back("Query Atom Type: \t" + queryAtomType);
				queryAtomType.erase(queryAtomType.begin(), queryAtomType.end());

				ans.push_back("Aligned Tag:     \t" + aliTag);
				aliTag.erase(aliTag.begin(), aliTag.end());

				ans.push_back("Aligned Distance:\t" + aliDis);
				aliDis.erase(aliDis.begin(), aliDis.end());

				ans.push_back("Templ Atom Type: \t" + templAtomType);
				templAtomType.erase(templAtomType.begin(), templAtomType.end());

				ans.push_back("Templ Atom Index:\t" + templAtomInd);
				templAtomInd.erase(templAtomInd.begin(), templAtomInd.end());

				ans.push_back("    ");
			}
		}

		if (queryAtomInd.size() > 0){
			ans.push_back("Query Atom Index:\t" + queryAtomInd);
			queryAtomInd.erase(queryAtomInd.begin(), queryAtomInd.end());
			
			ans.push_back("Query Atom Type: \t" + queryAtomType);
			queryAtomType.erase(queryAtomType.begin(), queryAtomType.end());
			
			ans.push_back("Aligned Tag:     \t" + aliTag);
			aliTag.erase(aliTag.begin(), aliTag.end());
			
			ans.push_back("Aligned Distance:\t" + aliDis);
			aliDis.erase(aliDis.begin(), aliDis.end());

			ans.push_back("Templ Atom Type: \t" + templAtomType);
			templAtomType.erase(templAtomType.begin(), templAtomType.end());

			ans.push_back("Templ Atom Index:\t" + templAtomInd);
			templAtomInd.erase(templAtomInd.begin(), templAtomInd.end());

			ans.push_back("    ");
		}
		
		return ans;
	}

private:
	void align(){
		int i = 0, j = 0, k = 0;		

		query_fragments = divideLigandToFragment(query);
		templ_fragments = divideLigandToFragment(templ);
		
		RigidLSalign** *rigidLSalignMtx = new RigidLSalign**[query_fragments.size()];
		double** rigidLSscoreMtx = new double*[query_fragments.size()];
		for (i = 0; i < query_fragments.size(); i++){
			rigidLSalignMtx[i] = new RigidLSalign*[templ_fragments.size()];
			rigidLSscoreMtx[i] = new double[templ_fragments.size()];

			for (j = 0; j < templ_fragments.size(); j++){
				Ligand* p_query_frag = new Ligand("QUERYfrag", query_fragments[i].atompts, query_fragments[i].atomtypes, query_fragments[i].bondInfoes);
				Ligand* p_templ_frag = new Ligand("TEMPLfrag", templ_fragments[j].atompts, templ_fragments[j].atomtypes, templ_fragments[j].bondInfoes);

				rigidLSalignMtx[i][j] = new RigidLSalign(*p_query_frag, *p_templ_frag, d0, 3, 1);
				rigidLSscoreMtx[i][j] = rigidLSalignMtx[i][j]->getLSscore();

				// release memory
                delete p_query_frag;
				p_query_frag = NULL;
                delete p_templ_frag;
				p_templ_frag = NULL;
			}
			
		}

		// using dynamic programming to search a good fragment alignment
		// 	double** scores, int row_num, int col_num, double gap_cost
		int* ali = runNWDP(rigidLSscoreMtx, query_fragments.size(), templ_fragments.size(), -1.0);
		alisco = 0.0;
		
		vector<Point3D> rottedQueryPts;
		vector<Fragment> unRottedFragment;
		double** lastRotMtx = NULL;
		double* lastTransVec = NULL;
		for (i = 0; i < query_fragments.size(); i++){
			unRottedFragment.push_back(query_fragments[i]);

			if (-1 != ali[i]){
					// map the new ali to old index 
				int row = i;
				int col = ali[i];
				alisco += rigidLSscoreMtx[row][col] 
					* (query_fragments[row].atominds.size()>templ_fragments[col].atominds.size() 
						? query_fragments[row].atominds.size() : templ_fragments[col].atominds.size());
				RigidLSalign* p_rigid_lsalign = rigidLSalignMtx[row][col];

				vector<AlignedPair> alipairVecOnNewIndex = p_rigid_lsalign->getFinalAlignedPair();

				for (j = 0; j < alipairVecOnNewIndex.size(); j++){
					int ali_row_on_new_index = alipairVecOnNewIndex[j].rowInd;
					int ali_col_on_new_index = alipairVecOnNewIndex[j].colInd;

					int ali_row_on_old_index = query_fragments[row].atominds[ali_row_on_new_index];
					int ali_col_on_old_index = templ_fragments[col].atominds[ali_col_on_new_index];

					AlignedPair alipair;
					alipair.rowInd = ali_row_on_old_index;
					alipair.colInd = ali_col_on_old_index;

					finalAlignedPair.push_back(alipair);
				}// end for_j

				if (NULL != lastRotMtx){
					for (j = 0; j < 3; j++){
						delete[] lastRotMtx[j];
					}
					delete[] lastRotMtx;
					lastRotMtx = NULL;
					delete[] lastTransVec;
					lastTransVec = NULL;
				}
				lastRotMtx = p_rigid_lsalign->getRotMtx();
				lastTransVec = p_rigid_lsalign->getTransVec();
				
				for (j = 0; j < unRottedFragment.size(); j++){
					// TODO using lastROTMTX to rot unRottedFragment
					for (k = 0; k < unRottedFragment[j].atompts.size(); k++){
						double nx = lastTransVec[0] + lastRotMtx[0][0]*unRottedFragment[j].atompts[k].getX() 
													+ lastRotMtx[0][1]*unRottedFragment[j].atompts[k].getY()
													+ lastRotMtx[0][2]*unRottedFragment[j].atompts[k].getZ();

						double ny = lastTransVec[1] + lastRotMtx[1][0]*unRottedFragment[j].atompts[k].getX() 
													+ lastRotMtx[1][1]*unRottedFragment[j].atompts[k].getY()
													+ lastRotMtx[1][2]*unRottedFragment[j].atompts[k].getZ();

						double nz = lastTransVec[2] + lastRotMtx[2][0]*unRottedFragment[j].atompts[k].getX() 
													+ lastRotMtx[2][1]*unRottedFragment[j].atompts[k].getY()
													+ lastRotMtx[2][2]*unRottedFragment[j].atompts[k].getZ();

						Point3D new_pt(nx, ny, nz);
						rottedQueryPts.push_back(new_pt);
					}
				}
				unRottedFragment.erase(unRottedFragment.begin(), unRottedFragment.end());
			
			
			}//end if
		}// end for_i

		// process the end unaligned segment
		for (i = 0; i < unRottedFragment.size(); i++){
			double nx = lastTransVec[0] + lastRotMtx[0][0]*unRottedFragment[j].atompts[k].getX() 
										+ lastRotMtx[0][1]*unRottedFragment[j].atompts[k].getY()
										+ lastRotMtx[0][2]*unRottedFragment[j].atompts[k].getZ();
			double ny = lastTransVec[1] + lastRotMtx[1][0]*unRottedFragment[j].atompts[k].getX() 
										+ lastRotMtx[1][1]*unRottedFragment[j].atompts[k].getY()
										+ lastRotMtx[1][2]*unRottedFragment[j].atompts[k].getZ();

			double nz = lastTransVec[2] + lastRotMtx[2][0]*unRottedFragment[j].atompts[k].getX() 
										+ lastRotMtx[2][1]*unRottedFragment[j].atompts[k].getY()
										+ lastRotMtx[2][2]*unRottedFragment[j].atompts[k].getZ();
			Point3D new_pt(nx, ny, nz);
			rottedQueryPts.push_back(new_pt);
		}

		p_query_sup = new Ligand(query);
		p_query_sup->setPoints(rottedQueryPts);

		// release memory
		delete[] ali;
		ali = NULL;
		
		if (NULL != lastRotMtx){
			for (i = 0; i < 3; i++){
				delete[] lastRotMtx[i];
			}
			delete[] lastRotMtx;
			delete[] lastTransVec;
			lastRotMtx = NULL;
			lastTransVec = NULL;
		}

		for (i = 0; i < query_fragments.size(); i++){
			for (j = 0; j < templ_fragments.size(); j++){
				delete rigidLSalignMtx[i][j];
			}
			delete[] rigidLSalignMtx[i];
			delete[] rigidLSscoreMtx[i];
		}
		delete[] rigidLSalignMtx;
		rigidLSalignMtx = NULL;
		delete[] rigidLSscoreMtx;
		rigidLSscoreMtx = NULL;
	}	
	
	vector<Fragment> divideLigandToFragment(Ligand lig){
		int i = 0, j = 0;		
		
		vector<Fragment> ans;

		ShortestRoad* p_sroad = shortestRoadBetweenTwoAtoms(lig);
		vector<int>  lig_main_backbone = p_sroad->getLongestShortestRoad8Floyd();
		
		bool* isSelected = new bool[lig.size()];
		for (i = 0; i < lig.size(); i++) // initial
			isSelected[i] = false;		

		for (i = 1; i < lig_main_backbone.size(); i++){
			int aInd = lig_main_backbone[i-1];
			int bInd = lig_main_backbone[i];
			bool isSI = judgeIsSingleBondInTwoAtoms(lig, aInd, bInd);
			if (isSI){
				vector<vector<int> > leftRight = 
					leftRightPartOnBond(p_sroad, lig, aInd, bInd);
				if (leftRight.size() == 2){
					// here, we divide the fragment from left to right
					vector<int> left = leftRight[0];
					vector<int> right = leftRight[1];
					
					int selectable_num_on_left = 0;
					for (j = 0; j < left.size(); j++){
						if (! isSelected[left[j]])
							selectable_num_on_left++;
					}
					
					if (selectable_num_on_left < 3 || right.size() < 3){
						continue;
					}
					
					// generate a fragment for ligand
					Fragment frag;
					frag.link_ind = aInd;

					for (j = 0; j < left.size(); j++){
						if (! isSelected[left[j]]){
							frag.atominds.push_back(left[j]);
							frag.atompts.push_back(lig.getIthPoints(left[j]));
							frag.atomtypes.push_back(lig.getAtomType(left[j]));

							isSelected[left[j]] = true;	
						}
					}

					vector<BondInfo> allBondInfo  = lig.getBondInfoes();
					for (j = 0; j < allBondInfo.size(); j++){
						BondInfo bi = allBondInfo[j];
						int AnewInd = indexOf(frag.atominds, bi.aInd);
						if (-1 == AnewInd) 
							continue;
						int BnewInd = indexOf(frag.atominds, bi.bInd);
						if (-1 == BnewInd) 
							continue;

						BondInfo newBI;
						newBI.aInd = AnewInd;			
						newBI.bInd = BnewInd;
						newBI.bt = bi.bt;
						newBI.bondDisPow2 = bi.bondDisPow2;
						newBI.bondDis = bi.bondDis;
						
						frag.bondInfoes.push_back(newBI);
					}			

			
					
					ans.push_back(frag);
				}
			}
		}

		if (ans.size() == 0){
			// may the ligand does not contain a rotatable single bond
			// so we just return itself
			Fragment frag;
			frag.link_ind = -1;
			
			vector<Point3D> atompts = lig.getPoints();
			vector<string> atomtypes = lig.getAtomTypes();
			for (i = 0; i < atompts.size(); i++){
				frag.atompts.push_back(atompts[i]);
				frag.atomtypes.push_back(atomtypes[i]);
				frag.atominds.push_back(i);
			}			
			
			vector<BondInfo> atombis = lig.getBondInfoes();
			for (i = 0; i < atombis.size(); i++){
				frag.bondInfoes.push_back(atombis[i]);
			}
			
			ans.push_back(frag);
		}else{
			Fragment frag;
			frag.link_ind = -1;

			// generate the last fragment
			for (i = 0; i < lig.size(); i++){
				if (! isSelected[i]){
					frag.atominds.push_back(i);
					frag.atompts.push_back(lig.getIthPoints(i));
					frag.atomtypes.push_back(lig.getAtomType(i));	
					
					isSelected[i] = true;
				}
			}

			vector<BondInfo> allBondInfo  = lig.getBondInfoes();
			for (i = 0; i < allBondInfo.size(); i++){
				BondInfo bi = allBondInfo[i];
				int AnewInd = indexOf(frag.atominds, bi.aInd);
				if (-1 == AnewInd) 
					continue;
				int BnewInd = indexOf(frag.atominds, bi.bInd);
				if (-1 == BnewInd) 
					continue;

				BondInfo newBI;
				newBI.aInd = AnewInd;			
				newBI.bInd = BnewInd;
				newBI.bt = bi.bt;
				newBI.bondDisPow2 = bi.bondDisPow2;
				newBI.bondDis = bi.bondDis;
				
				frag.bondInfoes.push_back(newBI);
			}			
			
			ans.push_back(frag);
		}

		delete[] isSelected;
		isSelected = NULL;
		delete p_sroad;
		p_sroad = NULL;

		return ans;
	}

	ShortestRoad* shortestRoadBetweenTwoAtoms(Ligand lig){
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
		return p_sroad;
	}
	
	vector<vector<int> > leftRightPartOnBond(ShortestRoad* p_sroad, Ligand lig, int firstBondAtomInd, int secondBondAtomInd){ 
		int i = 0, endInd = 0;
		int node_num = lig.size();		                                                  
		int*** allRoad = p_sroad->getAllRoad8Floyd();
		                                                                                                                           
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
		}
		                                                                                                                           
		// here we could consider whether the single bond is rotatable bond or not                                                 
		bool isRotatable = true;                                                                                                
		for (i = 0; i < node_num; i++){                                                                                         
			if (isFirstPart[i] && isSecondPart[i]){                                                                                  
				isRotatable = false;                                                                                                   
				break;                                                                                                                 
			}                                                                                                                        
		} 
		
		vector<vector<int> > ans;           
		if (isRotatable){                                                                                    
			ans.push_back(firstPartAtomInds);
			ans.push_back(secondPartAtomInds);  
		}  
		
		delete[] isFirstPart;
		isFirstPart = NULL;
		delete[] isSecondPart;
		isSecondPart = NULL;
		return ans;                                                                                                                
	}  
	
	/***********************************************                                                                                                                          
	 * Do not consider H atom.
	 ***********************************************/ 
	bool judgeIsSingleBondInTwoAtoms(Ligand lig, int aInd, int bInd){
		BondType** btMtx = lig.getBondTypeMtx();	// not new memory, you could not release it.

		string aAtomType = lig.getAtomType(aInd);
		toUpperString(aAtomType); 
		string bAtomType = lig.getAtomType(bInd);
		toUpperString(bAtomType);   
		
		if (SINGLE == btMtx[aInd][bInd]                                                                           
				&& 0 != strcmp(aAtomType.c_str(), "H")                                                              
				&& 0 != strcmp(bAtomType.c_str(), "H")){                                                                            
			return true;                                                                                    
		}

		return false;
	}	

};

#endif
