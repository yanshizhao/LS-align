#ifndef LIGANDS_H
#define LIGANDS_H

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <fstream>
#include "Point3D.h"
#include "Str.h"
#include "DijkstraSRoad.h"
using namespace std;

enum BondType{
		// MOL2 original definition
		SINGLE,  			// 1 = single
		DOUBLE,  			// 2 = double
		TRIPLE,  			// 3 = triple
		AMIDE,  			// 4 = am (amide)
		AROMATIC, 			// 5 = ar (aromatic)
		DUMMY,  			// 6 = du (dummy) 
		UNKNOWN,  			// 7 = un (unknown) (cannot be determined from the parameter tables)
		NOTCONNECTED, 			// 8 = nc (not connected)
		
		// J. Hu defined, only use in aligned, not for load Mol2
		SINGLE_RING_THREE,	// means the single bond in four atom ring
		SINGLE_RING_FOUR,	// means the single bond in four atom ring
		SINGLE_RING_FIVE,	// means the single bond in four atom ring
		SINGLE_RING_SIX,	// means the single bond in four atom ring
		
		DOUBLE_RING_THREE,	// means the single bond in four atom ring
		DOUBLE_RING_FOUR,	// means the single bond in four atom ring
		DOUBLE_RING_FIVE,	// means the single bond in four atom ring
		DOUBLE_RING_SIX,	// means the single bond in four atom ring
		
		TRIPLE_RING_THREE,	// means the single bond in four atom ring
		TRIPLE_RING_FOUR,	// means the single bond in four atom ring
		TRIPLE_RING_FIVE,	// means the single bond in four atom ring
		TRIPLE_RING_SIX,	// means the single bond in four atom ring
		
		AMIDE_RING_THREE,	// means the single bond in four atom ring
		AMIDE_RING_FOUR,	// means the single bond in four atom ring
		AMIDE_RING_FIVE,	// means the single bond in four atom ring
		AMIDE_RING_SIX,	// means the single bond in four atom ring
		
		AROMATIC_RING_THREE,	// means the single bond in four atom ring
		AROMATIC_RING_FOUR,	// means the single bond in four atom ring
		AROMATIC_RING_FIVE,	// means the single bond in four atom ring
		AROMATIC_RING_SIX,	// means the single bond in four atom ring
		
		DUMMY_RING_THREE,	// means the single bond in four atom ring
		DUMMY_RING_FOUR,	// means the single bond in four atom ring
		DUMMY_RING_FIVE,	// means the single bond in four atom ring
		DUMMY_RING_SIX,	// means the single bond in four atom ring
		
		UNKNOWN_RING_THREE,	// means the single bond in four atom ring
		UNKNOWN_RING_FOUR,	// means the single bond in four atom ring
		UNKNOWN_RING_FIVE,	// means the single bond in four atom ring
		UNKNOWN_RING_SIX,	// means the single bond in four atom ring
};

typedef struct {
	int aInd;			// the first atom index in the original ligand
	int bInd;			// the second atom index in the original ligand
	BondType bt;					// the bond type
	double bondDisPow2;		// the distance^2 between two bind atoms
	double bondDis;
}BondInfo;

class Ligand{
private:
	string ligname;
	vector<Point3D> points;								// the coordinate of each atom
	vector< vector<BondType> > atomBondTypes;			// record the bondtypes of each atom
	
	vector<string> atomtypes;			// the atom type is Upper word
	vector<BondInfo> bondInfoes;		// the bond	information in the ligand, without Jun definition bond type
	BondType** btMtx;					// record the i - j bond type
	
	double** bondDisMtx;				// not-connect is setted as -1
public:
	Ligand(string ligname, vector<Point3D> points, vector<string> atomtypes, vector<BondInfo> bondInfoes){
		this->ligname = ligname;
		this->points = points;
		this->atomtypes = atomtypes;
		this->bondInfoes = bondInfoes;
		
		btMtx = NULL;
		bondDisMtx = NULL;
	}
	
	/************************************************	
	* Allocate new memory for points
	* no new memory for atomtypes and bondInfoes
	************************************************/
	Ligand(const Ligand& lig){
		int i = 0;
		for (i = 0; i < lig.points.size(); i++){
			Point3D pt(lig.points[i]);
			this->points.push_back(pt);
		}
		
		this->ligname = lig.ligname;
		this->atomtypes = lig.atomtypes;
		this->bondInfoes = lig.bondInfoes;
		
		if (NULL != lig.btMtx){
			btMtx = new BondType*[points.size()];
			for (int i = 0; i < points.size(); i++){
				btMtx[i] = new BondType[points.size()];
				for (int j = 0; j < points.size(); j++){
					btMtx[i][j] = lig.btMtx[i][j];
				}
			}
		}else{
			btMtx = NULL;
		}
		
		if (NULL != lig.bondDisMtx){
			bondDisMtx = new double*[points.size()];
			for (int i = 0; i < points.size(); i++){
				bondDisMtx[i] = new double[points.size()];
				for (int j = 0; j < points.size(); j++){
					bondDisMtx[i][j] = lig.bondDisMtx[i][j];
				}
			}
		}else{
			bondDisMtx = NULL;
		}
		
	}	
	
	void operator=(const Ligand& lig){
		for (int i = 0; i < lig.points.size(); i++){
			Point3D pt(lig.points[i]);
			this->points.push_back(pt);
		}
		
		this->ligname = lig.ligname;
		this->atomtypes = lig.atomtypes;
		this->bondInfoes = lig.bondInfoes;
		
		if (NULL != lig.btMtx){
			btMtx = new BondType*[points.size()];
			for (int i = 0; i < points.size(); i++){
				btMtx[i] = new BondType[points.size()];
				for (int j = 0; j < points.size(); j++){
					btMtx[i][j] = lig.btMtx[i][j];
				}
			}
		}else{
			btMtx = NULL;
		}
		
		if (NULL != lig.bondDisMtx){
			bondDisMtx = new double*[points.size()];
			for (int i = 0; i < points.size(); i++){
				bondDisMtx[i] = new double[points.size()];
				for (int j = 0; j < points.size(); j++){
					bondDisMtx[i][j] = lig.bondDisMtx[i][j];
				}
			}
		}else{
			bondDisMtx = NULL;
		}
	}
	
	Ligand(){btMtx = NULL;}
	
	~Ligand(){
		if (NULL != btMtx){
			for (int i = 0; i < points.size(); i++){
				delete[] btMtx[i];
			}
			delete[] btMtx;
			
			btMtx = NULL;
		}
		
		if (NULL != bondDisMtx){
			for (int i = 0; i < points.size(); i++){
				delete[] bondDisMtx[i];
			}
			delete[] bondDisMtx;
			
			bondDisMtx = NULL;
		}
	}
	
	BondType** getBondTypeMtx(){
		int i;
		
		if (NULL == btMtx){
			btMtx = new BondType*[points.size()];
			for (i = 0; i < points.size(); i++){
				btMtx[i] = new BondType[points.size()];
				for (int j = 0; j < points.size(); j++){
					btMtx[i][j] = NOTCONNECTED;
				}
			}
			
			for (i = 0; i < bondInfoes.size(); i++){
				int aInd = bondInfoes[i].aInd;
				int bInd = bondInfoes[i].bInd;
				BondType bt = bondInfoes[i].bt;
				
				if (aInd == bInd){
					cout << "getBondTypeMtx says " << aInd << "\t=\t" << bInd << endl;
				}
				
				btMtx[aInd][bInd] = bt;
				btMtx[bInd][aInd] = bt;
			}
		}
		
		return btMtx;
	}
	
	double** getBondDistMtx(){
		int i;
		
		if (NULL == bondDisMtx){
			bondDisMtx = new double*[points.size()];
			for (i = 0; i < points.size(); i++){
				bondDisMtx[i] = new double[points.size()];
				for (int j = 0; j < points.size(); j++){
					bondDisMtx[i][j] = -1.0;
				}
			}
			
			for (i = 0; i < bondInfoes.size(); i++){
				int aInd = bondInfoes[i].aInd;
				int bInd = bondInfoes[i].bInd;
				double dist = bondInfoes[i].bondDis;
				
				if (aInd == bInd){
					cout << "getBondTypeMtx says " << aInd << "\t=\t" << bInd << endl;
				}
				
				bondDisMtx[aInd][bInd] = dist;
				bondDisMtx[bInd][aInd] = dist;
			}
		}
		
		return bondDisMtx;
	}
	
	vector<Point3D> getPoints(){
		return points;
	}
	
	void setPoints(vector<Point3D> pts){
		vector<Point3D>().swap(this->points);
		
		this->points = pts;
	}
	
	Point3D getIthPoints(int ind){
		return points[ind];
	}
	
	vector<string> getAtomTypes(){
		return atomtypes;
	}
	
	string getAtomType(int ind){
		return atomtypes[ind];
	}
	
	vector<BondInfo> getBondInfoes(){
		return bondInfoes;
	}
	
	int size(){
		return points.size();
	}
	
	string getLigname(){
		return ligname;
	}
	
	/*********************************************************
	 * With Jun Definition bond types
	 * vector<BondType> ithAtomBt = ans[i]; 
	 * ans.size() == points.size()
	 *********************************************************/
	vector< vector<BondType> > getAtomBondTypes(){
		if (0 == atomBondTypes.size()){
			int i = 0, j = 0;
			
			vector< vector<int> > rings;
			getBondDistMtx();
			for (i = 0; i < bondInfoes.size(); i++){
				BondInfo bi = bondInfoes[i];
				DijkstraSRoad* dsr = new DijkstraSRoad(bondDisMtx, size());
				vector<int> road = dsr->getRoad(bi.aInd, bi.bInd, true);

				if (road.size() > 2 && road.size() < 7){ // judge is ring
					// it is a ring

					bool isRingExisted = false;
					for (j = 0; j < rings.size(); j++){
						if (isSame(rings[j], road)){
							isRingExisted = true;
						}
					}

					if (!isRingExisted){
						rings.push_back(road);
					}
				}

				delete dsr;
				dsr = NULL;
			}// end for_i
			
			getBondTypeMtx();
			for (i = 0; i < size(); i++){
				vector<BondType> ithAtomBts;

				// set mol2 bond type
				for (j = 0; j < size(); j++){
					if (NOTCONNECTED == btMtx[i][j])
						continue;

					ithAtomBts.push_back(btMtx[i][j]);
				}

				for (j = 0; j < rings.size(); j++){
					vector<int> ring = rings[j];

					int ind = indexOf(ring, i);
					if (-1 == ind){
						continue;
					}else{
						BondType prior_obi = btMtx[ring[ind-1<0 ? (int)ring.size()-1 : ind-1]][ring[ind]];
						BondType prior_junBi = toJunDefineBondType(prior_obi, ring.size());
						ithAtomBts.push_back(prior_junBi);

						BondType rail_obi = btMtx[ring[ind]][ring[ind+1>=ring.size() ? 0 : ind+1]];
						BondType rail_junBi = toJunDefineBondType(rail_obi, ring.size());
						ithAtomBts.push_back(rail_junBi);
					}
				}// end for_j

				atomBondTypes.push_back(ithAtomBts);
			}// end for _i
				
		}// end if 

		return atomBondTypes;
	}

	string toString(){
		int i = 0;
		
		string ans;
		ans += "ligname : " + ligname + "\n";
		
		for(i = 0; i < size(); i++){
			ans += points[i].toString() + "\t" + atomtypes[i] + "\n";
		}
		
		for(i = 0; i < bondInfoes.size(); i++){
			char tmp[1000];
			sprintf(tmp, "%d - %d    %d\n", bondInfoes[i].aInd, bondInfoes[i].bInd, bondInfoes[i].bt);
			ans += tmp;
		}
		
		return ans;
	}
	
    /*******************************************************
	* format this ligand as pdb file content
	*******************************************************/
	vector<string> format2PDB(){
		vector<string> ans;
		for (int i = 0; i < points.size(); i++){
			//HETATM 7590  PG  ATP B1800       1.475  -6.506  12.945  1.00 36.69           P  
			char buf[1000];
			sprintf(buf, "HETATM %4d      JUN 88888    %8.3f%8.3f%8.3f  1.00  0.00         %3s", 
				(i+1), points[i][0], points[i][1], points[i][2], atomtypes[i].c_str());			
			string str(buf);
			ans.push_back(str);
		}
		ans.push_back("TER");
		
		return ans;
	}
	
	vector<string> format2PDB(string ligtype){
		int i = 0;
		vector<string> ans;
		char lt[4];
		for (i = 0; i < 3; i++){
			if (i < ligtype.size()){
				lt[i] = ligtype[i];
			}else{
				lt[i] = ' ';
			}
		}
		lt[3] = '\0';
		
		for (i = 0; i < points.size(); i++){
			//HETATM 7590  PG  ATP B1800       1.475  -6.506  12.945  1.00 36.69           P  
			char buf[1000];
			sprintf(buf, "HETATM %4d      %3s 88888    %8.3f%8.3f%8.3f  1.00  0.00         %3s", 
				(i+1), lt, points[i][0], points[i][1], points[i][2], atomtypes[i].c_str());			
			string str(buf);
			ans.push_back(str);
		}
		ans.push_back("TER");
		
		return ans;
	}
	
	vector<string> format2PDB(string ligtype, int atom_start_index, int mol_no){
		int i = 0;
		vector<string> ans;
		char lt[4];
		for (i = 0; i < 3; i++){
			if (i < ligtype.size()){
				lt[i] = ligtype[i];
			}else{
				lt[i] = ' ';
			}
		}
		lt[3] = '\0';
		
		for (i = 0; i < points.size(); i++){
			//HETATM 7590  PG  ATP B1800       1.475  -6.506  12.945  1.00 36.69           P  
			char buf[1000];
			sprintf(buf, "HETATM %4d      %3s %5d    %8.3f%8.3f%8.3f  1.00  0.00         %3s", 
				(i+atom_start_index), lt, mol_no, points[i][0], points[i][1], points[i][2], atomtypes[i].c_str());			
			string str(buf);
			ans.push_back(str);
		}
		ans.push_back("TER");
		
		return ans;
	}
	
	void save8PDBFormat(string savepdb){
		ofstream fout(savepdb.c_str());
		vector<string> contents = format2PDB();		
		
		for (int i = 0; i < contents.size(); i++){
			fout << contents[i] << endl;
		}
		
		fout.close();
	}

private:
	BondType toJunDefineBondType(BondType obi, int ring_size){
		if (SINGLE == obi){
			if (3 == ring_size){
				return (SINGLE_RING_THREE);
			}else if (4 == ring_size){
				return (SINGLE_RING_FOUR);
			}else if (5 == ring_size){
				return (SINGLE_RING_FIVE);
			}else if (6 == ring_size){
				return (SINGLE_RING_SIX);
			}
		}else if (DOUBLE == obi){
			if (3 == ring_size){
				return (DOUBLE_RING_THREE);
			}else if (4 == ring_size){
				return (DOUBLE_RING_FOUR);
			}else if (5 == ring_size){
				return (DOUBLE_RING_FIVE);
			}else if (6 == ring_size){
				return (DOUBLE_RING_SIX);
			}
		}else if (TRIPLE == obi){
			if (3 == ring_size){
				return (TRIPLE_RING_THREE);
			}else if (4 == ring_size){
				return (TRIPLE_RING_FOUR);
			}else if (5 == ring_size){
				return (TRIPLE_RING_FIVE);
			}else if (6 == ring_size){
				return (TRIPLE_RING_SIX);
			}
		}else if (AMIDE == obi){
			if (3 == ring_size){
				return (AMIDE_RING_THREE);
			}else if (4 == ring_size){
				return (AMIDE_RING_FOUR);
			}else if (5 == ring_size){
				return (AMIDE_RING_FIVE);
			}else if (6 == ring_size){
				return (AMIDE_RING_SIX);
			}
		}else if (AROMATIC == obi){
			if (3 == ring_size){
				return (AROMATIC_RING_THREE);
			}else if (4 == ring_size){
				return (AROMATIC_RING_FOUR);
			}else if (5 == ring_size){
				return (AROMATIC_RING_FIVE);
			}else if (6 == ring_size){
				return (AROMATIC_RING_SIX);
			}
		}else if (DUMMY == obi){
			if (3 == ring_size){
				return (DUMMY_RING_THREE);
			}else if (4 == ring_size){
				return (DUMMY_RING_FOUR);
			}else if (5 == ring_size){
				return (DUMMY_RING_FIVE);
			}else if (6 == ring_size){
				return (DUMMY_RING_SIX);
			}
		}else if (UNKNOWN == obi){
			if (3 == ring_size){
				return (UNKNOWN_RING_THREE);
			}else if (4 == ring_size){
				return (UNKNOWN_RING_FOUR);
			}else if (5 == ring_size){
				return (UNKNOWN_RING_FIVE);
			}else if (6 == ring_size){
				return (UNKNOWN_RING_SIX);
			}
		}
		
		return NOTCONNECTED;
	}
};

/******************************************
* In this class file, we just load mol2 format
* file, no matter the file contains how many ligands
******************************************/
class Ligands{
private:
	string mol2path;
	vector<Ligand> ligs;
	
	bool isLoad_H_atom;
	
public:
	Ligands(string mol2path, bool isLoad_H_atom){
		this->mol2path = mol2path;
		this->isLoad_H_atom = isLoad_H_atom;
		
		load();
	}
	
	~Ligands(){
		vector<Ligand>().swap(ligs);
	}	
	
	vector<Ligand> getLigands(){
		return ligs;
	}
	
	Ligand operator[](int ind){
		return ligs[ind];
	}
	
	int size(){
		return ligs.size();
	}
private:
	
/*****************************************
* we load all the ligands in mol2path
	*****************************************/
	void load(){
		vector<vector<string> > mol2StrInfos = loadAllLigStrInfoMol2();
		
		for (int i = 0; i < mol2StrInfos.size(); i++){
			vector<string> mol2Info = mol2StrInfos[i];
			Ligand lig = parseOneLigMol2(mol2Info);
			
			ligs.push_back(lig);
		}
		
	}	
	
	Ligand parseOneLigMol2(vector<string> oneLigMol2Infos){
		string ligname;
		vector<Point3D> points;
		vector<string> atomtypes;
		vector<BondInfo> bondInfoes;
		
		int l = 0;
		while (l < oneLigMol2Infos.size() && oneLigMol2Infos[l].compare(0, 13, "@<TRIPOS>ATOM") != 0){
			if (oneLigMol2Infos[l].compare(0, 17, "@<TRIPOS>MOLECULE") == 0){
				ligname = oneLigMol2Infos[++l];
			}
			
			l++;
		}
		
		l++;
		while (l < oneLigMol2Infos.size()){
			if (oneLigMol2Infos[l].compare(0, 9, "@<TRIPOS>") == 0){
				break;
			}
			
			char cstr[50];
			double x = 0.0;
			strcpy(cstr, (oneLigMol2Infos[l].substr(17, 9)).c_str());
			sscanf(cstr, "%lf", &x);	
			
			double y = 0.0;
			strcpy(cstr, (oneLigMol2Infos[l].substr(27, 9)).c_str());
			sscanf(cstr, "%lf", &y);
			
			double z = 0.0;
			strcpy(cstr, (oneLigMol2Infos[l].substr(37, 9)).c_str());
			sscanf(cstr, "%lf", &z);
			
			Point3D pos(x, y, z);
			points.push_back(pos);
			

			//string atomtype = oneLigMol2Infos[l].substr(47, 1);
			//if (atomtype.compare(0, 1, " ") == 0)
			//	atomtype = oneLigMol2Infos[l].substr(48, 1);
			string atomtype = oneLigMol2Infos[l].substr(47, 2);
			if (atomtype.compare(0, 1, " ") == 0)
				atomtype = oneLigMol2Infos[l].substr(48, 2);
			atomtype = eraseAll(atomtype, " .1234567890~!@#$%^&*()_+", 25);

			toUpperString(atomtype);
			atomtypes.push_back(atomtype);
			
			l++;
		}
		
		int ll = l;
		while (ll < oneLigMol2Infos.size() && oneLigMol2Infos[ll].compare(0, 13, "@<TRIPOS>BOND") != 0){
			ll++;
		}
		
		ll++;
		while (ll < oneLigMol2Infos.size()){
			if (oneLigMol2Infos[ll].compare(0, 9, "@<TRIPOS>") == 0){
				break;
			}
			
			vector<string> llc = split(oneLigMol2Infos[ll], ' ', '\t');
			if (llc.size() < 4){
				ll++;
				continue;
			}
			
			int atom1Ind = 0;
			sscanf(llc[1].c_str(), "%d", &atom1Ind);
			atom1Ind--;
			
			int atom2Ind = 0;
			sscanf(llc[2].c_str(), "%d", &atom2Ind);
			atom2Ind--;
			
			BondType bondtype;	
			if (llc[3].size() == 1 && llc[3].compare(0, 1, "1")==0){
				bondtype = SINGLE;
			}else if (llc[3].size() == 1 && llc[3].compare(0, 1, "2")==0){
				bondtype = DOUBLE;
			}else if (llc[3].size() == 1 && llc[3].compare(0, 1, "3")==0){
				bondtype = TRIPLE;
			}else if (llc[3].size() == 2 && llc[3].compare(0, 2, "am")==0){
				bondtype = AMIDE;		
			}else if (llc[3].size() == 2 && llc[3].compare(0, 2, "ar")==0){
				bondtype = AROMATIC;	
			}else if (llc[3].size() == 2 && llc[3].compare(0, 2, "du")==0){
				bondtype = DUMMY;		
			}else if (llc[3].size() == 2 && llc[3].compare(0, 2, "un")==0){
				bondtype = UNKNOWN;		
			}else /*if (llc[3].size() == 2 && llc[3].compare(0, 2, "nc")==0)*/{
				bondtype = NOTCONNECTED;		
			}
			
			BondInfo bondInfo;
			bondInfo.aInd = atom1Ind;
			bondInfo.bInd = atom2Ind;
			bondInfo.bt = bondtype;
			bondInfo.bondDisPow2 = points[atom1Ind].distPow2(points[atom2Ind]);
			bondInfo.bondDis = sqrt(bondInfo.bondDisPow2);
			
			bondInfoes.push_back(bondInfo);
			
			ll++;
		}
		
		if (!isLoad_H_atom){			
			vector<string>::iterator at_it;
			vector<Point3D>::iterator pt_it;
			int *oldMapNewInd = new int[atomtypes.size()];
			int oldInd = 0, NewInd = 0;
			
			for (at_it = atomtypes.begin(), pt_it = points.begin(); at_it != atomtypes.end();){
				if (0 == at_it->compare(0, at_it->size(), "H")
					|| 0 == at_it->compare(0, at_it->size(), "h")){
					at_it = atomtypes.erase(at_it);
					pt_it = points.erase(pt_it);
					
					oldMapNewInd[oldInd++] = -1;
				}else{
					at_it++;
					pt_it++;
					
					oldMapNewInd[oldInd++] = NewInd++;
				}
			}
			
			
			vector<BondInfo>::iterator bi_it;
			for (bi_it = bondInfoes.begin(); bi_it != bondInfoes.end();){
				int aInd = bi_it->aInd;
				int bInd = bi_it->bInd;
				
				if (-1 == oldMapNewInd[aInd]
					|| -1 == oldMapNewInd[bInd]){					
					bi_it = bondInfoes.erase(bi_it);
				}else{ 
					bi_it->aInd = oldMapNewInd[aInd];
					bi_it->bInd = oldMapNewInd[bInd];
					
					bi_it++;
				}
			}
			
			
			delete[] oldMapNewInd;
			oldMapNewInd = NULL;
		}
		
		Ligand ans(ligname, points, atomtypes, bondInfoes);
		return ans;
	}
	
	/***********************************************
	* load all the ligand information in the mol2 file
	* using the string container
	***********************************************/
	vector<vector<string> > loadAllLigStrInfoMol2(){
		vector<vector<string> > ans;
		
		bool isFirst = true;
		vector<string> mol2Sb;
		ifstream fin(mol2path.c_str());
		string line;
		getline(fin, line);
		while (fin.good()){
			line = eraseAll(line, '\r');	// Since the differnece between windows (\r\n) and Linux (\n), so we should do this step.
			if (line.compare(0, 17, "@<TRIPOS>MOLECULE") == 0){
				if (isFirst){
					isFirst = false;
				}else{
					ans.push_back(mol2Sb);
					mol2Sb.erase(mol2Sb.begin(), mol2Sb.end());
				}
			}
			
			mol2Sb.push_back(line);
			getline(fin, line);
			if (!fin.good()){
				mol2Sb.push_back(line);
			}
		}
		fin.close();
		
		if (!isFirst)
			ans.push_back(mol2Sb);
		return ans;
	}
};

#endif
