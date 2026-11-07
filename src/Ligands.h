#ifndef LIGANDS_H
#define LIGANDS_H

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <fstream>
#include <sstream>
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
	vector<string> atomNames;			// original atom names, same index as points
	vector<BondInfo> bondInfoes;		// the bond	information in the ligand, without Jun definition bond type
	BondType** btMtx;					// record the i - j bond type
	
	double** bondDisMtx;				// not-connect is setted as -1
public:
	Ligand(string ligname, vector<Point3D> points, vector<string> atomtypes, vector<BondInfo> bondInfoes,
		vector<string> atomNames = vector<string>()){
		this->ligname = ligname;
		this->points = points;
		this->atomtypes = atomtypes;
		this->atomNames = atomNames;
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
		this->atomNames = lig.atomNames;
		this->bondInfoes = lig.bondInfoes;
		
		if (NULL != lig.btMtx){
			btMtx = new BondType*[lig.points.size()];
			for (int i = 0; i < lig.points.size(); i++){
				btMtx[i] = new BondType[lig.points.size()];
				for (int j = 0; j < lig.points.size(); j++){
					btMtx[i][j] = lig.btMtx[i][j];
				}
			}
		}else{
			btMtx = NULL;
		}
		
		if (NULL != lig.bondDisMtx){
			bondDisMtx = new double*[lig.points.size()];
			for (int i = 0; i < lig.points.size(); i++){
				bondDisMtx[i] = new double[lig.points.size()];
				for (int j = 0; j < lig.points.size(); j++){
					bondDisMtx[i][j] = lig.bondDisMtx[i][j];
				}
			}
		}else{
			bondDisMtx = NULL;
		}
		
	}	
	
	Ligand& operator=(const Ligand& lig){
		if (this == &lig) return *this;

		resetCaches();
		vector<Point3D>().swap(points);

		for (int i = 0; i < (int)lig.points.size(); i++){
			points.push_back(lig.points[i]);
		}
		ligname    = lig.ligname;
		atomtypes  = lig.atomtypes;
		atomNames  = lig.atomNames;
		bondInfoes = lig.bondInfoes;

		int n = (int)lig.points.size();
		if (NULL != lig.btMtx){
			btMtx = new BondType*[n];
			for (int i = 0; i < n; i++){
				btMtx[i] = new BondType[n];
				for (int j = 0; j < n; j++) btMtx[i][j] = lig.btMtx[i][j];
			}
		} else btMtx = NULL;

		if (NULL != lig.bondDisMtx){
			bondDisMtx = new double*[n];
			for (int i = 0; i < n; i++){
				bondDisMtx[i] = new double[n];
				for (int j = 0; j < n; j++) bondDisMtx[i][j] = lig.bondDisMtx[i][j];
			}
		} else bondDisMtx = NULL;

		return *this;
	}
	
	Ligand(){ btMtx = NULL; bondDisMtx = NULL; }
	
	~Ligand(){ resetCaches(); }
	
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
					cerr << "getBondTypeMtx says " << aInd << "\t=\t" << bInd << endl;
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
					cerr << "getBondTypeMtx says " << aInd << "\t=\t" << bInd << endl;
				}
				
				bondDisMtx[aInd][bInd] = dist;
				bondDisMtx[bInd][aInd] = dist;
			}
		}
		
		return bondDisMtx;
	}
	
	const vector<Point3D>& getPoints() const{
		return points;
	}
	
	void setPoints(vector<Point3D> pts){
		vector<Point3D>().swap(this->points);
		
		this->points = pts;
	}
	
	Point3D getIthPoints(int ind) const{
		return points[ind];
	}
	
	const vector<string>& getAtomTypes() const{
		return atomtypes;
	}
	
	string getAtomType(int ind){
		return atomtypes[ind];
	}
	
	string getAtomName(int ind) const{
		if (ind < 0 || ind >= (int)atomNames.size()){
			return "";
		}
		return atomNames[ind];
	}
	
	const vector<BondInfo>& getBondInfoes() const{
		return bondInfoes;
	}
	
	int size() const{
		return points.size();
	}
	
	string getLigname() const{
		return ligname;
	}
	
	/*********************************************************
	 * With Jun Definition bond types
	 * vector<BondType> ithAtomBt = ans[i]; 
	 * ans.size() == points.size()
	 *********************************************************/
	const vector< vector<BondType> >& getAtomBondTypes(){
		if ((int)atomBondTypes.size() != (int)points.size()){
			atomBondTypes.clear();
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
			atomBondTypes.assign(points.size(), vector<BondType>());
			for (i = 0; i < size(); i++){
				vector<BondType>& ithAtomBts = atomBondTypes[i];

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
	vector<string> format2PDB() const{
		return format2PDB("JUN", 1, 88888);
	}
	
	vector<string> format2PDB(string ligtype) const{
		return format2PDB(ligtype, 1, 88888);
	}
	
	vector<string> format2PDB(string ligtype, int atom_start_index, int mol_no) const{
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
			string elem = (atomtypes[i].size() > 2) ? atomtypes[i].substr(0, 2) : atomtypes[i];
			string nameField = fitPdbAtomName(getAtomName(i), elem);
			char buf[256];
			sprintf(buf, "HETATM %4d %s %3s %5d    %8.3f%8.3f%8.3f  1.00  0.00",
				(i+atom_start_index), nameField.c_str(), lt, mol_no,
				points[i][0], points[i][1], points[i][2]);
			string line(buf);
			line.resize(76, ' ');
			if (2 == elem.size()){
				line += elem;
			}else if (1 == elem.size()){
				line += (" " + elem);
			}else{
				line += "  ";
			}
			ans.push_back(line);
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
	void resetCaches(){
		if (NULL != btMtx){
			for (int i = 0; i < (int)points.size(); i++) delete[] btMtx[i];
			delete[] btMtx;
			btMtx = NULL;
		}
		if (NULL != bondDisMtx){
			for (int i = 0; i < (int)points.size(); i++) delete[] bondDisMtx[i];
			delete[] bondDisMtx;
			bondDisMtx = NULL;
		}
		vector< vector<BondType> >().swap(atomBondTypes);
	}

	/*************************************************************
	* PDB atom-name field (columns 13-16).  Four-character names,
	* two-letter elements and names starting with a digit are
	* left-justified; other names are right-justified.  Oversized
	* names are truncated to four characters with a warning.
	*************************************************************/
	static string fitPdbAtomName(const string& rawName, const string& elem){
		string name = rawName.empty() ? elem : rawName;
		if (name.size() > 4){
			cerr << "LSalign: WARNING: atom name \"" << name
				<< "\" exceeds 4 chars, truncated for PDB output." << endl;
			name = name.substr(0, 4);
		}
		bool twoLetterElem = (2 == (int)elem.size());
		bool startsDigit = (!name.empty() && name[0] >= '0' && name[0] <= '9');
		bool leftJustify = (4 == (int)name.size()) || twoLetterElem || startsDigit;
		char buf[8];
		if (leftJustify){
			sprintf(buf, "%-4s", name.c_str());
		}else{
			sprintf(buf, " %-3s", name.c_str());
		}
		return string(buf);
	}
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

/*************************************************************
* LigandParser -- file-format parser interface.  Only mol2 is
* implemented; PDB/SDF parsers are reserved for a later ticket.
*************************************************************/
class LigandParser{
public:
	virtual ~LigandParser(){}
	virtual vector<Ligand> parse(const string& path, bool isLoadH) = 0;
	virtual const char* name() const = 0;
};

class Mol2Parser : public LigandParser{
public:
	const char* name() const { return "mol2"; }

	vector<Ligand> parse(const string& path, bool isLoadH){
		filePath = path;
		isLoad_H_atom = isLoadH;

		vector<Ligand> ans;
		vector<int> blockStartLines;
		vector<vector<string> > mol2StrInfos = loadAllLigStrInfoMol2(blockStartLines);
		for (int i = 0; i < (int)mol2StrInfos.size(); i++){
			ans.push_back(parseOneLigMol2(mol2StrInfos[i], blockStartLines[i]));
		}
		return ans;
	}
private:
	string filePath;
	bool isLoad_H_atom;

	static string trimWhitespace(const string& str){
		size_t begin = str.find_first_not_of(" \t");
		if (string::npos == begin){
			return "";
		}
		size_t end = str.find_last_not_of(" \t");
		return str.substr(begin, end - begin + 1);
	}

	static bool isFiniteCoord(double v){
		return v == v && v <= 1e308 && v >= -1e308;    // rejects NaN and +-Inf
	}

	static bool isPureNumber(const string& str){
		string trimmed = trimWhitespace(str);
		if (trimmed.empty()){
			return false;
		}
		char* endp = NULL;
		strtod(trimmed.c_str(), &endp);
		return (NULL != endp) && ('\0' == *endp);
	}

	static void splitWhitespace(const string& line, vector<string>& tokens){
		tokens.clear();
		istringstream iss(line);
		string token;
		while (iss >> token){
			tokens.push_back(token);
		}
	}

	static string reduceToElement(const string& atomType){
		string ans = eraseAll(atomType, '\t');
		ans = eraseAll(ans, " .1234567890~!@#$%^&*()_+", 25);
		toUpperString(ans);
		return ans;
	}

	static string elementFromAtomName(const string& atomName){
		string elem;
		for (int i = 0; i < (int)atomName.size(); i++){
			char ch = atomName[i];
			bool isAlpha = (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z');
			if (isAlpha){
				elem += ch;
				if (2 == (int)elem.size()){
					break;
				}
			}else if (!elem.empty()){
				break;
			}
		}
		toUpperString(elem);
		return elem;
	}

	/*************************************************************
	* Fixed columns first (all shipped data uses that layout); when
	* the slices are not pure numbers fall back to whitespace tokens.
	*************************************************************/
	bool parseAtomRecord(const string& record, string& atomName, string& atomType,
		double& x, double& y, double& z) const {
		atomName = "";
		atomType = "";

		if (record.size() >= 48){
			string sx = record.substr(17, 9);
			string sy = record.substr(27, 9);
			string sz = record.substr(37, 9);
			if (isPureNumber(sx) && isPureNumber(sy) && isPureNumber(sz)){
				x = strtod(trimWhitespace(sx).c_str(), NULL);
				y = strtod(trimWhitespace(sy).c_str(), NULL);
				z = strtod(trimWhitespace(sz).c_str(), NULL);
				atomName = trimWhitespace(record.substr(8, 8));
				atomType = record.substr(47, 2);
				if (atomType.compare(0, 1, " ") == 0){
					atomType = record.substr(48, 2);
				}
				atomType = reduceToElement(atomType);
				if (atomType.empty()){
					atomType = elementFromAtomName(atomName);
				}
				return true;
			}
		}

		vector<string> tokens;
		splitWhitespace(record, tokens);

		if (tokens.size() >= 5
			&& isPureNumber(tokens[2]) && isPureNumber(tokens[3]) && isPureNumber(tokens[4])){
			atomName = tokens[1];
			x = strtod(tokens[2].c_str(), NULL);
			y = strtod(tokens[3].c_str(), NULL);
			z = strtod(tokens[4].c_str(), NULL);
			atomType = (tokens.size() >= 6) ? reduceToElement(tokens[5])
				: elementFromAtomName(atomName);
			return true;
		}

		// Four fields (id x y z) carry no name/type; coordinates are usable.
		if (4 == tokens.size()
			&& isPureNumber(tokens[1]) && isPureNumber(tokens[2]) && isPureNumber(tokens[3])){
			x = strtod(tokens[1].c_str(), NULL);
			y = strtod(tokens[2].c_str(), NULL);
			z = strtod(tokens[3].c_str(), NULL);
			return true;
		}

		return false;
	}

	Ligand parseOneLigMol2(vector<string> oneLigMol2Infos, int blockStartLine){
		string ligname;
		vector<Point3D> points;
		vector<string> atomtypes;
		vector<string> atomNames;
		vector<BondInfo> bondInfoes;
		
		int l = 0;
		while (l < oneLigMol2Infos.size() && oneLigMol2Infos[l].compare(0, 13, "@<TRIPOS>ATOM") != 0){
			if (oneLigMol2Infos[l].compare(0, 17, "@<TRIPOS>MOLECULE") == 0){
				if (l + 1 < (int)oneLigMol2Infos.size()){
					ligname = oneLigMol2Infos[++l];
				}
			}
			
			l++;
		}
		
		if (l < (int)oneLigMol2Infos.size()){
			l++;
		}
		while (l < oneLigMol2Infos.size()){
			if (oneLigMol2Infos[l].compare(0, 9, "@<TRIPOS>") == 0){
				break;
			}
			
			// Empty strings separate sections; every other line is a record.
			if (oneLigMol2Infos[l].empty()){
				l++;
				continue;
			}

			string atomName;
			string atomtype;
			double x = 0.0;
			double y = 0.0;
			double z = 0.0;
			if (!parseAtomRecord(oneLigMol2Infos[l], atomName, atomtype, x, y, z)){
				if (oneLigMol2Infos[l].size() < 48){
					cerr << "LSalign: " << filePath << ":" << blockStartLine + l
						<< ": ATOM record too short for atom type column (ligand '"
						<< ligname << "')\n     " << oneLigMol2Infos[l] << endl;
				}else{
					cerr << "LSalign: " << filePath << ":" << blockStartLine + l
						<< ": cannot parse ATOM record (ligand '"
						<< ligname << "')\n     " << oneLigMol2Infos[l] << endl;
				}
				exit(EXIT_FAILURE);
			}

			if (!isFiniteCoord(x) || !isFiniteCoord(y) || !isFiniteCoord(z)){
				cerr << "LSalign: " << filePath << ":" << blockStartLine + l
					<< ": non-finite coordinate (ligand '" << ligname << "')" << "\n     "
					<< oneLigMol2Infos[l] << endl;
				exit(EXIT_FAILURE);
			}

			Point3D pos(x, y, z);
			points.push_back(pos);
			atomtypes.push_back(atomtype);
			atomNames.push_back(atomName);
			
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

			int nAtoms = (int)points.size();
			if (atom1Ind < 0 || atom1Ind >= nAtoms
				|| atom2Ind < 0 || atom2Ind >= nAtoms || atom1Ind == atom2Ind){
				cerr << "LSalign: " << filePath << ":" << blockStartLine + ll
					<< ": BOND record refers to atom index out of range (ligand '"
					<< ligname << "', atoms " << (atom1Ind+1) << "-" << (atom2Ind+1)
					<< ", valid 1-" << nAtoms << ")\n     " << oneLigMol2Infos[ll] << endl;
				exit(EXIT_FAILURE);
			}
			
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
			if (!(bondInfo.bondDisPow2 >= 0.0 && bondInfo.bondDisPow2 <= 1e308)){
				cerr << "LSalign: " << filePath << ":" << blockStartLine + ll
					<< ": non-finite bond length (ligand '" << ligname
					<< "', atoms " << (atom1Ind+1) << "-" << (atom2Ind+1) << ")" << "\n     "
					<< oneLigMol2Infos[ll] << endl;
				exit(EXIT_FAILURE);
			}
			bondInfo.bondDis = sqrt(bondInfo.bondDisPow2);
			
			bondInfoes.push_back(bondInfo);
			
			ll++;
		}
		
		if (!isLoad_H_atom){			
			vector<string>::iterator at_it;
			vector<Point3D>::iterator pt_it;
			vector<string>::iterator an_it;
			int *oldMapNewInd = new int[atomtypes.size()];
			int oldInd = 0, NewInd = 0;
			
			for (at_it = atomtypes.begin(), pt_it = points.begin(), an_it = atomNames.begin(); at_it != atomtypes.end();){
				if (0 == at_it->compare(0, at_it->size(), "H")
					|| 0 == at_it->compare(0, at_it->size(), "h")){
					at_it = atomtypes.erase(at_it);
					pt_it = points.erase(pt_it);
					an_it = atomNames.erase(an_it);
					
					oldMapNewInd[oldInd++] = -1;
				}else{
					at_it++;
					pt_it++;
					an_it++;
					
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
		
		if (points.empty()){
			cerr << "LSalign: warning: ligand '" << ligname << "' in " << filePath
				<< " has no atoms; its scores will be 0" << endl;
		}
		Ligand ans(ligname, points, atomtypes, bondInfoes, atomNames);
		return ans;
	}

	/***********************************************
	* load all the ligand information in the mol2 file
	* using the string container
	***********************************************/
	vector<vector<string> > loadAllLigStrInfoMol2(vector<int>& blockStartLines){
		vector<vector<string> > ans;
		
		bool isFirst = true;
		vector<string> mol2Sb;
		ifstream fin(filePath.c_str());
		string line;
		int lineNo = 1;			// 1-based file line number of `line`
		int curBlockStartLine = 1;	// file line number of mol2Sb[0]
		getline(fin, line);
		while (fin.good()){
			line = eraseAll(line, '\r');	// Since the differnece between windows (\r\n) and Linux (\n), so we should do this step.
			if (line.compare(0, 17, "@<TRIPOS>MOLECULE") == 0){
				if (isFirst){
					isFirst = false;
				}else{
					ans.push_back(mol2Sb);
					blockStartLines.push_back(curBlockStartLine);
					mol2Sb.erase(mol2Sb.begin(), mol2Sb.end());
				}
			}
			
			// mol2Sb[0] may precede the first MOLECULE (comment header); this
			// keeps the index->file-line mapping exact for every block.
			curBlockStartLine = lineNo - (int)mol2Sb.size();
			
			mol2Sb.push_back(line);
			getline(fin, line);
			lineNo++;
			if (!fin.good()){
				mol2Sb.push_back(line);
			}
		}
		fin.close();
		
		if (!isFirst){
			ans.push_back(mol2Sb);
			blockStartLines.push_back(curBlockStartLine);
		}
		return ans;
	}
};

/******************************************
* Ligands -- thin loader delegating to a format parser.
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

		LigandParser* parser = createLigandParser(mol2path);
		ligs = parser->parse(mol2path, isLoad_H_atom);
		delete parser;
	}

	
	~Ligands(){
		vector<Ligand>().swap(ligs);
	}	
	
	vector<Ligand> getLigands(){
		return ligs;
	}
	
	const Ligand& operator[](int ind) const{
		return ligs[ind];
	}
	
	int size(){
		return ligs.size();
	}
private:
	static LigandParser* createLigandParser(const string& path){
		string ext;
		size_t dot = path.find_last_of('.');
		if (string::npos != dot){
			ext = path.substr(dot + 1);
			toLowerString(ext);
		}
		if (ext.empty() || ext == "mol2"){
			return new Mol2Parser();
		}
		cerr << "LSalign: warning: unsupported file extension "
			<< ext << ", falling back to mol2 parser." << endl;
		return new Mol2Parser();
	}
};

#endif
