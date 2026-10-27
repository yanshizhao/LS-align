#ifndef RIGIDLSALIGN_H
#define RIGIDLSALIGN_H

#include "basic_fun.h"
#include "Ligands.h"
#include "AtomFixInfo.h"
#include "EnhancedGS.h"
#include "Point3D.h"
#include "Kabsch.h"
#include <cmath>
#include "Evaluation.h"
#include <map>
using namespace std;

//#define OU_SCO_WEIGHT 0.5
//#define MASS_SCO_WEIGHT 0.05
//#define JUNBOND_SCO_WEIGHT 0.45

#define OU_SCO_WEIGHT 0.45
#define MASS_SCO_WEIGHT 0.1
#define JUNBOND_SCO_WEIGHT 0.45

#define OU_SCO_WEIGHT_FOR_FINAL_SCO 0.45
#define MASS_SCO_WEIGHT_FOR_FINAL_SCO 0.1
#define JUNBOND_SCO_WEIGHT_FOR_FINAL_SCO 0.45


#define delta_m0_pow2 100		// 10 * 10
#define delta_vdwR0_pow2 0.36		// 0.6 * 0.6, 1.2 is the var der waal radii of atom H

class RigidLSalign{
private:
	Ligand	query;		// this Ligand will be rotated to make an good alignment with the templ
	Ligand	templ;		// this is a template Ligand
	double 	d0_pow2;		// the value is the square of d0 (angstrom^2)
	int 	rough_sliding_step;	// the rough sliding step default is 5
	int 	detail_sliding_step;	// the detail sliding step default is 1

	double alisco;					// 1/alinum * sum_{ali_{i}} (1 / (1 + (d_{ali{i}})^2 / (d0)^2))

	AlignInfo* p_AliInfo;			// the last alignment result
	double rotMtx[3][3];			// the last rotation matrix
	double transVec[3];				// the last transform vector

	vector<AlignedPair> finalAlignedPair;	// the final aligned information
	vector<double> finalAlignedPairDist;	// record the distance of the final aligned pairs after superposed

	// temp information
	vector<vector<double> > massDiffScoMtx;
	vector<vector<double> > vdwRDiffScoMtx;
	vector<vector<double> > bondInfoScoMtx;

	vector<vector<double> > junBondInfoScoMtx;
	vector<vector<double> > junBondInfoScoMtx_for_final_sco;
public: 
	RigidLSalign(Ligand query, Ligand templ){
		this->query = query;
		this->templ = templ;
		
		this->d0_pow2 = 4.0;	// 2*2
		this->rough_sliding_step = 5;	// default is 5
		this->detail_sliding_step = 1;
		this->p_AliInfo = NULL;
		
		alisco = -1;

		align();
	}

	RigidLSalign(const RigidLSalign& rLSalign){
		this->query = rLSalign.query;
		this->templ = rLSalign.templ;
		this->d0_pow2 = rLSalign.d0_pow2;
		this->rough_sliding_step = rLSalign.rough_sliding_step;
		this->p_AliInfo = new AlignInfo(*rLSalign.p_AliInfo);

		int i = 0, j = 0;
		for (i = 0; i < 3; i++){
			for (j = 0; j < 3; j++){
				this->rotMtx[i][j] = rLSalign.rotMtx[i][j];
			}
			this->transVec[i] = rLSalign.transVec[i];
		}
		
		for (i = 0; i < rLSalign.finalAlignedPair.size(); i++){
			this->finalAlignedPair.push_back( rLSalign.finalAlignedPair[i]);
		}

		for (i = 0; i < rLSalign.finalAlignedPairDist.size(); i++){
			this->finalAlignedPairDist.push_back( rLSalign.finalAlignedPairDist[i]);
		}

		alisco = rLSalign.alisco;

		for (i = 0; i < rLSalign.massDiffScoMtx.size(); i++){
			vector<double> tmp_mass_row;
			for (j = 0; j < rLSalign.massDiffScoMtx[i].size(); j++){
				tmp_mass_row.push_back(rLSalign.massDiffScoMtx[i][j]);
			}

			massDiffScoMtx.push_back(tmp_mass_row);
		}

		for (i = 0; i < rLSalign.vdwRDiffScoMtx.size(); i++){
			vector<double> tmp_vdwR_row;
			for (j = 0; j < rLSalign.vdwRDiffScoMtx[i].size(); j++){
				tmp_vdwR_row.push_back(rLSalign.vdwRDiffScoMtx[i][j]);
			}

			vdwRDiffScoMtx.push_back(tmp_vdwR_row);
		}

		for (i = 0; i < rLSalign.bondInfoScoMtx.size(); i++){
			vector<double> tmp_bondInfo_row;
			for (j = 0; j < rLSalign.bondInfoScoMtx[i].size(); j++){
				tmp_bondInfo_row.push_back(rLSalign.bondInfoScoMtx[i][j]);
			}

			bondInfoScoMtx.push_back(tmp_bondInfo_row);
		}

		for (i = 0; i < rLSalign.junBondInfoScoMtx.size(); i++){
			vector<double> tmp_junhbondInfo_row;
			for (j = 0; j < rLSalign.junBondInfoScoMtx[i].size(); j++){
				tmp_junhbondInfo_row.push_back(rLSalign.junBondInfoScoMtx[i][j]);
			}

			junBondInfoScoMtx.push_back(tmp_junhbondInfo_row);
		}

		for (i = 0; i < rLSalign.junBondInfoScoMtx_for_final_sco.size(); i++){
			vector<double> tmp_junhbondInfo_for_final_sco_row;
			for (j = 0; j < rLSalign.junBondInfoScoMtx_for_final_sco[i].size(); j++){
				tmp_junhbondInfo_for_final_sco_row.push_back(rLSalign.junBondInfoScoMtx_for_final_sco[i][j]);
			}

			junBondInfoScoMtx_for_final_sco.push_back(tmp_junhbondInfo_for_final_sco_row);
		}
	}

	RigidLSalign(Ligand query, Ligand templ, double d0, 
			int rough_sliding_step, int detail_sliding_step){
		this->query = query;
		this->templ = templ;
		this->d0_pow2 = d0*d0;
		this->rough_sliding_step = rough_sliding_step;
		this->detail_sliding_step = detail_sliding_step;
		
		this->p_AliInfo = NULL;


		alisco = -1;

		align();
	}

	RigidLSalign(Ligand query, Ligand templ, double d0){
		this->query = query;
		this->templ = templ;
		this->d0_pow2 = d0*d0;
		this->rough_sliding_step = 5;	// default
		this->detail_sliding_step = 1;	// default
		
		this->p_AliInfo = NULL;

		alisco = -1;

		align();
	}

	~RigidLSalign(){
		int i = 0;

		if (NULL == p_AliInfo){
			delete p_AliInfo;
			p_AliInfo = NULL;
		}

		vector<AlignedPair>().swap(finalAlignedPair);
		vector<double>().swap(finalAlignedPairDist);

		for (i = 0; i < massDiffScoMtx.size(); i++){
			massDiffScoMtx[i].clear();
			vector<double>().swap(massDiffScoMtx[i]);
		}

		massDiffScoMtx.clear();
		vector<vector<double> >().swap(massDiffScoMtx);
		for (i = 0; i < vdwRDiffScoMtx.size(); i++){
			vdwRDiffScoMtx[i].clear();
			vector<double>().swap(vdwRDiffScoMtx[i]);
		}

		vdwRDiffScoMtx.clear();
		vector<vector<double> >().swap(vdwRDiffScoMtx);

		for (i = 0; i < bondInfoScoMtx.size(); i++){
			bondInfoScoMtx[i].clear();
			vector<double>().swap(bondInfoScoMtx[i]);
		}

		bondInfoScoMtx.clear();
		vector<vector<double> >().swap(bondInfoScoMtx);
	}

	vector<string> toStrVec(){
		if (0 == finalAlignedPair.size()){
			parseFinalAliInfo();
		}

		vector<string> ans;
		
		char buf[2000];
		sprintf(buf, "The atom number of the Query Ligand is %d", query.size());
		string firstLine(buf);
		ans.push_back(firstLine);
		
		sprintf(buf, "The atom number of the Templ Ligand is %d", templ.size());
		string secondLine(buf);
		ans.push_back(secondLine);
		
		sprintf(buf, "The LS-score in the LS-align is %.3f", p_AliInfo->sco/(templ.size()>query.size()?templ.size():query.size()));
		string thirdLine(buf);
		ans.push_back(thirdLine);

		sprintf(buf, "The AlignedRMSD in the LS-align is %.3f", getAlignedRMSD());
		string fourthLine(buf);
		ans.push_back(fourthLine);

		sprintf(buf, "The JaccardRatio in the LS-align is %.3f", getJaccardRatio());
		string fifthLine(buf);
		ans.push_back(fifthLine);

		sprintf(buf, "The RMSD_lb in the LS-align is %.3f", getRMSD_lb());
		string sixthLine(buf);
		ans.push_back(sixthLine);

		sprintf(buf, "The Aligned number in the LS-align is %d", finalAlignedPairDist.size());
		string seventhLine(buf);
		ans.push_back(seventhLine);
		
		for (int i = 0; i < finalAlignedPair.size(); i++){
			AlignedPair ithPair = finalAlignedPair[i];
			double dis = finalAlignedPairDist[i];
			sprintf(buf, "%d - %d : %.3f", ithPair.rowInd, ithPair.colInd, dis);
			string line(buf);
			ans.push_back(line);
		}
		
		return ans;
	}

	double** getRotMtx(){	// need release the return ans
		double** u = new double*[3];
		for (int i = 0; i < 3; i++){
			u[i] = new double[3];
			for (int j = 0; j < 3; j++){
				u[i][j] = rotMtx[i][j];
			}
		}
		
		return u;
	}

	double* getTransVec(){
		double* t = new double[3];
		for (int i = 0; i < 3; i++){
			t[i] = transVec[i];
		}
		
		return t;
	}

	vector<string> formatRotMtx(){
		vector<string> ans;
		
		char buf[1000];
		sprintf(buf, "i%18s %15s %15s %15s", "t[i]", "u[i][0]", "u[i][1]", "u[i][2]");
		ans.push_back(string(buf));
		for (int k = 0; k < 3; k++){
			sprintf(buf, "%d%18.10f %15.10f %15.10f %15.10f", k, transVec[k], rotMtx[k][0], rotMtx[k][1], rotMtx[k][2]);
			ans.push_back(string(buf));
		}
		ans.push_back(string("\nCode for rotating Structure Query from (x,y,z) to (X,Y,Z):"));
		ans.push_back(string("for(k=0; k<L; k++)"));
		ans.push_back(string("{"));
		ans.push_back(string("   X[k] = t[0] + u[0][0]*x[k] + u[0][1]*y[k] + u[0][2]*z[k]"));
		ans.push_back(string("   Y[k] = t[1] + u[1][0]*x[k] + u[1][1]*y[k] + u[1][2]*z[k]"));
		ans.push_back(string("   Z[k] = t[2] + u[2][0]*x[k] + u[2][1]*y[k] + u[2][2]*z[k]"));
		ans.push_back(string("}"));

		return ans;
	}

	vector<string> formatAliInfo(){
		if (0 == finalAlignedPair.size()){
			parseFinalAliInfo();
		}
		/* double d0 = sqrt(d0_pow2); */

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
			double dis = finalAlignedPairDist[i];

			char tag[2];
			tag[0] = '.';
			if (dis < 1.0) // if (dis < d0/2.0)
				tag[0] = ':';
			tag[1] = '\0';
			
			string tAT = templ.getAtomType(ithPair.colInd);
			int tAI = ithPair.colInd + 1;		// index start from 1

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

	Ligand getRottedQueryLigand(){
		Ligand ans(query);
		vector<Point3D> rotPts = rotateAndTrans(ans.getPoints(), transVec, rotMtx);
		ans.setPoints(rotPts);
		
		return ans;
	}
	
	double get_d0(){
		return sqrt(d0_pow2);
	}

	double getLSscore(){
		return p_AliInfo->sco/(templ.size()>query.size()?templ.size():query.size());
	}
	
	double getLSscore8QuerySize(){
		return p_AliInfo->sco / query.size();
	}

	double getLSscore8TemplSize(){
		return p_AliInfo->sco / templ.size();
	}

	vector<AlignedPair> getFinalAlignedPair(){
		if (0 == finalAlignedPair.size()){
			parseFinalAliInfo();
		}	
		
		return finalAlignedPair;
	}

	vector<double> getFinalAlignedPairDist(){
		if (0 == finalAlignedPairDist.size()){
			parseFinalAliInfo();
		}	
		
		return finalAlignedPairDist;
	}

	double getJaccardRatio(){
		if (0 == finalAlignedPairDist.size()){
			parseFinalAliInfo();
		}
		double ans = 1.0 * finalAlignedPairDist.size() / ((int)templ.size() + (int)query.size() - (int)finalAlignedPairDist.size());

		return ans;
	}

	double getAlignedRMSD(){
		if (0 == finalAlignedPairDist.size()){
			parseFinalAliInfo();
		}
		
		double ans = 0.0;		
		for (int i = 0; i < finalAlignedPairDist.size(); i++){
			ans += finalAlignedPairDist[i] * finalAlignedPairDist[i];
		}
		ans = sqrt(ans / finalAlignedPairDist.size());		
		
		return ans;
	}

	double getRMSD_lb(){
		vector<Point3D> rotQueryPts = rotateAndTrans(query.getPoints(), transVec, rotMtx);
		return RMSD_lb(rotQueryPts, templ.getPoints());
	}
	
	double getRMSD_lb_min(){
		vector<Point3D> rotQueryPts = rotateAndTrans(query.getPoints(), transVec, rotMtx);
		return RMSD_lb_min(rotQueryPts, templ.getPoints());
	}

	double getAliSco(){
		if (-1 == alisco){
			parseFinalAliInfo();
		}

		return alisco;
	}

private:
	void parseFinalAliInfo(){
		alisco = 0.0;
		int alinum = 0;

		// junBondInfoScoMtx_for_final_sco
		double sco = 0.0;

		vector<AlignedPair> alignedPaired = p_AliInfo->getAliPairs();
		for (int i = 0; i < alignedPaired.size(); i++){
			AlignedPair ithPair = alignedPaired[i];
			
			Point3D rotQueryPt = rotateAndTrans(query.getIthPoints(ithPair.rowInd), transVec, rotMtx);
			double disPow2 = rotQueryPt.distPow2(templ.getIthPoints(ithPair.colInd));
			double ithsco = OU_SCO_WEIGHT_FOR_FINAL_SCO * (1.0 / (1.0 + disPow2/d0_pow2))
				+ MASS_SCO_WEIGHT_FOR_FINAL_SCO * massDiffScoMtx[ithPair.rowInd][ithPair.colInd]
				+ JUNBOND_SCO_WEIGHT_FOR_FINAL_SCO * junBondInfoScoMtx_for_final_sco[ithPair.rowInd][ithPair.colInd];

			sco += ithsco;

			if (disPow2 > d0_pow2)			
				continue;
			finalAlignedPair.push_back(ithPair);
			finalAlignedPairDist.push_back(sqrt(disPow2));

			alinum++;
			alisco += ithsco;
		}

		if (0 != alinum){
			alisco /= alinum;
		}
		
		p_AliInfo->sco = sco;
	}	

	void align(){
		calculatemassDiffScoMtx();
		calculateVdwRDiffScoMtx();
		calculateBondInfoScoMtx();
		calculateJunBondInfoScoMtx();

		double maxInitLSsco = -999999;
		AlignInfo* p_finalInitAli = NULL;
		{
			// ===================================================================
			// using rough sliding_step to search a mass to obatin the initial 
			// alignment
			// ===================================================================
			EnhancedGS massInitEgs(massDiffScoMtx);
			AlignInfo* p_massInitAli = massInitEgs.getPAliInfo();
			p_massInitAli->sco = 0.0;
			
			double u[3][3];
			double* t = new double[3];
			AlignInfo* p_massBestAli = nonsequentialIterAlign(p_massInitAli, rough_sliding_step, 3, true, u, t, false, 10);
			if (maxInitLSsco < p_massBestAli->sco){
				maxInitLSsco = p_massBestAli->sco;
				
				for (int h = 0; h < 3; h++){
					transVec[h] = t[h];
					for (int d = 0; d < 3; d++){
						rotMtx[h][d] = u[h][d];
					}
				}
				
				if (NULL != p_finalInitAli){
					delete p_finalInitAli;
					p_finalInitAli = NULL;
				}
				p_finalInitAli = p_massBestAli;
			}else{
				delete p_massBestAli;
				p_massBestAli = NULL;
			}
		
			delete[] t;
			t = NULL;
		}
		
		{
			// ===================================================================
			// using rough sliding_step (5) to search a atom VDW radius difference 
			// to obatin the initial alignment.
			// ===================================================================
			EnhancedGS vdwrInitEgs(vdwRDiffScoMtx);
			AlignInfo* p_vdwrInitAli = vdwrInitEgs.getPAliInfo();
			p_vdwrInitAli->sco = 0.0;			

			double u[3][3];
			double* t = new double[3];
			AlignInfo* p_vdwrBestAli = nonsequentialIterAlign(p_vdwrInitAli, rough_sliding_step, 3, true, u, t, false, 10);
			if (maxInitLSsco < p_vdwrBestAli->sco){
				maxInitLSsco = p_vdwrBestAli->sco;
				
				for (int h = 0; h < 3; h++){
					transVec[h] = t[h];
					for (int d = 0; d < 3; d++){
						rotMtx[h][d] = u[h][d];
					}
				}
				
				if (NULL != p_finalInitAli){
					delete p_finalInitAli;
					p_finalInitAli = NULL;
				}
				p_finalInitAli = p_vdwrBestAli;
			}else{
				delete p_vdwrBestAli;
				p_vdwrBestAli = NULL;
			}
			
			delete[] t;
			t = NULL;
		}

		{
			// ===================================================================
			// using rough sliding_step (5) to search a bond information matrix
			// to obatin the initial alignment.
			// ===================================================================
			EnhancedGS biInitEgs(bondInfoScoMtx);
			AlignInfo* p_biInitAli = biInitEgs.getPAliInfo();
			p_biInitAli->sco = 0.0;			

			double u[3][3];
			double* t = new double[3];
			AlignInfo* p_biBestAli = nonsequentialIterAlign(p_biInitAli, rough_sliding_step, 3, true, u, t, false, 10);
			if (maxInitLSsco < p_biBestAli->sco){
				maxInitLSsco = p_biBestAli->sco;
				
				for (int h = 0; h < 3; h++){
					transVec[h] = t[h];
					for (int d = 0; d < 3; d++){
						rotMtx[h][d] = u[h][d];
					}
				}
				
				if (NULL != p_finalInitAli){
					delete p_finalInitAli;
					p_finalInitAli = NULL;
				}
				p_finalInitAli = p_biBestAli;
			}else{
				delete p_biBestAli;
				p_biBestAli = NULL;
			}
			
			delete[] t;
			t = NULL;
		}
		
		{
			// ===================================================================
			// using rough sliding_step (5) to search a jun's definition bond information matrix
			// to obatin the initial alignment.
			// ===================================================================
			EnhancedGS junhBiInitEgs(bondInfoScoMtx);
			AlignInfo* p_junBiInitAli = junhBiInitEgs.getPAliInfo();
			p_junBiInitAli->sco = 0.0;			

			double u[3][3];
			double* t = new double[3];
			AlignInfo* p_junBiBestAli = nonsequentialIterAlign(p_junBiInitAli, rough_sliding_step, 3, true, u, t, false, 10);
			if (maxInitLSsco < p_junBiBestAli->sco){
				maxInitLSsco = p_junBiBestAli->sco;
				
				for (int h = 0; h < 3; h++){
					transVec[h] = t[h];
					for (int d = 0; d < 3; d++){
						rotMtx[h][d] = u[h][d];
					}
				}
				
				if (NULL != p_finalInitAli){
					delete p_finalInitAli;
					p_finalInitAli = NULL;
				}
				p_finalInitAli = p_junBiBestAli;
			}else{
				delete p_junBiBestAli;
				p_junBiBestAli = NULL;
			}
			
			delete[] t;
			t = NULL;
		}
		
		this->p_AliInfo = nonsequentialIterAlign(p_finalInitAli, detail_sliding_step, 6, true, this->rotMtx, this->transVec, true, 1);
		
		parseFinalAliInfo();

		delete p_finalInitAli;
		p_finalInitAli = NULL;
	}
		
	AlignInfo* nonsequentialIterAlign(AlignInfo* p_initAliInfo, int sliding_step, int max_win_type_num, 
						bool isNeedRotMtxAndTransVec, double rotMtx[3][3], double transVec[3], bool isDirectlyInitAli, int max_iter_num){
		int iter = 0;
		AlignInfo* p_bestAliInfo = NULL;
		if (isDirectlyInitAli){
			p_bestAliInfo = new AlignInfo(*p_initAliInfo);
		}else{
			p_bestAliInfo = nonsequentialUpdateAliInfo(p_initAliInfo, sliding_step, max_win_type_num, isNeedRotMtxAndTransVec, rotMtx, transVec);
			iter = 1;
		}
		 
		for (; iter < max_iter_num; iter++){
			double u[3][3];
			double* t = new double[3];
			
			AlignInfo* p_iterAliInfo = nonsequentialUpdateAliInfo(p_bestAliInfo, sliding_step, max_win_type_num, isNeedRotMtxAndTransVec, u, t);
			
			if (p_iterAliInfo->sco > p_bestAliInfo->sco){
				delete p_bestAliInfo;
				p_bestAliInfo = NULL;
				p_bestAliInfo = p_iterAliInfo;
				p_iterAliInfo = NULL;
				
				if (isNeedRotMtxAndTransVec){
					for (int h = 0; h < 3; h++){
						transVec[h] = t[h];
						for (int d = 0; d < 3; d++){
							rotMtx[h][d] = u[h][d];
						}
					}
				}
			}else{
				delete p_iterAliInfo;
				p_iterAliInfo = NULL;			
			
				delete[] t;
				t = NULL;
				break;
			}
			
			// release and t
			delete[] t;
			t = NULL;
		}
		
		return p_bestAliInfo;
	}
	
	AlignInfo* nonsequentialUpdateAliInfo(AlignInfo* p_oldAliInfo, int sliding_step, int max_win_type_num, 
						bool isNeedRotMtxAndTransVec, double rotMtx[3][3], double transVec[3]){
		vector<AlignedPair> oldAli = p_oldAliInfo->getAliPairs();
		
		int n_ali = oldAli.size();
		int L_min = n_ali > 4 ? 4 : n_ali;
		vector<int> L_init_vec;
		if (L_min <= 4){
			L_init_vec.push_back(L_min);
		}else{
			int win_type_num = 0;
			for (int i = 0; ; i++){
				int fenmu = (int)pow(2.0, i);
	
				int L_init = n_ali / fenmu;
				if (L_init > L_min && win_type_num < max_win_type_num){
					L_init_vec.push_back(L_init);
				}else{
					L_init_vec.push_back(L_min);
					break;
				}
				
				win_type_num++;
			}
		}
		
		double maxScore = -999999;
		AlignInfo* p_corrAliInfo = NULL;
		
		for (int i = 0; i < L_init_vec.size(); i++) {
			int L_init = L_init_vec[i];		// the ali atom number
			int iL_max = n_ali - L_init + 1;	// slide number 
			for (int iL = 0; iL < iL_max; iL+=sliding_step) {	// each sliding
				vector<Point3D> queryali;
				vector<Point3D> templali;
				
				for (int ii = 0; ii < L_init; ii++) {
					int k = iL + ii;
					
					queryali.push_back(query.getIthPoints(oldAli[k].rowInd));
					templali.push_back(templ.getIthPoints(oldAli[k].colInd));
				}
				
				double u[3][3];
				double* t = new double[3];
				Kabsch(queryali, templali, t, u);
				
				vector<Point3D> queryPoses_T = rotateAndTrans(query.getPoints(), t, u);
				AlignInfo* pAliInfo = LSscoreFunc(queryPoses_T, templ.getPoints());	// the return new a new memory
				
				if (pAliInfo->sco > maxScore){
					maxScore = pAliInfo->sco;
					
					if (isNeedRotMtxAndTransVec){
						for (int h = 0; h < 3; h++){
							transVec[h] = t[h];
							for (int d = 0; d < 3; d++){
								rotMtx[h][d] = u[h][d];
							}
						}				
					}
										
					if (NULL != p_corrAliInfo){
						delete p_corrAliInfo;
						p_corrAliInfo = NULL;
					}
					p_corrAliInfo = pAliInfo;
					pAliInfo = NULL;
				}else{
					delete pAliInfo;
					pAliInfo = NULL;
				}
				
				// release memory
				delete[] t;
				t = NULL;
			}
		}
		
		return p_corrAliInfo;
	}
	
	AlignInfo* LSscoreFunc(vector<Point3D> query_rotted_pts, vector<Point3D> templ_pts){
		vector<vector<double> > itemScoreMtx;
		
		for (int i = 0; i < query_rotted_pts.size(); i++){
			vector<double> qRow;
			for (int j = 0; j < templ_pts.size(); j++){
				double disPow2 = query_rotted_pts[i].distPow2(templ_pts[j]);
				double disSco_ij = 1.0/(1.0 + disPow2 / d0_pow2);
			    
				/*
				  OU_SCO_WEIGHT       几何重叠的权重
				  MASS_SCO_WEIGHT     质量相似的权重
				  JUNBOND_SCO_WEIGHT  键环境相似的权重
				*/
				qRow.push_back(OU_SCO_WEIGHT*disSco_ij 
					+ MASS_SCO_WEIGHT*massDiffScoMtx[i][j]
					+ JUNBOND_SCO_WEIGHT * junBondInfoScoMtx[i][j]);
			}

			itemScoreMtx.push_back(qRow);
		}
		
		EnhancedGS egs(itemScoreMtx);
		AlignInfo* tmpans = egs.getPAliInfo();	// the tmpans point is include in EnhancedGS
		AlignInfo* pans = new AlignInfo(*tmpans);
		return pans;
	}
	
	void calculatemassDiffScoMtx(){
		AtomRelativeMass arm;		

		for (int i = 0; i < query.size(); i++){
			string qAtomtype = query.getAtomType(i);	// the atom type is upperletter
			double qrm = arm[qAtomtype];			// get the query atom relative mass
			
			vector<double> qRow;
			for (int j = 0; j < templ.size(); j++){
				string tAtomtype = templ.getAtomType(j);	// the atom type is upperletter
				double trm = arm[tAtomtype];			// get the templ atom relative mass

				double absDiff = abs(qrm - trm);
				double absDiffSco = 1.0/(1.0 + (absDiff*absDiff / delta_m0_pow2));
				qRow.push_back(absDiffSco);
			}

			massDiffScoMtx.push_back(qRow);
		}
	}

	void calculateVdwRDiffScoMtx(){
		AtomVarDerWaalRadius avdwr;

		for (int i = 0; i < query.size(); i++){
			string qAtomtype = query.getAtomType(i);	// the atom type is upperletter
			double qvdwr = avdwr[qAtomtype];			// get the query atom var der waal radiue
			
			vector<double> qRow;
			for (int j = 0; j < templ.size(); j++){
				string tAtomtype = templ.getAtomType(j);	// the atom type is upperletter
				double tvdwr = avdwr[tAtomtype];			// get the templ atom var der waal radiue

				double absDiff = abs(qvdwr - tvdwr);
				double absDiffSco = 1.0/(1.0 + (absDiff*absDiff / delta_vdwR0_pow2));
				qRow.push_back(absDiffSco);
			}

			vdwRDiffScoMtx.push_back(qRow);
		}
	}

	void calculateJunBondInfoScoMtx(){
		vector< vector<BondType> > queryJunAtomBondTypeMtx = query.getAtomBondTypes();
		vector< vector<BondType> > templJunAtomBondTypeMtx = templ.getAtomBondTypes();

		vector<vector<double> > _junBondInfoScoMtx;
		double maxValEachElement = -9999999.0;
		double minValEachElement =  9999999.0;
		for (int i = 0; i < query.size(); i++){
			vector<BondType> ithQAtom = queryJunAtomBondTypeMtx[i];

			double ithQAtomScore_for_final_sco = 0.0;
			for (int l = 0; l < ithQAtom.size(); l++){
				ithQAtomScore_for_final_sco += score4EachBondTypePair(ithQAtom[l], ithQAtom[l]);
			}

			vector<double> qRow;
			vector<double> qRow_for_final_sco;
			for (int j = 0; j < templ.size(); j++){
				vector<BondType> jthTAtom = templJunAtomBondTypeMtx[j];				

				double jthTAtomScore_for_final_sco = 0.0;
				for (int m = 0; m < jthTAtom.size(); m++){
					jthTAtomScore_for_final_sco += score4EachBondTypePair(jthTAtom[m], jthTAtom[m]);
				}

				vector<vector<double> > itemScoMtx;
				for (int l = 0; l < ithQAtom.size(); l++){
					vector<double> row;
					for (int m = 0; m < jthTAtom.size(); m++){
						double lm_sco = score4EachBondTypePair(ithQAtom[l], jthTAtom[m]);

						row.push_back(lm_sco);
					}
					itemScoMtx.push_back(row);
				}
				
			
				EnhancedGS egs(itemScoMtx);
				AlignInfo* p_Ali = egs.getPAliInfo();
				double overlapScore = (p_Ali->sco) / (1.0 + abs(1.0 * ((int)ithQAtom.size() - (int)jthTAtom.size())));
				qRow.push_back(overlapScore);

				if (maxValEachElement < overlapScore)
					maxValEachElement = overlapScore;
				if (minValEachElement > overlapScore)
					minValEachElement = overlapScore;

				double overlapScore_for_final_sco = p_Ali->sco;
				double jaccardRatio_for_final_sco = 1.0*overlapScore_for_final_sco / 
					(ithQAtomScore_for_final_sco + jthTAtomScore_for_final_sco - overlapScore_for_final_sco);
				qRow_for_final_sco.push_back(jaccardRatio_for_final_sco);
			}// end for j

			_junBondInfoScoMtx.push_back(qRow);
			junBondInfoScoMtx_for_final_sco.push_back(qRow_for_final_sco);
		}// end for i

		for (int jun = 0; jun < _junBondInfoScoMtx.size(); jun++){
			vector<double> qRow;
			for (int dan = 0; dan < _junBondInfoScoMtx[jun].size(); dan++){
				double element = _junBondInfoScoMtx[jun][dan];
				if (maxValEachElement > minValEachElement){
					 element = (element - minValEachElement) / (maxValEachElement - minValEachElement);
				}
				qRow.push_back(element);
			}
			junBondInfoScoMtx.push_back(qRow);
		}
	}

	/***************************************************************************************
	 * This function is a sub function for "calculateJunBondInfoScoMtx" function
	 ***************************************************************************************/
	double score4EachBondTypePair(BondType bt1, BondType bt2){
		double ans = 0.0;
		if (bt1 == bt2 && bt1 != NOTCONNECTED){
			switch (bt1){
			case UNKNOWN:				ans = 0.8;	break;
			case SINGLE:				ans = 1.0;	break;
			case DOUBLE:				ans = 1.2;	break;
			case TRIPLE:				ans = 1.4;	break;
			case DUMMY:					ans = 1.6;	break;
			case AMIDE:					ans = 1.8;	break;
			case AROMATIC:				ans = 2.0;	break;
			
			// Three bond ring
			case UNKNOWN_RING_THREE:	ans = 2.8;	break;
			case SINGLE_RING_THREE:		ans = 3.0;	break;
			case DOUBLE_RING_THREE:		ans = 3.2;	break;
			case TRIPLE_RING_THREE:		ans = 3.4;	break;
			case AMIDE_RING_THREE:		ans = 3.6;	break;
			case AROMATIC_RING_THREE:	ans = 3.8;	break;
			case DUMMY_RING_THREE:		ans = 4.0;	break;

			case UNKNOWN_RING_FOUR:		ans = 4.8;	break;
			case SINGLE_RING_FOUR:		ans = 5.0;	break;
			case DOUBLE_RING_FOUR:		ans = 5.2;	break;
			case TRIPLE_RING_FOUR:		ans = 5.4;	break;
			case AMIDE_RING_FOUR:		ans = 5.6;	break;
			case AROMATIC_RING_FOUR:	ans = 5.8;	break;
			case DUMMY_RING_FOUR:		ans = 6.0;	break;
			
			case UNKNOWN_RING_FIVE:		ans = 6.8;	break;
			case SINGLE_RING_FIVE:		ans = 7.0;	break;
			case DOUBLE_RING_FIVE:		ans = 7.2;	break;
			case TRIPLE_RING_FIVE:		ans = 7.4;	break;
			case AMIDE_RING_FIVE:		ans = 7.6;	break;
			case AROMATIC_RING_FIVE:	ans = 7.8;	break;
			case DUMMY_RING_FIVE:		ans = 8.0;	break;

			case UNKNOWN_RING_SIX:		ans = 8.8;	break;
			case SINGLE_RING_SIX:		ans = 9.0;	break;
			case DOUBLE_RING_SIX:		ans = 9.2;	break;
			case TRIPLE_RING_SIX:		ans = 9.4;	break;
			case AMIDE_RING_SIX:		ans = 9.6;	break;
			case AROMATIC_RING_SIX:		ans = 9.8;	break;
			case DUMMY_RING_SIX:		ans = 10.0;	break;

			default:
				ans = 0.6;
			}
		}else{
			ans = 0.5;
		}

		return ans;
	}

	void calculateBondInfoScoMtx(){
		int i = 0;

		vector<BondInfo> queryBondOrder = query.getBondInfoes();
		vector<BondInfo> templBondOrder = templ.getBondInfoes();		

		map<int, vector<BondType> > queryBondTypeInfo;
		map<int, vector<double> > queryBondDisInfo;
		for (i = 0; i < queryBondOrder.size(); i++){
			BondInfo ithBI = queryBondOrder[i];
			
			// fill BondType information into map
			vector<BondType> firstTmp;
			map<int, vector<BondType> >::iterator it1 = queryBondTypeInfo.find(ithBI.aInd);			
			if (it1 != queryBondTypeInfo.end()){
				firstTmp = it1->second;
			}
			firstTmp.push_back(ithBI.bt);
			queryBondTypeInfo[ithBI.aInd] = firstTmp;

			vector<BondType> secondTmp;
			map<int, vector<BondType> >::iterator it2 = queryBondTypeInfo.find(ithBI.bInd);			
			if (it2 != queryBondTypeInfo.end()){
				secondTmp = it2->second;
			}
			secondTmp.push_back(ithBI.bt);
			queryBondTypeInfo[ithBI.bInd] = secondTmp;

			// fill distance pow2 information into map
			vector<double> firstDis;
			map<int, vector<double> >::iterator dit1 = queryBondDisInfo.find(ithBI.aInd);
			if (dit1 != queryBondDisInfo.end()){
				firstDis = dit1->second;
			}
			firstDis.push_back(ithBI.bondDis);
			queryBondDisInfo[ithBI.aInd] = firstDis;
			
			vector<double> secondDis;
			map<int, vector<double> >::iterator dit2 = queryBondDisInfo.find(ithBI.bInd);
			if (dit2 != queryBondDisInfo.end()){
				secondDis = dit2->second;
			}
			secondDis.push_back(ithBI.bondDis);
			queryBondDisInfo[ithBI.bInd] = secondDis;
		}		

		map<int, vector<BondType> > templBondTypeInfo;
		map<int, vector<double> > templBondDisInfo;
		for (i = 0; i < templBondOrder.size(); i++){
			BondInfo ithBI = templBondOrder[i];
			
			// fill BondType information into map
			vector<BondType> firstTmp;
			map<int, vector<BondType> >::iterator it1 = templBondTypeInfo.find(ithBI.aInd);			
			if (it1 != templBondTypeInfo.end()){
				firstTmp = it1->second;
			}
			firstTmp.push_back(ithBI.bt);
			templBondTypeInfo[ithBI.aInd] = firstTmp;

			vector<BondType> secondTmp;
			map<int, vector<BondType> >::iterator it2 = templBondTypeInfo.find(ithBI.bInd);			
			if (it2 != templBondTypeInfo.end()){
				secondTmp = it2->second;
			}
			secondTmp.push_back(ithBI.bt);
			templBondTypeInfo[ithBI.bInd] = secondTmp;

			// fill distance pow2 information into map
			vector<double> firstDis;
			map<int, vector<double> >::iterator dit1 = templBondDisInfo.find(ithBI.aInd);
			if (dit1 != templBondDisInfo.end()){
				firstDis = dit1->second;
			}
			firstDis.push_back(ithBI.bondDis);
			templBondDisInfo[ithBI.aInd] = firstDis;
			
			vector<double> secondDis;
			map<int, vector<double> >::iterator dit2 = templBondDisInfo.find(ithBI.bInd);
			if (dit2 != templBondDisInfo.end()){
				secondDis = dit2->second;
			}
			secondDis.push_back(ithBI.bondDis);
			templBondDisInfo[ithBI.bInd] = secondDis;
		}

		for (i = 0; i < query.size(); i++){
			vector<double> qRow;
			for (int j = 0; j < templ.size(); j++){
				qRow.push_back(-9999.0);
			}
			bondInfoScoMtx.push_back(qRow);
		}

		for (i = 0; i < query.size(); i++){
			map<int, vector<BondType> >::iterator qitBt = queryBondTypeInfo.find(i);	
			map<int, vector<double> >::iterator qitDis = queryBondDisInfo.find(i);
			if (qitBt == queryBondTypeInfo.end() || qitDis == queryBondDisInfo.end()){
				continue;
			}
			vector<BondType> queryBondTypes = qitBt->second;
			vector<double> queryBondDis = qitDis->second;
			
			for (int j = 0; j < templ.size(); j++){
				map<int, vector<BondType> >::iterator titBt = templBondTypeInfo.find(j);	
				map<int, vector<double> >::iterator titDis = templBondDisInfo.find(j);
				if (titBt == templBondTypeInfo.end() || titDis == templBondDisInfo.end()){
					continue;
				}
				vector<BondType> templBondTypes = titBt->second;
				vector<double> templBondDis = titDis->second;

				double sco = similar4TwoAtom8BondInfo(queryBondTypes, queryBondDis,
									templBondTypes, templBondDis);
				bondInfoScoMtx[i][j] = sco;
			}
		}
	}
	
	/***************************************************************************************
	 * This function is a sub function for "calculateBondInfoScoMtx" function
	 ***************************************************************************************/
 	double similar4TwoAtom8BondInfo(vector<BondType> aAtomBondTypes, vector<double> aAtomBondDis, 
			vector<BondType> bAtomBondTypes, vector<double> bAtomBondDis){
		vector<vector<double> > itemScoMtx;
		int i  = 0, j = 0;
		for (i = 0; i < aAtomBondTypes.size(); i++){
			BondType a_bt = aAtomBondTypes[i];
			double a_bdis = aAtomBondDis[i];

			vector<double> qRow;
			for (j = 0; j < bAtomBondTypes.size(); j++){
				BondType b_bt = bAtomBondTypes[j];
				double b_bdis = bAtomBondDis[j];
				
				double disDiff = abs(a_bdis - b_bdis);
				double ij_sco = 1.0 / (1.0 + pow(disDiff / 1.2, 2.0));
				if (a_bt == b_bt && a_bt != NOTCONNECTED){
					switch (a_bt){
					case SINGLE:
						ij_sco *= 1.2;
						break;
					case DOUBLE:
						ij_sco *= 1.5;
						break;
					case TRIPLE:
					case DUMMY:
						ij_sco *= 1.7;
						break;
					case AROMATIC:
						ij_sco *= 2.0;
						break;
					default:
						ij_sco *= 1.1;
					}
				}

				qRow.push_back(ij_sco);
			}
			itemScoMtx.push_back(qRow);
		}

		EnhancedGS egs(itemScoMtx);
		AlignInfo* p_Ali = egs.getPAliInfo();
		double aliSco = p_Ali->sco;
		
		double simSco = aliSco / (1 + abs(1.0 * ((int)aAtomBondTypes.size() - (int)bAtomBondTypes.size())));
		return simSco;
	}
};

#endif
