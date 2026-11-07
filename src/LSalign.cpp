#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include "Point3D.h"
#include "Ligands.h"
#include "basic_fun.h"
#include "Str.h"
#include "RigidLSalign.h"
#include "AtomFixInfo.h"
#include "Sort.h"
#include "EnhancedGS.h"
#include <iterator>
#include "Kabsch.h"
#include "ShortestRoad.h"
#include "LigRotamerDetector.h"
#include "FlexibleLSalign.h"
#include "FragmentLSalign.h"
#include "time.h"
#include "DijkstraSRoad.h"
using namespace std;

char version[100];

void print_help(char* arg){
	cout <<endl;
	cout << " *********************************************************************************************************************" << endl
		 << " * LS-align (Version "<<version<<"): An atom-level ligand structural alignment algorithm                             *" << endl
		 << " * Based on statistics:                                                                                              *" << endl
		 << " *          0.00 < PC-score < 0.46, random structural similarity                                                     *" << endl
		 << " *          0.46 < PC-score < 1.00, in about the same sub-structure                                                  *" << endl
		 << " * Reference: J. Hu, Z. Liu, D.J. Yu and Y. Zhang, LS-align: an atom-level, flexible ligand structural alignment     *" << endl
		 << " *            algorithm for efficient virtual screening, Bioinformatics, submitted.                                  *" << endl
		 << " *                                                                                                                   *" << endl
		 << " * Please email your comments and suggestions to Jun Hu (junh_cs@126.com) or Yang Zhang (zhng@umich.edu)             *" << endl
		 << " *********************************************************************************************************************" << endl;	
	cout << endl
		 << " Usage: " << arg << " QUERYs.mol2 TEMPLs.mol2 [Options]" << endl << endl
		 << " Options:" << endl
		 << "       -rf   use Rigid-LS-align of Flexi-LS-align algorithm" << endl
		 << "                 0 (default): use Rigid-LS-align " << endl
		 << "                 1          : use Flexi-LS-align " << endl
		 << "       -md   Only if \"-rf\" is be set 1, it is useful. It will control whether using more detail" << endl
		 << "             rotated angles of single bond to generate more rotamers of query ligand(s) or not." << endl
		 << "                 0 (default): use Raw rotated angle vector"<< endl
		 << "                 1          : use Detail rotated angle vector" << endl
		 << "       -d0   Set the d0 value of PC-score (suggestted value is 1.0 or 2.0 angstrom)" << endl
		 << "             warning: its value range from [0.1, 4] angstrom. " << endl
		 << "             default: d0 = 0.55 * pow(Nmin - 9, 1.0/3.0) + 0.15. " << endl
		 << "       -oq   If the value is setted, LS-align will only align the oq-th query ligand in QUERYs.mol2" << endl
		 << "             to all or ot-th (if \"-ot\" is setted) template ligand(s) in TEMPLS.mol2" << endl    
		 << "             NOTICE: If the value larger than the number of ligands in QUERYs.mol2 or less than 1,"<<endl
		 << "                     we set its value as 1" << endl
		 << "       -ot   If the value is setted, LS-align will only align the all or oq-th query (if \"-oq\" is setted)"<<endl
		 << "             ligand(s) in QUERYs.mol2 to the ot-th template ligand in TEMPLS.mol2" << endl
		 << "             NOTICE: If the value larger than the number of ligands in TEMPLs.mol2 or less than 1,"<<endl
		 << "                     we set its value as 1" << endl
		 << "       -o    output the first rotated query structure which can superpose the first templ structure." << endl
         << "             if the \"-oq\" are \"-ot\" are both not setted." <<endl << endl

		 << "             output the oq-th rotated query structure which can superpose the first templ structure" << endl
         << "             if the \"-oq\" is setted and \"-ot\" is not setted." <<endl << endl

		 << "             output the first rotated query structure which can superpose the ot-th templ structure" << endl
         << "             if the \"-oq\" is not setted and \"-ot\" is setted." <<endl << endl

		 << "             Or, output the the oq-th rotated query structure which can superpose the ot-th templ structure,"<< endl
		 << "             if the \"-oq\" and \"-ot\" are both setted."<< endl << endl

		 << "       -a    output the alignment information between first query structure with the first templ structure" << endl
		 << "             if the \"-oq\" are \"-ot\" are both not setted." <<endl << endl

		 << "             output the alignment information between the oq-th query structure with the first templ structure." << endl
         << "             if the \"-oq\" is setted and \"-ot\" is not setted." <<endl << endl

		 << "             output the alignment information between the first query structure with the ot-th templ structure." << endl
         << "             if the \"-oq\" is not setted and \"-ot\" is setted." <<endl<<endl

		 << "             Or, output the alignment information between the oq-th query structure with the ot-th templ structure,"<< endl
		 << "             if the \"-oq\" and \"-ot\" are both setted."<< endl << endl

		 << "       -m    output the rotation matrix which can rotate the first query structure to first templ structure" << endl
		 << "             if the \"-oq\" are \"-ot\" are both not setted." <<endl<< endl

		 << "             output the rotation matrix between the oq-th query structure with the first templ structure." << endl
         << "             if the \"-oq\" is setted and \"-ot\" is not setted." <<endl<<endl

		 << "             output the rotation matrix between the first query structure with the ot-th templ structure." << endl
         << "             if the \"-oq\" is not setted and \"-ot\" is setted." <<endl<<endl
		
		 << "             Or, output the rotation matrix between the oq-th query structure with the ot-th templ structure,"<< endl
		 << "             if the \"-oq\" and \"-ot\" are both setted."<< endl<<endl

		 << "       -acc  Do you want to take some time to search more accuracy alignment. "<<endl
		 << "                 0 (default) : You do not."<<endl
		 << "                 1           : You do."<<endl<<endl

		 << "       -H    Do you want to load the hydrogen atom(s) in the mol2 file. "<<endl
		 << "                 0 (default) : You do not. (It will run faster!)"<<endl
		 << "                 1           : You do. (It will run slower!)"<<endl<<endl

		 << "       -v    print the version of LS-align" << endl << endl
		 << "       -h    print this help" << endl << endl
		 << " Example usages:" << endl
		 << "        "<< arg <<" QUERYs.mol2 TEMPLs.mol2" << endl
		 << "        "<< arg <<" QUERYs.mol2 TEMPLs.mol2 -d0 3.0" << endl
		 << "        "<< arg <<" QUERYs.mol2 TEMPLs.mol2 -rf 1 -d0 3.0" << endl
		 << "        "<< arg <<" QUERYs.mol2 TEMPLs.mol2 -oq 2 -ot 3 -o QUERYs.sup" << endl
		 << "        "<< arg <<" QUERYs.mol2 TEMPLs.mol2 -a align.txt" << endl
	     << "        "<< arg <<" QUERYs.mol2 TEMPLs.mol2 -m matrix.txt" << endl
		 << "        "<< arg <<" -h"<< endl
		 << "        "<< arg <<" -v" << endl << endl
		 << " NOTICE: " << endl
		 << "         The input mol2 file can contain lots of ligands with mol2 format." << endl << endl;

	cout << "LICENCE: FREE FOR ACADEMIC USE. IF YOU WANT, YOU CAN CHANGE THIS CODE." << endl << endl;
	exit(EXIT_SUCCESS);
}

// ---------------------------------------------------------------------------//
// ResultTable: single definition of every result column.  Header and data rows
// are produced here so a column change cannot drift apart.  The header is
// printed lazily on the first row (or by finish() when there is no row) so data
// already computed stays visible if a run dies halfway.
// ---------------------------------------------------------------------------//
class ResultTable{
public:
	struct Col{
		const char* h1;		// header, first line
		const char* h2;		// header, second line ("" = blank field)
		int width;		// field width shared by headers and data
		const char* dfmt;	// data conversion; its width must equal 'width'
	};
	static const Col kCols[10];

	// Right-align a ligand name in 18 characters (the original inline loops).
	static void formatName18(const string& name, char out[19]);

	ResultTable() : m_headerPrinted(false) {}

	void emitRow(const char* q18, const char* t18,
	             double pcQ, double pcT, double pvalQ, double pvalT,
	             double jaccard, double rmsd, int nq, int nt, bool isLoad_H);

	void finish(clock_t start, clock_t middle, clock_t stop, bool isLoad_H);

private:
	void printHeader(bool isLoad_H);
	bool m_headerPrinted;
};

const ResultTable::Col ResultTable::kCols[10] = {
	{"QEURY_NAME",   "", 18, "%18s"},
	{"TEMPL_NAME",   "", 18, "%18s"},
	{"PC-score8Q",   "", 10, "%10.4f"},
	{"PC-score8T",   "", 10, "%10.4f"},
	{"Pval_PC8Q",    "", 10, "%10.6f"},
	{"Pval_PC8T",    "", 10, "%10.6f"},
	{"JaccardR",     "", 10, "%10.4f"},
	{"RMSD_lb",      "", 10, "%10.4f"},
	{"  QUERY_SIZE ", "", 14, "%14d"},
	{"  TEMPL_SIZE ", "", 14, "%14d"}
};

void ResultTable::formatName18(const string& name, char out[19]){
	int i;
	for (i = 0; i < 18; i++){
		if (i < (int)name.size())
			out[17-i] = name[(int)name.size()-1-i];
		else
			out[17-i] = ' ';
	}
	out[18] = '\0';
}

void ResultTable::printHeader(bool isLoad_H){
	char row[1000];
	char fmt[256];
	int len = 0;
	int i;

	cout << "" << endl;
	for (i = 0; i < 10; i++){
		if (i > 0)
			fmt[len++] = ' ';
		len += sprintf(fmt + len, "%%%ds", kCols[i].width);
	}
	fmt[len] = '\0';
	sprintf(row, fmt,
	        kCols[0].h1, kCols[1].h1, kCols[2].h1, kCols[3].h1, kCols[4].h1,
	        kCols[5].h1, kCols[6].h1, kCols[7].h1, kCols[8].h1, kCols[9].h1);
	cout << row << endl;

	const char* h2Last = isLoad_H ? "(With \"H\")" : "(Without \"H\")";
	sprintf(row, fmt,
	        kCols[0].h2, kCols[1].h2, kCols[2].h2, kCols[3].h2, kCols[4].h2,
	        kCols[5].h2, kCols[6].h2, kCols[7].h2, h2Last, h2Last);
	cout << row << endl;
}

void ResultTable::emitRow(const char* q18, const char* t18,
                          double pcQ, double pcT, double pvalQ, double pvalT,
                          double jaccard, double rmsd, int nq, int nt, bool isLoad_H){
	if (!m_headerPrinted){
		printHeader(isLoad_H);
		m_headerPrinted = true;
	}

	char row[1000];
	char fmt[256];
	int len = 0;
	int i;
	for (i = 0; i < 10; i++){
		if (i > 0)
			fmt[len++] = ' ';
		len += sprintf(fmt + len, "%s", kCols[i].dfmt);
	}
	fmt[len] = '\0';
	sprintf(row, fmt, q18, t18, pcQ, pcT, pvalQ, pvalT, jaccard, rmsd, nq, nt);
	cout << row << endl;
}

void ResultTable::finish(clock_t start, clock_t middle, clock_t stop, bool isLoad_H){
	if (!m_headerPrinted){
		printHeader(isLoad_H);
		m_headerPrinted = true;
	}
	cout << endl;
	cout << "Using " << (1.0*(middle-start)/CLOCKS_PER_SEC) << "-s to load mol2 files." << endl;
	cout << "Using " << (1.0*(stop-middle)/CLOCKS_PER_SEC) << "-s to do LS-align. " << endl;
	cout << "Totally Using " << (1.0*(stop-start)/CLOCKS_PER_SEC) << "-s." << endl << endl;
}

int main(int argc, char* args[]){
	strcpy(version, "J201704171741");
	clock_t start, middle, finish;
	start = clock();

	if (argc < 2){
		print_help(args[0]);
	}

	string queryMol2s;
	string templMol2s;
	bool isQuerysProvide = false;
	bool isTemplsProvide = false;
	
	int is_Rigid_Flexi_Fragm = 0;	// default using Rigid-LS-align, 0 mean Rigid, 1 : Flexible, 2 : Fragment

	bool rawOrDetailSBRotAngles = false; // false mean raw angles, true mean detail angles

	double d0 = 2.0;
	bool isDefult_d0 = true;

	bool isOnlyRunQthLigSetted = false;
	int OnlyRunQthLig = 0;
	bool isOnlyRunTthLigSetted = false;
	int OnlyRunTthLig = 0;
	bool isSaveSupPath = false;
	string saveSupPath;	

	bool isSaveAlignInfo = false;
	string saveAlignPath;

	bool isSaveRotMtx = false;
	string saveRotMtxPath;

	bool isNeedMoreCorrectAlign = false; // TODO

	bool isLoad_H = false;

	int i = 0, j = 0, mol2Ind = 0;
	// --------------------------------------------------------------------------//
	// ----------------          Load Input Parameters      ---------------------//
	// --------------------------------------------------------------------------//
	for (i = 1; i < argc; i++){
		if (0 == strcmp(args[i], "-rf") && i < argc-1){
			int rf = atoi(args[++i]);
			if (rf == 0) is_Rigid_Flexi_Fragm = 0;
			else if (rf == 1) is_Rigid_Flexi_Fragm = 1;
			else is_Rigid_Flexi_Fragm = 2;
		}else if (0 == strcmp(args[i], "-md") && i < argc-1){
			int md = atof(args[++i]);
			if (md != 0) rawOrDetailSBRotAngles = true;
		}else if (0 == strcmp(args[i], "-d0") && i < argc-1){
			d0 = atof(args[++i]);
			if (d0 < 0.1) d0 = 0.1;
			else if (d0 > 4) d0 = 4;
			isDefult_d0 = false;
		}else if (0 == strcmp(args[i], "-oq") && i < argc-1){
			isOnlyRunQthLigSetted = true;
			OnlyRunQthLig = atoi(args[++i]) - 1;
		}else if (0 == strcmp(args[i], "-ot") && i < argc-1){
			isOnlyRunTthLigSetted = true;
			OnlyRunTthLig = atoi(args[++i]) - 1;
		}else if (0 == strcmp(args[i], "-o") && i < argc-1){
			isSaveSupPath = true;
			saveSupPath = args[++i];
		}else if (0 == strcmp(args[i], "-a") && i < argc-1){
			isSaveAlignInfo = true;
			saveAlignPath = args[++i];
		}else if (0 == strcmp(args[i], "-m") && i < argc-1){
			isSaveRotMtx = true;
			saveRotMtxPath = args[++i];
		}else if (0 == strcmp(args[i], "-acc") && i < argc-1){
			int acc = atof(args[++i]);
			if (acc != 0) isNeedMoreCorrectAlign = true;
		}else if (0 == strcmp(args[i], "-H") && i < argc-1){
			int H = atof(args[++i]);
			if (H != 0) isLoad_H = true;
		}else if (0 == strcmp(args[i], "-v")){
			cout <<endl;
			cout << " *********************************************************************************************************************" << endl
				 << " * LS-align (Version "<<version<<"): An atom-level ligand structural alignment algorithm                             *" << endl
				 << " * Based on statistics:                                                                                              *" << endl
				 << " *          0.00 < PC-score < 0.46, random structural similarity                                                     *" << endl
				 << " *          0.46 < PC-score < 1.00, in about the same sub-structure                                                  *" << endl
				 << " * Reference: J. Hu, Z. Liu, D.J. Yu and Y. Zhang, LS-align: an atom-level, flexible ligand structural alignment     *" << endl
				 << " *            algorithm for efficient virtual screening, Bioinformatics, submitted.                                  *" << endl
				 << " *                                                                                                                   *" << endl
				 << " * Please email your comments and suggestions to Jun Hu (junh_cs@126.com) or Yang Zhang (zhng@umich.edu)             *" << endl
				 << " *********************************************************************************************************************" << endl;
			exit(EXIT_FAILURE);
		}else if (0 == strcmp(args[i], "-h")){
			print_help(args[0]);
		}else{
			if (mol2Ind < 2){
				fstream _file;
				_file.open(args[i], ios::in);
				if(!_file){
					cout << args[i] << " is not existed!" << endl;
					cout << "Please check it and input a correct mol2 file path!" << endl;
					exit(EXIT_FAILURE);
				}
			}			

			if (mol2Ind == 0){
				isQuerysProvide = true;
				queryMol2s = args[i];      
			}else if (mol2Ind == 1){
				isTemplsProvide = true;
				templMol2s = args[i];      
			}
			mol2Ind++;
		}
	}

	if (!isQuerysProvide){
		cout << "PLEASE PROVIDE QUERYs.mol2 FILE!!!" << endl << endl;
		print_help(args[0]);
	}
	if (!isTemplsProvide){
		cout << "PLEASE PROVIDE TEMPLs.mol2 FILE!!!" << endl << endl;
		print_help(args[0]);
	}
	
	// --------------------------------------------------------------------------//
	// ----------------             Load Mol2 Ligands       ---------------------//
	// --------------------------------------------------------------------------//
	// load query and templ mol2 ligands 
	Ligands queryligs(queryMol2s, isLoad_H);
	Ligands templligs(templMol2s, isLoad_H);

//	// TESE DELETE START
//	is_Rigid_Flexi_Fragm = 1;
//	isOnlyRunQthLigSetted = true;
//	OnlyRunQthLig = 3;
//	//Ligands queryligs("D:/delete/LSalign/query.mol2", isLoad_H);
//	//Ligands templligs("D:/delete/LSalign/templ.mol2", isLoad_H);
//	Ligands queryligs("D:/delete/LSalign/xtal-lig.mol2", isLoad_H);
//	Ligands templligs("D:/delete/LSalign/ace_ligands.mol2", isLoad_H);
//	// TESTE DELETE END

	if (queryligs.size() <= 0){
		cout << "PLEASE PROVIDE QUERYs.mol2 FILE!!!" << endl << endl;
		print_help(args[0]);
	}
	if (templligs.size() <= 0){
		cout << "PLEASE PROVIDE TEMPLs.mol2 FILE!!!" << endl << endl;
		print_help(args[0]);
	}

	if (isOnlyRunQthLigSetted && (queryligs.size() <= OnlyRunQthLig || OnlyRunQthLig < 0))
		OnlyRunQthLig = 0;
	if (isOnlyRunTthLigSetted && (templligs.size() <= OnlyRunTthLig || OnlyRunTthLig < 0))
		OnlyRunTthLig = 0;

	middle = clock();

	// --------------------------------------------------------------------------//
	// ----------------             Run LS-Align Algorithm    ------------------//
	// --------------------------------------------------------------------------//
	// Using LS-align to superpose Querys and Templs
	ResultTable table;

	for (i = 0; i < queryligs.size(); i++){
		if (isOnlyRunQthLigSetted && i != OnlyRunQthLig){
			continue;
		}

		for (j = 0; j < templligs.size(); j++){
			if (isOnlyRunTthLigSetted && j != OnlyRunTthLig){
				continue;
			}

			char ithQnameChArr[19];
			char jthTnameChArr[19];
			ResultTable::formatName18(queryligs[i].getLigname(), ithQnameChArr);
			ResultTable::formatName18(templligs[j].getLigname(), jthTnameChArr);
			
			if (isDefult_d0){
				//----------------------------------------------------------------------------//
				//-------  Set the default d0 value based on the following equation     ------//
				//-------  d0 = 0.55 * pow(Nmin - 9, 1.0/3.0) + 0.15;     ------//
				//----------------------------------------------------------------------------//
				int Nq = queryligs[i].size();
				int Nt = templligs[j].size();
				
				int Nmin = Nq<Nt ? Nq : Nt;
				
				if (Nmin < 9){
					d0 = 0.1;
				}else{
					// d0 = 0.55 * pow(Nmin - 9, 1.0/3.0) + 0.15;
					d0 = 0.55 * pow(Nmin - 9, 1.0/3.0) + 0.15;
				}
				
				if (d0 < 0.1)
					d0 = 0.1;
				if (d0 > 4.0)
					d0 = 4.0;
			}

			//----------------------------------------------------------------------------//
			//-------  Calcualte the u and delta for obtaining P-value              ------//
			//----------------------------------------------------------------------------//
			double Nq = 1.0*queryligs[i].size();
			double Nt = 1.0*templligs[j].size();
			double u_T     = -0.563358 +  0.152792 * log(Nq) +  0.147930 *log(Nt);
			double delta_T =  0.110170 + -0.007444 * log(Nq) + -0.008468 *log(Nt);
			if (delta_T == 0.0) delta_T = 0.00000000000000001;

			double u_Q     = -0.563358 +  0.152792 * log(Nt) +  0.147930 *log(Nq);
			double delta_Q =  0.110170 + -0.007444 * log(Nt) + -0.008468 *log(Nq);
			if (delta_Q == 0.0) delta_Q = 0.00000000000000001;
		

			if (0 == is_Rigid_Flexi_Fragm || 1 == is_Rigid_Flexi_Fragm){
				RigidLSalign* p_ali = NULL;
				if (0 == is_Rigid_Flexi_Fragm){	// Flexible LS-align
					p_ali = new RigidLSalign(queryligs[i], templligs[j], d0);
				}else{
					FlexibleLSalign* p_fLSalign = new FlexibleLSalign(queryligs[i], templligs[j], d0, 
						rawOrDetailSBRotAngles, true, isNeedMoreCorrectAlign);
					p_ali = p_fLSalign->align(10);

					delete p_fLSalign;
					p_fLSalign = NULL;
				}
//				double LSscore = p_ali->getLSscore();
				double LSscore8QuerySize = p_ali->getLSscore8QuerySize();
				double LSscore8TemplSize = p_ali->getLSscore8TemplSize();
				double Pval_8Q = 1.0 - exp( -exp(- ((LSscore8QuerySize - u_Q)/delta_Q)));
				double Pval_8T = 1.0 - exp( -exp(- ((LSscore8TemplSize - u_T)/delta_T)));
				double JaccardR = p_ali->getJaccardRatio();
//				double alisco = p_ali->getAliSco();
				double RMSD_lb_val = p_ali->getRMSD_lb();
				double RMSD_lb_out = (RMSD_lb_val < 0.0) ? std::numeric_limits<double>::quiet_NaN() : RMSD_lb_val;
//				double RMSD_lb_mim = p_ali->getRMSD_lb_min();
				
				table.emitRow(ithQnameChArr, jthTnameChArr,
					LSscore8QuerySize, LSscore8TemplSize, Pval_8Q, Pval_8T, JaccardR, RMSD_lb_out, queryligs[i].size(), templligs[j].size(), isLoad_H);
				
				// save the superposition between the Qth-query ligand with the Tth-templ ligand
				if (isSaveSupPath 
					&& i == OnlyRunQthLig
					&& j == OnlyRunTthLig){
						Ligand rotQuery = p_ali->getRottedQueryLigand();
						vector<string> rotQVstr = rotQuery.format2PDB("QUE");

						vector<string> TVstr = templligs[j].format2PDB("TEM", 1001, 99999);

						ofstream fout(saveSupPath.c_str());	
						int k = 0;
						for (k = 0; k < rotQVstr.size(); k++){
							fout << rotQVstr[k] << endl;
						}
						for (k = 0; k < TVstr.size(); k++){
							fout << TVstr[k] << endl;
						}
						fout.close();
				} // end if (isSaveSupPath 

				// save the alignment information between the Qth-query ligand with the Tth-templ ligand
				if (isSaveAlignInfo
					&& i == OnlyRunQthLig
					&& j == OnlyRunTthLig){
					vector<string> aliInfos = p_ali->formatAliInfo();

					ofstream fout(saveAlignPath.c_str());	
					int k = 0;
					for (k = 0; k < aliInfos.size(); k++){
						fout << aliInfos[k] << endl;
					}
					fout.close();	
				}// end if (isSaveAlignInfo

				if (isSaveRotMtx
					&& i == OnlyRunQthLig
					&& j == OnlyRunTthLig){
					vector<string> rotMtxInfo = p_ali->formatRotMtx();

					ofstream fout(saveRotMtxPath.c_str());	
					int k = 0;
					for (k = 0; k < rotMtxInfo.size(); k++){
						fout << rotMtxInfo[k] << endl;
					}
					fout.close();	
				}// end if (isSaveRotMtx

				delete p_ali;
				p_ali = NULL;
			}
		}// end for (j = 0; j < templligs.size(); j++){
	}// end for (i = 0; i < queryligs.size(); i++){

	finish = clock();
	table.finish(start, middle, finish, isLoad_H);

	return 1;
}
