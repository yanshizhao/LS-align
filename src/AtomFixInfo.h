#ifndef ATOMFIXINFO_H
#define ATOMFIXINFO_H

#include <string.h>
#include <map>
#include "Str.h"
using namespace std;

class AtomVarDerWaalRadius{

private :
	map<string, double> vdwRMap;
	

public:
	AtomVarDerWaalRadius(){
		
		// the var der waal radius information is come from "Van der Waals Radii of Elements" (S. S. Batsanov)
		string varDerWaalRadiusStr = "H   1.2 \nHe 1.4 \nLi 2.2 \nBe 1.9 \nB 1.8 \nC 1.7 \nN 1.6 \nO 1.55 \nF 1.5 \nNe 1.5 \nNa 2.4 \nMg 2.2 \nAl 2.1 \nSi 2.1 \nP 1.95 \nS 1.8 \nCl 1.8 \nAr 1.9 \nK 2.8 \nCa 2.4 \nSc 2.3 \nTi 2.15 \nV 2.05 \nCr 2.05 \nMn 2.05 \nFe 2.05 \nCo 2.0 \nNi 2.0 \nCu 2.0 \nZn 2.1 \nGa 2.1 \nGe 2.1 \nAs 2.05 \nSe 1.9 \nBr 1.9 \nKr 2.0 \nRb 2.9 \nSr 2.55 \nY 2.4 \nZr 2.3 \nNb 2.152 \nMo 2.1 \nTc 2.05 \nRu 2.05 \nRh 2.0 \nPd 2.05 \nAg 2.1 \nCd 2.2 \nIn 2.2 \nSn 2.25 \nSb 2.2 \nTe 2.1 \nI 2.1 \nXe 2.2 \nCs 3.0 \nBa 2.7 \nLa 2.5 \nCe 2.35 \nPr 2.39 \nNd 2.29 \nPm 2.36 \nSm 2.29 \nEu 2.33 \nGd 2.37 \nTb 2.21 \nDy 2.29 \nHo 2.16 \nEr 2.35 \nTm 2.27 \nYb 2.42 \nLu 2.21 \nHf 2.25 \nTa 2.2 \nW 2.1 \nRe 2.05 \nOs 2.0 \nIr 2.0 \nPt 2.05 \nAu 2.1 \nHg 2.05 \nTl 2.2 \nPb 2.3 \nBi 2.3 \nPo 2.29 \nAt 2.36 \nRn 2.43 \nFr 2.56 \nRe 2.43 \nRa 2.43 \nAc 2.6 \nTh 2.4 \nPa 2.43 \nU 2.3 \nNp 2.21 \nPu 2.56 \nAm 2.56 \nCm 2.56 \nBk 2.56 \nCf 2.56 \nEs 2.56 \nFm 2.56 \nMd \nNo \nLr \nRf \nDb \nSg \nBh \nHs \nMt \nDs \nRg \nCn \nUut \nUuq \nUup \nUuh \nUus \nUuo ";

		vector<string> atomInfos = split(varDerWaalRadiusStr, '\n');
		for (int i = 0; i < atomInfos.size(); i++){
			vector<string> atomRadii = split(atomInfos[i], ' ');
			if (2 == atomRadii.size()){
				string key(atomRadii[0]);
				toUpperString(key);
				
				double value = 0.0;
				sscanf(atomRadii[1].c_str(), "%lf", &value);
				
				vdwRMap[key] = value;
			}
		}
	}
	
	double operator[](string atomtype){
		string at = eraseAll(atomtype, ' ');
		toUpperString(at);

		return vdwRMap[at];
	}
};

class AtomRelativeMass{
private:
	map<string, double> relativeMassMap;
public:
	AtomRelativeMass(){
		string relativeMassInfo = "H  1.007947 \nHe 4.0026022 \nLi 6.9412 \nBe 9.0121823 \nB 10.8117 \nC 12.01078 \nN 14.00672 \nO 15.99943 \nF 18.99840325 \nNe 20.17976 \nNa 22.989769282\nMg 24.30506 \nAl 26.98153868 \nSi 28.08553 \nP 30.9737622 \nS 32.0655 \nCl 35.4532 \nAr 39.9481 \nK 39.09831 \nCa 40.0784 \nSc 44.9559126 \nTi 47.8671 \nV 50.94151 \nCr 51.99616 \nMn 54.9380455 \nFe 55.8452 \nCo 58.9331955 \nNi 58.69342 \nCu 63.5463 \nZn 65.4094 \nGa 69.7231 \nGe 72.641 \nAs 74.921602 \nSe 78.963 \nBr 79.9041 \nKr 83.7982 \nRb 85.46783 \nSr 87.621 \nY 88.905852 \nZr 91.2242 \nNb 92.906382 \nMo 95.942 \nTc 97.9072 \nRu 101.072 \nRh 102.905502 \nPd 106.421 \nAg 107.86822 \nCd 112.4118 \nIn 114.8183 \nSn 118.7107 \nSb 121.7601 \nTe 127.603 \nI 126.904473 \nXe 131.2936 \nCs 132.90545192\nBa 137.3277\nLa 138.905477 \nCe 140.1161 \nPr 140.907652 \nNd 144.2423 \nPm 145 \nSm 150.362 \nEu 151.9641 \nGd 157.253 \nTb 158.925352 \nDy 162.5001 \nHo 164.930322 \nEr 167.2593 \nTm 168.934212 \nYb 173.043 \nLu 174.9671 \nHf 178.492 \nTa 180.947882 \nW 183.841 \nRe 186.2071 \nRa 186.2071 \nOs 190.233 \nIr 192.2173 \nPt 195.0849 \nAu 196.9665694\nHg 200.592 \nTl 204.38332 \nPb 207.21 \nBi 208.980401 \nPo 208.9824 \nAt 209.9871 \nRn 222.0176 \nFr 223 \nRe 226 \nAc 227 \nTh 232.038062 \nPa 231.035882 \nU 238.028913 \nNp 237 \nPu 244 \nAm 243 \nCm 247 \nBk 247 \nCf 251 \nEs 252 \nFm 257 \nMd 258 \nNo 259 \nLr 262 \nRf 261 \nDb 262 \nSg 266 \nBh 264 \nHs 277 \nMt 268 \nDs 281 \nRg 272 \nCn 285 \nUut 284 \nUuq 289 \nUup 288 \nUuh 292 \nUus 291 \nUuo 293 ";

		vector<string> atomInfos = split(relativeMassInfo, '\n');
		for (int i = 0; i < atomInfos.size(); i++){
			vector<string> atomRadii = split(atomInfos[i], ' ');
			if (2 == atomRadii.size()){
				string key(atomRadii[0]);
				toUpperString(key);
				
				double value = 0.0;
				sscanf(atomRadii[1].c_str(), "%lf", &value);
				
				relativeMassMap[key] = value;
			}
		}
	}

	double operator[](string atomtype){
		string at = eraseAll(atomtype, ' ');
		toUpperString(at);

		return relativeMassMap[at];
	}
};




#endif
