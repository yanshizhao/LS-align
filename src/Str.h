#ifndef STR_H
#define STR_H

#include <vector>
#include <string>
#include <iterator>
#include <cctype>  
#include <string>  
#include <algorithm>  
using namespace std;

string trim(const string &str)
{
	string s(str.c_str());
	if( !s.empty() )
    	{
        	s.erase(0,s.find_first_not_of(" "));
        	s.erase(s.find_last_not_of(" ") + 1);

    	}
		
	return s;
}

string eraseAll(const string &str, char ch){
    
    string s(str.c_str());
    int index = 0;
    if( !s.empty())
    {
        while( (index = s.find(ch, index)) != string::npos)
        {
            s.erase(index,1);
        }
    }
		
		return s;
}

string eraseAll(const string &str, const char* arr, int len){
    
    string s(str.c_str());
    for (int i = 0; i < len; i++){
		s = eraseAll(s, arr[i]);
	}
		
	return s;
}

/***********************************************
 * if there two continue ch word, we just split one time
 * e.g. str = "a--cdd-d", ch = '-', ans = {"a", "cdd", "d"}
 ***********************************************/
vector<string> split(const string& str, char ch){
	vector<string> ans;
	string tmp;
	for (int i = 0; i < str.size(); i++){
		if (str[i] == ch){
			if (0 != tmp.size()){
				string ttt(tmp.c_str());
				ans.push_back(ttt);
			}
			
			tmp.erase(tmp.begin(), tmp.end());
		}else{
			tmp += (str[i]);
		}
	}
	
	if (0 != tmp.size()){
		string ttt(tmp.c_str());
		ans.push_back(ttt);
	}
	
	return ans;
}

/***********************************************
 * if there two continue ch word, we just split one time
 * e.g. str = "a--cdd-d", ch = '-', ans = {"a", "cdd", "d"}
 ***********************************************/
vector<string> split(const string& str, char ch1, char ch2){
	vector<string> ans;
	string tmp;
	for (int i = 0; i < str.size(); i++){
		if (str[i] == ch1 || str[i] == ch2){
			if (0 != tmp.size()){
				string ttt(tmp.c_str());
				ans.push_back(ttt);
			}
			
			tmp.erase(tmp.begin(), tmp.end());
		}else{
			tmp += (str[i]);
		}
	}
	
	if (0 != tmp.size()){
		string ttt(tmp.c_str());
		ans.push_back(ttt);
	}
	
	return ans;
}

void toUpperString(string &str)  
{  
   transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper);  
}
  
void toLowerString(string &str)  
{  
   transform(str.begin(), str.end(), str.begin(), (int (*)(int))tolower);  
}  


#endif
