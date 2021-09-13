#include "file_io.h"

using namespace std;

//--------------------------------------------------
vector<string> separate_line(string line)
{
	string s="";
	vector<string> sv;
	for(string::iterator i=line.begin(); i!=line.end(); i++)
	{
		if(*i!=',')
		{
			s += *i;
		}
		else
		{
			sv.push_back(s);
			s="";
		}
	}
	sv.push_back(s);

	return sv;
}

//--------------------------------------------------------------
void make_directory(string name)
{
	string narrow_string(name);
	wstring wide_string = wstring(narrow_string.begin(), narrow_string.end());
	const wchar_t* name_wchar = wide_string.c_str();
	#if defined(_WIN32)
		_wmkdir(name_wchar);
	#else 
		mkdir(name.c_str(), 0700); 
	#endif
}