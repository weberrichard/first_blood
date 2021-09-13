#ifndef FILE_IO_H
#define FILE_IO_H

#include <string>
#include <vector>
#include <sys/stat.h> // mkdir

using namespace std;

// separating string to vector<string>, comma
vector<string> separate_line(string line);

// make new directory, works for windows and linux
void make_directory(string name);

#endif