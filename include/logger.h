#pragma once
#include<string.h>
#include<iostream>
#include<fstream>
#include<queue>

class Logger{

private:

std::queue<std::string>logs;

public:

void log(std::string inputString, int priority);
bool isEmpty();
std::string getLog();

};
