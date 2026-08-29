#include"logger.h"

void Logger::log(std::string inputstring, int priority){
        logs.push(inputstring);

        if(priority==1){
                while(!logs.empty()){
                        std::string output=logs.front();
                        std::cout<<output<<std::endl;
                        logs.pop();
                };
        };

};

std::string Logger::getLog(){
        std::string output=logs.front();
        logs.pop();
        return output;
};

bool Logger::isEmpty(){
        return logs.empty();
};
