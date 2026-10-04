#include"vars.h"

//Since this file defines a template class, function definitions MUST be in the .h file itself, this .cpp is left empty

wp findMatElement(int iCell,int jCell, int iEntry, int jEntry, Grid<boundMatRow>&mat){

        std::vector<boundMatEntry> dependencies;
        Point p;
        wp data,ans;

        dependencies=mat(iCell,jCell);

        ans=0.0;

        for(boundMatEntry item:dependencies){
               p=item.first; 
               data=item.second;
               
               if(p.first==iEntry && p.second==jEntry){
                         ans=data;
                         break;
               };
        };

        return ans;

};

void changeMatElement(int iCell,int jCell, int iEntry, int jEntry, wp val, Grid<boundMatRow>&mat){

        std::vector<boundMatEntry> dependencies;
        Point p;
        wp data;

        dependencies=mat(iCell,jCell);

        bool isFound=false;
        int counter=0;

        for(boundMatEntry item:dependencies){
               p=item.first; 
               data=item.second;
               
               if(p.first==iEntry && p.second==jEntry){
                         isFound=true;
                         mat(iCell,jCell)[counter].second+=val;
                         break;
               };

               counter=counter+1;
        };

        if(!isFound){
                boundMatEntry newEntry={{iEntry,jEntry},val};
                mat(iCell,jCell).push_back(newEntry);
        };
};

