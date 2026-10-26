#include"utils.h"

wp sign(wp in){
        if(in==0.0) return in;
        return in/fabs(in);
};

