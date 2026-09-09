#include"mesh.h"

Mesh::~Mesh(){
delete[] data;
logger.log("Debug: destroying mesh",1);
};
