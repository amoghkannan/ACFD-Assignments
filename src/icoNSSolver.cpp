#include"icoNSSolver.h"

icoNSSolver::icoNSSolver(){

        std::ifstream infile(meshFile);
        int imx, jmx;

        infile>>imx>>jmx;

        mesh = new Mesh(Box(),{NODE,NODE});

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        infile>>(*mesh)(i,j).x>>(*mesh)(i,j).y;
                };
        };

        mesh->setGhostNodes();
        mesh->calcAreas();
        mesh->calcVolumes();

        infile.close();

        std::ifstream infile1(settingsFile);

        infile1.close();

        (*this).setVar("u",{CELL,CELL});
        (*this).setVar("v",{CELL,CELL});
        (*this).setVar("p",{CELL,CELL});
        (*this).setVar("u*",{CELL,CELL});
        (*this).setVar("v*",{CELL,CELL});
        (*this).setVar("p*",{CELL,CELL});

};
