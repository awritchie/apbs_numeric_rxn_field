//
//  main.cpp
//  apbs_numeric_rxn_field
//
//  Created by Andrew Ritchie on 5/19/14.
//  Copyright (c) 2014 Andrew Ritchie. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>

int main(int argc, const char * argv[])
{
    const char *pqr, *solv, *ref;
    if (argc == 4 || argc == 5)
    {
        pqr = argv[1];
        solv = argv[2];
        ref = argv[3];
    }
    else
    {
        std::cerr << "\nUsage: <pqr> <apbs atom potentials in solvent> <apbs atom potentials in vacuum> <opt: fit to polynomial order N>\n" << std::endl;
        std::exit(1);
    }
    
    /* Read dummy atom coords from PQR file */
    std::string pqrline,throwAway,atom;
    std::vector<float> coords(3,0);
    std::vector<int> DumLineNumbers;
    int n = 0;
    std::vector<std::vector<float> > coordArray;
    std::ifstream file (pqr);
    if (file.is_open())
    {
        while (file.good())
        {
            getline(file, pqrline);
            if (not pqrline.empty())
            {
                std::stringstream linestream(pqrline);
                linestream >> throwAway >> throwAway >> atom >> throwAway >> throwAway >> coords[0] >> coords[1] >> coords[2];
                if (std::strncmp(atom.c_str(),"DUM",4) == 0)
                {
                    coordArray.push_back(coords);
                    DumLineNumbers.push_back(n);
                }
            }
            n++;
        }
    }
    int nDum = (int)coordArray.size();
    
    /* Read potentials from APBS atom potentials in solvent */
    std::string solvline;
    std::vector<float> solvpot(nDum,0);
    std::ifstream solvfile (solv);
    int i = 0, dumN = 0;
    if (solvfile.is_open())
    {
        while (solvfile.good())
        {
            getline(solvfile,solvline);
            if (not solvline.empty())
            {
                std::stringstream solvstream(solvline);
                solvstream >> throwAway;
                if (std::strncmp(throwAway.c_str(),"#",1) != 0)
                {
                    if (i == DumLineNumbers[dumN])
                    {
                        std::stringstream pot(solvline);
                        pot >> solvpot[dumN];
                        dumN++;
                    }
                    i++;
                }
            }
        }
    }
    
    /* Read potentials from APBS atom potentials in vacuum */
    std::string refline;
    std::vector<float> refpot(nDum,0);
    std::ifstream reffile (ref);
    i = 0;
    dumN = 0;
    if (reffile.is_open())
    {
        while (reffile.good())
        {
            getline(reffile,refline);
            if (not refline.empty())
            {
                std::stringstream refstream(refline);
                refstream >> throwAway;
                if (std::strncmp(throwAway.c_str(),"#",1) != 0)
                {
                    if (i == DumLineNumbers[dumN])
                    {
                        std::stringstream pot(refline);
                        pot >> refpot[dumN];
                        dumN++;
                    }
                    i++;
                }
            }
        }
    }
    
    /* Rxn potential is (potential in solvent) - (potential in vacuum) */
    std::vector<float> rxnpot(nDum,0);
    for (int i=0; i<nDum; i++)
    {
        rxnpot[i] = solvpot[i] - refpot[i];
    }
    
    /* 
       Assuming equally spaced dummy atoms, the middle 2 values are the closest to the
       midpoint, and will therefore have the smallest error in the first derivative. 
    */
    int n0 = DumLineNumbers[0], n1 = DumLineNumbers[nDum-1];
    while (n0+1 != n1 && n0 != n1-1 && n0+1 != n1-1 && n1+1 > n0-1)
    {
        n0++;
        n1--;
    }

    /* Field is -dV/dR */
    float dx2, dy2, dz2, dR, dV;
    dx2 = (coordArray[n0][0]-coordArray[n1][0])*(coordArray[n0][0]-coordArray[n1][0]);
    dy2 = (coordArray[n0][1]-coordArray[n1][1])*(coordArray[n0][1]-coordArray[n1][1]);
    dz2 = (coordArray[n0][2]-coordArray[n1][2])*(coordArray[n0][2]-coordArray[n1][2]);
    dR = sqrt(dx2+dy2+dz2);
    dV = rxnpot[n1] - rxnpot[n0];
    
    float num_der = -dV/dR;
    std::cout << num_der << " " << dR;
    
    if (argc == 5 )
    {
        using namespace Eigen;
        int nPoly;
        try
        {
            const char *fitOrder = argv[4];
            nPoly = atoi(fitOrder);
            if (nDum/2 < nPoly)
            {
                nPoly = nDum/2;
            }
        }
        catch(...)
        {
            return 0;
        }
        MatrixXd A(nDum,nPoly);
        VectorXd y(nDum);
        VectorXd coefficients(nPoly);
        for (int i=0; i<nDum; i++)
        {
            dx2 = (coordArray[i][0]-coordArray[0][0])*(coordArray[i][0]-coordArray[0][0]);
            dy2 = (coordArray[i][1]-coordArray[0][1])*(coordArray[i][1]-coordArray[0][1]);
            dz2 = (coordArray[i][2]-coordArray[0][2])*(coordArray[i][2]-coordArray[0][2]);
            dR = sqrt(dx2+dy2+dz2);
            for (int j=0; j<nPoly; j++)
            {
                A(i,j) = pow(dR,j);
            }
            y(i) = rxnpot[i];
        }

        coefficients = (A.transpose()*A).inverse()*A.transpose()*y;
        //std::cout << "PolyFit: ";
        //std::cout << coefficients(0) << "*x**" << i;
        float deriv = 0;
        float sderiv = 0;
        for (int i=0;i<nPoly;i++)
        {
            deriv += coefficients(i)*i*pow(dR/2.,i-1);
            sderiv += coefficients(i)*i*(i-1)*pow(dR/2.,i-2);
            //std::cout << coefficients(i) << "*x**" << i;
            //if (i < nPoly-1) { std::cout << " + ";}
        }
        int npoints = 51;
        float f_at_x[npoints];
        for (int j=0;j<npoints;j++)
        {
            f_at_x[j] = 0;
            for (int i=0; i<nPoly; i++)
            {
                if (not isnan(i*pow(j/(float)(npoints-1)*dR,i-1)))
                {
                    f_at_x[j] += coefficients(i)*i*pow(j/(float)(npoints-1)*dR,i-1);
                }
            }
        }
        std::cout << " " << -deriv << " " << -sderiv << " " << nPoly << " " << -deriv - num_der << std::endl;
        float max_abs_f = 0;
        float max_abs_f_r = 0;
        for (int i=0;i<npoints;i++)
        {
            std::cout << i/(float)(npoints -1)* dR << " " << -f_at_x[i] << "\n";
            if (fabs(f_at_x[i]) > fabs(max_abs_f))
            {
                max_abs_f = f_at_x[i];
                max_abs_f_r = i/(float)(npoints-1)*dR;
            }
        }
        std::cout << "Field max at " << max_abs_f_r << " of " << -max_abs_f << " KbT/eA\n";
    }
    else
    {
        std::cout << "\n";
    }
    return 0;
}

