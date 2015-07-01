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

#ifndef M_PI
#define M_PI 4*atan(1)
#endif

#define eps0 8.8541878e-12 //Vacuum permittivity C^2 s^2 /(kg m^3)
#define ec 1.6021773e-19 // Coulombs/e-
#define NA 6.0221367e23 // Avogadro's number, molecules/mol
#define kb 1.3806581e-23 // J/K
#define scale (ec/(4.*M_PI*eps0*1.E-10)) // kg m^2 A / (e- C s^2)
#define zmagic (ec*NA*1.E-3) // (C kJ/(E- mol J)
#define kbt (kb*300.*NA*1.E-3) // kJ/mol/kbT(300)
#define cfac (scale*zmagic/kbt) // kbT/(e- A)


float coulomb_field(const char *pqr, float *rArray, int nArray, float *fArray);

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
    while (n0+1 != n1-1 && n1+1 > n0-1)
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
        float deriv = 0, sderiv = 0;
        for (int i=0;i<nPoly;i++)
        {
            deriv += coefficients(i)*i*pow(dR/2.,i-1);
            sderiv += coefficients(i)*i*(i-1)*pow(dR/2.,i-2);
            //std::cout << coefficients(i) << "*x**" << i;
            //if (i < nPoly-1) { std::cout << " + ";}
        }

        // Look at various points along the  bond vector
        int npoints = 101;
        float f_at_x[npoints], r_at_x[npoints], ana_field[npoints];
        for (int j=0;j<npoints;j++)
        {
            f_at_x[j] = 0;
            r_at_x[j] = j/(float)(npoints-1)*dR;
            ana_field[j] = 0;
            for (int i=0; i<nPoly; i++)
            {
                if (not isnan(i*pow(j/(float)(npoints-1)*dR,i-1)))
                {
                    f_at_x[j] += coefficients(i)*i*pow(j/(float)(npoints-1)*dR,i-1);
                }
            }
            f_at_x[j] *= -1;
        }
        float midpoint_ana = coulomb_field(pqr,r_at_x,npoints,ana_field);
        float max_abs_rxn_f = 0;
        float max_abs_rxn_f_r = 0;
        for (int i=0;i<npoints;i++)
        {
            if (fabs(f_at_x[i]) > fabs(max_abs_rxn_f))
            {
                max_abs_rxn_f = f_at_x[i];
                max_abs_rxn_f_r = r_at_x[i];
            }
        }
        float max_abs_coul_f = 0;
        float max_abs_coul_f_r = 0;
        for (int i=0;i<npoints;i++)
        {
            if (fabs(ana_field[i]) > fabs(max_abs_coul_f))
            {
                max_abs_coul_f = ana_field[i];
                max_abs_coul_f_r = r_at_x[i];
            }
        }
        float max_abs_total_f = 0;
        float max_abs_total_f_r = 0;
        for (int i=0;i<npoints;i++)
        {
            if (fabs(ana_field[i]+f_at_x[i]) > fabs(max_abs_total_f))
            {
                max_abs_total_f = ana_field[i] + f_at_x[i];
                max_abs_total_f_r = r_at_x[i];
            }
        }
        std::cout << "; finite_rxn poly_rxn coulomb finite_total max_rxn max_coul max_total" << std::endl;
        std::cout << num_der << " " << -deriv << " " << midpoint_ana << " " << num_der + midpoint_ana << " " << max_abs_rxn_f << " " << max_abs_coul_f << " " << max_abs_total_f << std::endl;
    }
    else
    {
        std::cout << num_der << std::endl;
    }

    

    return 0;
}

float coulomb_field(const char *pqr, float *rArray, int nArray, float *fArray)
{
    std::string pqrline,throwAway,atom,resname;
    float charge = 0, midpoint = 0;
    std::vector<float> coords(3,0),CDcoords(3,0),NEcoords(3,0),bVector(3,0);
    std::vector<std::vector<float> > coordArray;
    std::vector<float> qArray;
    std::ifstream file (pqr);
    if (file.is_open())
    {
        while (file.good())
        {
            getline(file, pqrline);
            if (not pqrline.empty())
            {
                std::stringstream linestream(pqrline);
                linestream >> throwAway >> throwAway >> atom >> resname >> throwAway >> coords[0] >> coords[1] >> coords[2] >> charge;
                if (std::strncmp(atom.c_str(),"DUM",4) != 0)
                {
                    if (std::strncmp(resname.c_str(),"CNC",4) != 0 )
                    {
                        coordArray.push_back(coords);
                        qArray.push_back(charge);
                    }
                    else if (std::strncmp(atom.c_str(),"CD",3) == 0)
                    {
                        CDcoords[0] = coords[0];
                        CDcoords[1] = coords[1];
                        CDcoords[2] = coords[2];
                    }
                    else if (std::strncmp(atom.c_str(),"NE",3) == 0)
                    {
                        NEcoords[0] = coords[0];
                        NEcoords[1] = coords[1];
                        NEcoords[2] = coords[2];
                    }
                    else if (std::strncmp(atom.c_str(),"SG",3) != 0 && std::strncmp(atom.c_str(),"CB",3) != 0 && std::strncmp(atom.c_str(),"HB1",4) != 0 && std::strncmp(atom.c_str(),"HB2",4) != 0)
                    {
                        coordArray.push_back(coords);
                        qArray.push_back(charge);
                    }
                }
            }
        }
        float bL2 = 0;
        for (int i=0;i<3;i++)
        {
            bVector[i] = NEcoords[i] - CDcoords[i];
            bL2 += bVector[i]*bVector[i];
        }
        float bondLength = sqrt(bL2);
        for (int i=0;i<3;i++)
        {
            bVector[i] /= bondLength;
        }

        for (int i=0;i<nArray;i++)
        {
            float x = 0;
            float y = 0;
            float z = 0;
            float r = 0;
            float r2 = 0;
            float rr3 = 0;
            float pfield = 0;
            std::vector<float> field(3,0);
            std::vector<float> point(3,0);
            for (int j=0;j<3;j++)
            {
                point[j] = CDcoords[j] + rArray[i]*bVector[j];
            }
            for (int j=0; j<(int)qArray.size(); j++)
            {
                x = point[0] - coordArray[j][0];
                y = point[1] - coordArray[j][1];
                z = point[2] - coordArray[j][2];
                r2 = x*x+y*y+z*z;
                r = pow(r2,0.5);
                rr3 = 1.0 / (r*r2);
                field[0] += rr3 * qArray[j] * x;
                field[1] += rr3 * qArray[j] * y;
                field[2] += rr3 * qArray[j] * z;
            }
            // Project field along bond vector
            for (int j=0; j<3;j++)
            {
                pfield += field[j]*bVector[j] * cfac;
            }
            fArray[i] = pfield;
            if (rArray[i]*2 == bondLength)
            {
                midpoint = pfield;
            }
        }
    }
    return midpoint;
}


