//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once                                                             

#include "Global.h" 
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <string>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <sstream>
#include <cstring>

using namespace std;
using std::complex;
using std::vector;



int         numProcess;
int         myRank;
double TrackingTime=0.E0;


double cpuSecond()
{
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return((double) tp.tv_sec + (double) tp.tv_usec*1e-6 );
}

void sddsplot()
{
    //string command ="sddsplot -col=Turns,\'\(averAllBunchY,averAllBunchX\)\' result.sdds -graph=l,v -leg";
    string command ="python3.6 GrowthRate.py";
    char *cstr = new char[command.length()+1];
    strcpy(cstr,command.c_str());
    system(cstr);
    delete [] cstr;
}

void RMOutPutFiles()
{
    string command ="rm -rf *.sdds";
    char *cstr = new char[command.length()+1];
    strcpy(cstr,command.c_str());
    system(cstr);
    delete [] cstr;
}


double vectorDisNorm2(vector<double> a,vector<double> b)
{
    double sum=0;
    for(int i=0;i<a.size();i++)
    {
        sum += pow(a[i] - b[i],2); 
    }
    sum =sqrt(sum);

    return sum;
}


double Gaussrand(double rms, double aver,double randomIndex)
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
     
     srand(time(0)+randomIndex*10);
     
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
             
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
         
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
         
    phase = 1 - phase;
 
 
    X =  X * rms + aver;
    return X;
}


void gsl_matrix_mul(gsl_matrix *a,gsl_matrix *b,gsl_matrix *c)
{
    int dimAx = a->size1;
    int dimAy = a->size2;

    int dimBx = b->size1;
    int dimBy = b->size2; 

    int dimCx = c->size1;
    int dimCy = c->size2; 

    if(dimAy!= dimBx)
    {
        std::cerr<<"error of input matrix for multipilication"<<std::endl;
    }

    if( dimAx!= dimCx && dimBy != dimCy)
    {
        std::cerr<<"error of output matrix for multipilication"<<std::endl;
    }


      
    for (size_t i=0;i<a->size1;i++)
    {
        for (size_t j=0;j<b->size2;j++)
        {
            double sum=0.0;
            for (size_t k=0;k<b->size1;k++)
            {
                sum+=gsl_matrix_get(a,i,k)*gsl_matrix_get(b,k,j);
            }
            gsl_matrix_set(c,i,j,sum);
        }
    }
}

void PrintGSLMatrix(gsl_matrix *mat)
{
    int row = mat->size1;
    int col = mat->size2;
    for(int j=0;j<row;j++)
    {
        for(int k=0;k<col;k++) printf("%15.7f\t",mat->data[j * mat->tda+k]);
        cout<<endl;       
    }
    cout<<endl;
}

void gsl_matrix_inv(gsl_matrix *a)
{

    if(gsl_matrix_isnull(a)!=1);
    {
        size_t n=a->size1;
        size_t m=a->size2;
        gsl_matrix *temp1=gsl_matrix_calloc(n,n);
        gsl_matrix_memcpy(temp1,a);

        gsl_permutation *p=gsl_permutation_calloc(n);
        int sign=0;
        gsl_linalg_LU_decomp(temp1,p,&sign);
        gsl_matrix *inverse=gsl_matrix_calloc(n,n);

        gsl_linalg_LU_invert(temp1,p,inverse);
        gsl_matrix_memcpy(a,inverse);

        gsl_permutation_free(p);
        gsl_matrix_free(temp1);
        gsl_matrix_free(inverse);
    }
}

void GetInvMatrix(vector<vector<double>> &mat)
{
		
}


double get_det(gsl_matrix *A)
{
    double det=0.0; 
    
    if(gsl_matrix_isnull(A)) return det;
 
    int n = A->size1;
    gsl_permutation *p = gsl_permutation_calloc(n);
    gsl_matrix *tmpA = gsl_matrix_calloc(n, n);
    int signum;
    gsl_matrix_memcpy(tmpA, A);
    gsl_linalg_LU_decomp(tmpA, p, &signum);
    det = gsl_linalg_LU_det(tmpA, signum);
    gsl_permutation_free(p);
    gsl_matrix_free(tmpA);
    return det;
}

double gsl_get_trace(gsl_matrix * A)
{   
    int n = A->size1;
    double trace = 0.E0;
    for(int i=0;i<n;i++)
    {
        trace += gsl_matrix_get(A,i,i);
    }
    return trace;
}


int StringVecSplit(string str, vector<string> &strVec)
{
   if(strVec.size()>0) strVec.clear();
   str = str + " ";
   
   int position = str.find('=');
      
   string temp="";
    
   for(int i=0;i<position;i++)
   {
        if(!isblank(str[i]))
        {
            temp.push_back(str[i]);
        }
   }
    
   transform(temp.begin(), temp.end(), temp.begin(), ::tolower);
    
   strVec.push_back(temp);
   
   temp="";
   
    for(int i=position+1;i<str.size();i++)
    {
        /*
        if(str[i]!=' ')
        {
            temp.push_back(str[i]);
        }
        */

        if(!isblank(str[i]))
        {
            int  j=i;
            while(!isblank(str[j]))
            {
                j++;
            }
            strVec.push_back(str.substr(i,j - i));
            i=j;
        }

    }

    strVec.push_back(temp);

    return 0;
}



int StringSplit(string s,vector<string> &vec)
{

  if(vec.size()>0)
  vec.clear();

  for(int i=0;i<s.length();i++)
  {
    
    if(!isblank(s[i]))
    {
        int  j=i;
        while(!isblank(s[j]))
        {
            j++;
        }
        vec.push_back(s.substr(i,j - i));
        i=j;
    }

  }
  return 0;
}

int StringSplit2(string s, vector<string> &vec)
{

    vec.clear(); 
    stringstream iss(s);

	for(s; iss >> s; ) 
	    vec.push_back(s);
	
    return 0;
}

string convertToString(char* a, int size)
{
    string s(a);
    return s;
}
