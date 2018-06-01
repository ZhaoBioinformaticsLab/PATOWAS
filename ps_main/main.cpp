/**********************************************************************************************************
    This program is released under The GNU General Public License. You can redistribute it and/or
    modify it under the terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this program.
    If not, see <http://www.gnu.org/licenses/>.

***********************************************************************************************************/

/*     
    C/C++ Code written by Wenchao Zhang <wezhang@noble.org>
    Last Update: 05/25/2018   
*/
#include <iostream>
#include <boost/math/distributions/chi_squared.hpp>

#include <math.h>
#include <float.h>
#include "armadillo"

#ifndef math_max
# define math_max(a, b) (a < b)?(b):(a)
# define math_min(a, b) (a > b)?(b):(a)
#endif
//#define ARMA_USE_ATLAS

#ifdef  __cplusplus
extern "C" {
#endif
/* needed for isnan and isfinite, neither of which are used under C++ */

#define FINITE(x) (x <= DBL_MAX && x >= -DBL_MAX)
#define FALSE_ 0
#define TRUE_ 1

#define F77_CALL(x)	x
#define F77_NAME(x)    F77_CALL(x)

//#define _CRT_SECURE_NO_WARNINGS
/*The following marcos defines are copied from R_ext\BLAS.h*/
#ifndef BLAS_extern
#define BLAS_extern extern
#endif

BLAS_extern void   /* DSCAL - scale a one-dimensional array */
F77_NAME(dscal)(const int *n, const double *alpha, double *dx, const int *incx);

BLAS_extern double /* DDOT - inner product of x and y */
F77_NAME(ddot)(const int *n, const double *dx, const int *incx,
	       const double *dy, const int *incy);

BLAS_extern void   /* DCOPY - copy x to y */
F77_NAME(dcopy)(const int *n, const double *dx, const int *incx,
		double *dy, const int *incy);

BLAS_extern void   /* DAXPY - replace y by alpha*x + y */
F77_NAME(daxpy)(const int *n, const double *alpha,
		const double *dx, const int *incx,
		double *dy, const int *incy);

/*The following marcos defines are copied from R_ext\Linpack.h*/
extern void F77_NAME(dtrsl)(double*, int*, int*, double*, int*, int*);
extern void F77_NAME(dpofa)(double*, int*, int*, int*);
#ifdef	__cplusplus
}
#endif

using namespace arma;
using namespace std;
using namespace boost::math;
using boost::math::chi_squared;
using boost::math::cdf;

#include <boost/math/distributions/chi_squared.hpp>
// http://www.boost.org/doc/libs/1_49_0/libs/math/example/chi_square_std_dev_test.cpp


static int c__1 = 1;
static int c__11 = 11;

int Load_data(int Individual_Num, int Bin_Num, string Filename_Ka, string Filename_Kaa, string Filename_GZ, string Filename_PGA, string Filename_Pheno, double *kk_matrix, float *z_matrix, double *lambda, double *phenotype_vector);
int parse_cmd_line(int argc, char* argv[], string & Filename_Ka, string & Filename_Kaa, string & Filename_GZ, string & Filename_PGA, string & Filename_Pheno, /*long & Marker_Index_Low, long & Marker_Index_High,*/ string & Filename_Result, int & Bin_Num, int & Individual_Num);

int main(int argc, char *argv[])
{

    int Individual_Num ;
	int Bin_Num ;
    long Marker_Index_Low, Marker_Index_High;
	string Filename_Ka,  Filename_Kaa, Filename_GZ, Filename_PGA, Filename_Pheno, Filename_Result;
	double *kk_matrix;
	double *lambda;
	double *phenotype_vector;
    float *z_matrix;
    char delim =',';

	if(1==parse_cmd_line(argc, argv, Filename_Ka, Filename_Kaa, Filename_GZ, Filename_PGA, Filename_Pheno, /*Marker_Index_Low, Marker_Index_High,*/ Filename_Result, Bin_Num, Individual_Num))
		return 1;

/*
#ifdef Serial_Test
	cout <<"  Serial Test!"<<endl;
	string Filename_Result="..//Data_In_Out//LRT_Main_linux.txt";
#endif*/

    int PolyGenic_Num=3;

	kk_matrix       = new double [Individual_Num*(Individual_Num+1)];
	phenotype_vector= new double [Individual_Num];
	lambda          = new double [PolyGenic_Num];
	z_matrix        = new float  [Individual_Num*Bin_Num];

	if(1==Load_data(Individual_Num, Bin_Num, Filename_Ka, Filename_Kaa, Filename_GZ, Filename_PGA, Filename_Pheno, kk_matrix, z_matrix, lambda, phenotype_vector))
    {
       delete []kk_matrix;
	   delete []phenotype_vector;
	   delete []lambda;
	   delete []z_matrix;
       return 1;
    }
    else
	{
		// Load data successfully
		mat kk(Individual_Num, Individual_Num);
		long Trian_Matrix_Shift=Individual_Num*(Individual_Num+1)/2;
	    long Local_Shift =0;
	    for (int i=0;i< Individual_Num;i++)
	    {
	    	for (int j=0; j<=i; j++)
		    {
			   double tem =0.0;
			   for (int i_Para =0; i_Para<PolyGenic_Num-1 ; i_Para++) tem += lambda[i_Para]*kk_matrix[i_Para*Trian_Matrix_Shift+Local_Shift];
			   kk(i,j) =tem;
			   if(i!=j)
			   {
			      kk(j,i) =tem;
			   }
			   Local_Shift ++;
			}
		}

		// Begin to do eigen value decomposition
		vec eig_val;
        mat eig_vec;
        eig_sym(eig_val, eig_vec, kk);  //The armadil output the eigene_values in the increasing order, in R, the eigene values are output as decreasing order.

		mat X0(Individual_Num, 1);
		mat Y(Individual_Num, 1);
		mat Z(Individual_Num, Bin_Num);

		for (int i=0;i<Individual_Num; i++)
		{
			X0(i, 0)=1;
			Y(i, 0)=*(phenotype_vector+i);
			for (int j=0; j<Bin_Num; j++ )
			Z(i, j) = *(z_matrix+i*Bin_Num+j);
		}

        vec W(eig_val.n_elem);
        mat ZU=trans(eig_vec)*Z;
		mat YU=trans(eig_vec)*Y;
		mat XU=trans(eig_vec)*X0;

		for (int i=0; i<eig_val.n_elem; i++)
		{
		    W(i)=sqrt(1.0/(1.0+eig_val(i)));
		    ZU.row(i) = ZU.row(i)*W(i);
            YU.row(i) = YU.row(i)*W(i);
            XU.row(i) = XU.row(i)*W(i);
		}

		long Z_Size=Individual_Num*Bin_Num;

//      if((Marker_Index_Low<0)||(Marker_Index_Low>=Bin_Num))  Marker_Index_Low = 0;
        Marker_Index_Low = 0;
//      if((Marker_Index_High<0)||(Marker_Index_High>=Bin_Num))  Marker_Index_High =Bin_Num-1;
        Marker_Index_High =Bin_Num-1;

        ofstream Result_File; // used to output the Estimated 1D P-value arry
		Result_File.open(Filename_Result.c_str());
//		Result_File<< "<Marker #"<< delim<<"P-Value" <<endl;

		for (int i_marker=Marker_Index_Low; i_marker <= Marker_Index_High; i_marker++)
		{
			mat ZU_col = ZU.col(i_marker);
		//	mat x_join_rows= join_rows(XU,ZU.col(i_marker));

			mat X   = join_rows(XU, ZU.col(i_marker));
            mat XPX = inv(trans(X)*X);
            mat XPY = trans(X)*YU;
            mat beta= XPX*XPY;
            mat R   = YU- X*beta;
            mat R_Square =trans(R)*R;
            double sigma2 =accu(R_Square)*1.0/(Individual_Num-2);
        //    if (R_Square.n_elem==1)  sigma2 = R_Square(0,0)*1.0/(Individual_Num-2);
            mat covb   =XPX *sigma2;
            double b0  =beta(0,0);
            double b1  =beta(1,0);
            double v1  =covb(1,1);
            double wald=b1*b1/v1;

            chi_squared dist(1);

		    double p_main = 1-cdf(dist, wald);
//		    cout << i_marker+1 <<delim << p_main<<endl;
            Result_File<< i_marker+1 <<delim << p_main <<endl;

		}
            Result_File.close();
            delete []kk_matrix;
	        delete []phenotype_vector;
	        delete []lambda;
	        delete []z_matrix;
	}

   return 0;
}

int parse_cmd_line(int argc, char* argv[], string & Filename_Ka, string & Filename_Kaa, string & Filename_GZ, string & Filename_PGA, string & Filename_Pheno,  /* long & Marker_Index_Low, long & Marker_Index_High,*/  string & Filename_Result, int & Bin_Num, int & Individual_Num)
{
	if((argc!=2)&&(argc <15))
	{
		cerr <<"Please specific the correct parameters for mian effect P value Estimation!" <<endl;
		return 1;
	}
	else
	{
	   for (int i=1; i<argc;i++)
       {
	       string arg=argv[i];
	       if((arg=="-k")||(arg=="-K"))
	       {
              if(i+1<argc)
		      {
		         string Path = argv[i+1];
				 Filename_Ka = Path+"//" +"Ka.txt";
				 Filename_Kaa= Path+"//" +"Kaa.txt";
		      }
		      else
		      {
			     cerr<<"Parse cmd_line fail, Need to clearly specific the Path of kinshipmatrix files!" <<endl;
			     return 1; // Parse cmd_line fail
		      }
		    }

		    if((arg=="-z")||(arg=="-Z"))
	        {
              if(i+1<argc)
		      {
				 Filename_GZ = argv[i+1];
		      }
		      else
		      {
			     cerr<<"Parse cmd_line fail, Need to clearly specific the full name of Genotype Z(Additive effect) file!" <<endl;
			     return 1; // Parse cmd_line fail
		      }
		    }

			if((arg=="-c")||(arg=="-C"))
	        {
              if(i+1<argc)
		      {
				 Filename_PGA = argv[i+1];
		      }
		      else
		      {
			     cerr<<"Parse cmd_line fail, Need to clearly specific the full name of Polygenic Component Analysis result file!" <<endl;
			     return 1; // Parse cmd_line fail
		      }
		    }

            if((arg=="-p")||(arg=="-P"))
	        {
              if(i+1<argc)
		      {
		         Filename_Pheno = argv[i+1];
		      }
		      else
		      {
			     cerr<<"Parse cmd_line fail, Need to clearly specific the full name of phenotype file!" <<endl;
			     return 1; // Parse cmd_line fail
		      }
		   }


		  /*
		  if((arg=="-l")||(arg=="-L"))
          {
                if(i+1<argc)
		        {
		           Marker_Index_Low =atol(argv[i+1]);
		        }
		        else
		        {
			        cerr<<"Parse cmd_line fail, Need to clearly specific the Lower value of Marker index" <<endl;
			        return 1; // Parse cmd_line fail
		        }
           }

           if((arg=="-h")||(arg=="-H"))
           {
                if(i+1<argc)
		        {
		           Marker_Index_High =atol(argv[i+1]);
		        }
		        else
		        {
			        cerr<<"Parse cmd_line fail, Need to clearly specific the High value of Marker index!" <<endl;
			        return 1; // Parse cmd_line fail
		        }
	        }
	        */

           if((arg=="-o")||(arg=="-O"))
	       {
              if(i+1<argc)
		      {
		         Filename_Result = argv[i+1];
		      }
		      else
		      {
			     cerr<<"Parse cmd_line fail, Need to clearly specific the full name of output result file!" <<endl;
			     return 1; // Parse cmd_line fail
		      }
		    }

            if((arg=="-i")||(arg=="-I"))
	        {
                if(i+1<argc)
		        {
		           Individual_Num =atoi(argv[i+1]);
		        }
		        else
		        {
			        cerr<<"Parse cmd_line fail, Need to clearly specific the Individual Number!" <<endl;
			        return 1; // Parse cmd_line fail
		        }
	        }

			if((arg=="-b")||(arg=="-B"))
	        {
                if(i+1<argc)
		        {
		           Bin_Num =atoi(argv[i+1]);
		        }
		        else
		        {
			        cerr<<"Parse cmd_line fail, Need to clearly specific the Bin Num!" <<endl;
			        return 1; // Parse cmd_line fail
		        }
	        }

		    if((arg=="-u")||(arg=="-U"))
	        {
			     cout<<"Welcome to use this program to do Polygenic Component Analysis" <<endl;
				 cout << "The usuage of input parameter arguments are listed as followings:" <<endl;
				 cout << "-u or -U: Output this help usuage message" <<endl;
				 cout << "-k or -K: The path of the 6 avaliable kinship matrix files" <<endl;
				 cout << "-z or -Z: The full name of Genotype Z(Additive effect) file"<<endl;
				 cout << "-c or -C: The full name of Variance Component Analysis result file"<<endl;
                 cout << "-o or -O: The full name of the Output result" <<endl;
				 cout << "-p or -P: The full name of the phenotype file" <<endl;
				/*
				 cout << "-l or -L: The Lower value of Marker index allocatied for the current cpu node"<<endl;
				 cout << "-h or -H: The High value of Marker index allocatied for the current cpu node"<<endl;
				 */
				 cout << "-i or -I: The Individual number" <<endl;
				 cout << "-b or -B: The Bin number" <<endl;
			     return 1; // Parse cmd_line fail
		    }

	    }
	}

	return 0;
}


/* This function is used to load the 2D Kinship matrixes, 2D genotypic matrix, 2 lambda vector,  and 1D phenotype vector from the 6 K-matirx files, 1 genotype files, 1 polygenic analysis result file,  and the phenotype file*/
int Load_data(int Individual_Num, int Bin_Num, string Filename_Ka,  string Filename_Kaa, string Filename_GZ, string Filename_PGA, string Filename_Pheno, double *kk_matrix, float *z_matrix, double *lambda, double *phenotype_vector)
{
	ifstream K_File, Genotype_File, PGA_Result_File, Phenotype_File;
	string sLine;
	long Num_Count =0;
	K_File.open(Filename_Ka.c_str(), ios::in);
	while(getline(K_File, sLine))
	{
		if(sLine.empty()) ; // Ignore empty lines
		else
		{
            stringstream ss(sLine);
			vector <string> s_v;
			string item;
			char delim =',';
			char * End;
			string::size_type sz;

			while(getline(ss, item, delim))
			{
			   float value= atof(item.c_str());
			   *(kk_matrix+Num_Count) = value;
			   Num_Count++;
			}
		}
	}
	K_File.close();

	K_File.open(Filename_Kaa.c_str(), ios::in);
	while(getline(K_File, sLine))
	{
		if(sLine.empty()) ; // Ignore empty lines
		else
		{
            stringstream ss(sLine);
			vector <string> s_v;
			string item;
			char delim =',';
			char * End;
			string::size_type sz;

			while(getline(ss, item, delim))
			{
			   float value= atof(item.c_str());
			   *(kk_matrix+Num_Count) = value;
			   Num_Count++;
			}
		}
	}
	K_File.close();

	if (Num_Count!=Individual_Num*(Individual_Num+1))
	{
		cerr <<"Load Kinship Matrix failed!" <<endl;
		return 1;
	}

	//Load the genotype file
	Num_Count =0;
	int Bin_Count=0;
	Genotype_File.open(Filename_GZ.c_str(), ios::in);
	while(getline(Genotype_File, sLine))
	{
		if(sLine.empty()) ; // Ignore empty lines
		else
		{
            stringstream ss(sLine);
			vector <string> s_v;
			string item;
			char delim =','; //'\t';
			char * End;
			// string::size_type sz;

			while(getline(ss, item, delim))
			{
				s_v.push_back(item);
			}

			int Col_Num=s_v.size();
			for (int i=0; i<Col_Num;i++)
			{
			   float value= atof(s_v.at(i).c_str());
			   *(z_matrix+i*Bin_Num +Bin_Count) = value;
			   Num_Count++;
			}
		}
		Bin_Count++;
    }
	Genotype_File.close();

	if (Num_Count!=Individual_Num*Bin_Num)
	{
		cerr <<"Load Genotype Z_File failed!" <<endl;
		return 1;
	}

	//Load the polygenic analysis result file
	Num_Count =0;
	PGA_Result_File.open(Filename_PGA.c_str(),ios::in);
	getline(PGA_Result_File, sLine); // Read the first line.
    while(getline(PGA_Result_File, sLine))
	{
		if(sLine.empty()) ; // Ignore empty lines
		else
		{
		    stringstream ss(sLine);
			vector <string> s_v;
			string item;
			char delim =',';
			char * End;
			// string::size_type sz;

			while(getline(ss, item, delim))
			{
				s_v.push_back(item);
			}

			int Col_Num=s_v.size();
			double Residue= atof(s_v.back().c_str());
			for (int i=0; i<Col_Num;i++)
			{
			   float value= atof(s_v.at(i).c_str());
			   *(lambda+Num_Count) = 1.0*value/Residue;
			   Num_Count++;
			}
		}
	}

	PGA_Result_File.close();
	if(Num_Count!=3)
	{
		cerr <<" lamda Vector load failed!" <<endl;
		return 1;
	}

	//Load the phenotype file
	Num_Count =0;
	Phenotype_File.open(Filename_Pheno.c_str(), ios::in);
	while(getline(Phenotype_File, sLine))
	{
		if(sLine.empty()) ; // Ignore empty lines
		else
		{
		   float value= atof(sLine.c_str());
		   *(phenotype_vector+Num_Count) = value;
		   Num_Count++;
		}
	}
	Phenotype_File.close();

	if(Num_Count!=Individual_Num)
	{
		cerr <<"Phenotype Vector load failed!" <<endl;
		return 1;
	}

	return 0;
}

