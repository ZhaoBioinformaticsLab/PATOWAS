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
// km_calc.cpp : Defines the entry point for the console application.
//
// Copyright: Noble Foundation, 05/23/2016, Writted by Wenchao Zhang.

//#define Serial_Test

#include <iostream>
#include <stdio.h>
#include <iostream> // Cin/ Cout
#include <fstream>  // For Read/Write File
#include <string>
#include <numeric>
#include <vector>  // For using Vector
#include <algorithm> // use vector sort
#include <sstream>  // using sstream
#include <map>
#include <functional>
struct LUT_Entry
{
	int LUT_index;
	int Individual_i;
	int Individual_j;
};

using namespace std;

int parse_cmd_line(int argc, char* argv[], string & Z_File_Name, string & LUT_Allocation_File, int & LUT_Index_Low, int & LUT_Index_High, int & Individual_Num, int & Bin_Num);
double Ka_Calculation(char *LpZ_Vector_i, char *LpZ_Vector_j, int Bin_Num);
double Ka_Calculation(float *LpZ_Vector_i, float *LpZ_Vector_j, int Bin_Num);
double Kaa_Calculation(char *LpZ_Vector_i, char *LpZ_Vector_j, int Bin_Num);
double Kaa_Calculation(float *LpZ_Vector_i, float *LpZ_Vector_j, int Bin_Num);

bool Compare_LUT_by_LUT_Index(LUT_Entry val, LUT_Entry p)
{
	return val.LUT_index< p.LUT_index;
}

int main(int argc, char *argv[])
{
    string Z_File_Name;
	string LUT_Allocation_File_Name;
	int LUT_Index_Low, LUT_Index_High;
    int Individual_Num, Bin_Num;
	string OS_Path_Sep="//";

	// Parse the command line for the inputting
	if(1==parse_cmd_line(argc, argv, Z_File_Name, LUT_Allocation_File_Name, LUT_Index_Low, LUT_Index_High, Individual_Num, Bin_Num))
        return 1;

#ifdef Serial_Test
	cout <<" Serial Test!"<<endl;
#endif

	//Open and Read the Z Matrix from Z.txt.
	ifstream Z_File(Z_File_Name.c_str(), ios::in);

    if(!Z_File.is_open())
	{
		cerr<< Z_File_Name <<" Can't be accessed!"<<endl;
		return 1;
	}

    float *LpZ_Vector;
    LpZ_Vector =NULL;
    string sLine;
    /* Read the ZFile */
    if(Z_File)
    {
       LpZ_Vector=new float[Individual_Num*Bin_Num];

	   int Bin_Count=0;
//	   getline(Z_File, sLine); // By pass the header line
	   while(getline(Z_File, sLine))
	   {
		  if(sLine.empty()) ; // Ignore empty lines
		  else
		  {
            stringstream ss(sLine);
			vector <string> s_v;
			string item;
			char delim =','; //'\t';
		//	char * End;
			// string::size_type sz;

			while(getline(ss, item, delim))
			{
				s_v.push_back(item);
			}

			int Col_Num=s_v.size();
			for (int i=0; i<Col_Num;i++)
			{
			   float value= atof(s_v.at(i).c_str());
			   *(LpZ_Vector+i*Bin_Num +Bin_Count) = value;
			}
          }
		  Bin_Count++;
       }
	}


	// Open and read the LUT_Allocation_File
	vector <LUT_Entry> LUT_List;
	ifstream LUT_Allocation_File(LUT_Allocation_File_Name.c_str(), ios::in);
	if(!LUT_Allocation_File)
	{
		cerr<<LUT_Allocation_File_Name <<"Open Fail!"<<endl;
		return 1;
	}

	while(getline(LUT_Allocation_File, sLine))
	{
		if(sLine.empty()) ; // Ignore empty lines
		else
		{
			stringstream ss(sLine);
			vector <string> s_v;
			string item;
			char delim =',';
			LUT_Entry  LUT_Item;

			getline(ss, item, delim);
			LUT_Item.LUT_index =atoi(item.c_str());

			getline(ss, item, delim);
			LUT_Item.Individual_i =atoi(item.c_str());
			getline(ss, item, delim);
			LUT_Item.Individual_j =atoi(item.c_str());

			LUT_List.push_back(LUT_Item);
		}
	}

	double ka_diag_sum, kaa_diag_sum;
	ka_diag_sum=0;
	kaa_diag_sum=0;

	LUT_Entry  LUT_Item;
	vector <LUT_Entry>::iterator It_Low_LUT, It_High_LUT;

	bool Parallel_Mode= (LUT_Index_High>=LUT_Index_Low)&&(LUT_Index_Low>=0)&&(LUT_Index_High>0)&&((LUT_Index_High-LUT_Index_Low)<=LUT_List.size());

	if(Parallel_Mode)
	{
	   LUT_Item.LUT_index       = LUT_Index_Low;
	   It_Low_LUT               = lower_bound(LUT_List.begin(), LUT_List.end(), LUT_Item, Compare_LUT_by_LUT_Index);
	   LUT_Item.LUT_index       = LUT_Index_High;
	   It_High_LUT              = upper_bound(LUT_List.begin(), LUT_List.end(), LUT_Item, Compare_LUT_by_LUT_Index);
	}
	else
	{
	   It_Low_LUT               = LUT_List.begin();
	   It_High_LUT              = LUT_List.end();
	}

        //define and initialize the 6 kinship marix with 0.0
	double Ka, Kaa;
    Ka =0.1;
    Kaa=0.1;

#ifdef Serial_Test
    double *Ka_Array=new double[Individual_Num*(Individual_Num+1)/2];
    double *Kaa_Array=new double[Individual_Num*(Individual_Num+1)/2];
    long Index_shift =0;
#endif

	for (vector <LUT_Entry>::iterator It_LUT =It_Low_LUT; It_LUT!=It_High_LUT;It_LUT++)
	{
	    int i_n= It_LUT->Individual_i;
        int j_n =It_LUT->Individual_j;

        if(Z_File.is_open())
        {
            Ka = Ka_Calculation(LpZ_Vector+i_n*Bin_Num, LpZ_Vector+j_n*Bin_Num, Bin_Num);
            Kaa= Kaa_Calculation(LpZ_Vector+i_n*Bin_Num, LpZ_Vector+j_n*Bin_Num, Bin_Num);
        }

	    // Output the calculation result into each cpu's corresponding file
	    cout << i_n <<"," <<j_n<<"," << Ka<<"," <<Kaa << endl;

#ifdef Serial_Test
        Ka_Array[Index_shift] = Ka;
		Kaa_Array[Index_shift]= Kaa;
		Index_shift++;
#endif

	    if(i_n==j_n)
	    {
		  ka_diag_sum  +=Ka;
          kaa_diag_sum +=Kaa;
	     }
	}

        cout << -1 <<"," <<-1 <<"," << ka_diag_sum<<","<<kaa_diag_sum << endl;
        if(LpZ_Vector!=NULL) delete []LpZ_Vector;

#ifdef  Serial_Test
		ofstream File_Ka, File_Kaa;
		File_Ka.open("//root//Desktop//PATOWAS (Pipeline for Analysis of Trait by Omics Wide Association Studies)//Data_In_Out//Ka_linux.txt");
		File_Kaa.open("//root//Desktop//PATOWAS (Pipeline for Analysis of Trait by Omics Wide Association Studies)//Data_In_Out//Kaa_linux.txt");
        int Count_Shift =0;
		for (int i_row=0; i_row<Individual_Num; i_row++)
		{
			for (int i_col=0; i_col<=i_row; i_col++)
			{
				File_Ka  << Ka_Array[Count_Shift]*Individual_Num/ka_diag_sum;     //Ka_Array[i_col*Individual_Num+i_row-Trangle_Shift]*Individual_Num/ka_diag_sum;
				File_Kaa << Kaa_Array[Count_Shift]*Individual_Num/kaa_diag_sum;//  Kaa_Array[i_col*Individual_Num+i_row-Trangle_Shift]*Individual_Num/kaa_diag_sum;

				if(i_col!=i_row)
				{
					File_Ka <<",";
					File_Kaa <<",";
				}
				Count_Shift++;
                //Trangle_Shift += i_col+1;
			}
			File_Ka <<endl;
			File_Kaa <<endl;
		}
		File_Ka.close();
		File_Kaa.close();

		delete []Ka_Array;
		delete []Kaa_Array;
#endif

	return 0;
}

// This function is used to parse
int parse_cmd_line(int argc, char* argv[], string & Z_File_Name, string & LUT_Allocation_File, int & LUT_Index_Low, int & LUT_Index_High, int & Individual_Num, int & Bin_Num)
{
	if((argc!=2)&&(argc <13))
	{
		cerr <<"Please specific the full name of additive omics file and LUT allocation file, usually named as Z.txt and LUT_allication.csv!" <<endl;
		cerr <<"Please also specific the Individual_Num and Bin_Num !" <<endl;
		return 1;
	}
	else
	{
	   for (int i=1; i<argc;i++)
       {
	       string arg=argv[i];
	       if((arg=="-z")||(arg=="-Z"))
	       {
                     if(i+1<argc)
		     {
		       Z_File_Name =argv[i+1];
		     }
		     else
		     {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the full name of additive genotype file!" <<endl;
			   return 1; // Parse cmd_line fail
		      }
		    }

	        if((arg=="-a")||(arg=="-A"))
	        {
               if(i+1<argc)
		       {
		           LUT_Allocation_File =argv[i+1];
		       }
		       else
		       {
			       cerr<<"Parse cmd_line fail, Need to clearly specific the full name of LUT_Allocation File!" <<endl;
			       return 1; // Parse cmd_line fail
		       }
	       }

	       if((arg=="-l")||(arg=="-L"))
	       {
             if(i+1<argc)
		     {
		           LUT_Index_Low =atoi(argv[i+1]);
		     }
		     else
		     {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the lower value of the LUT_index for the current cpu node!" <<endl;
			   return 1; // Parse cmd_line fail
		     }
	       }

	       if((arg=="-h")||(arg=="-H"))
	       {
             if(i+1<argc)
		     {
		          LUT_Index_High =atoi(argv[i+1]);
		     }
		     else
		     {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the higher value of the LUT_index for the current cpu node!" <<endl;
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
			   cerr<<"Parse cmd_line fail, Need to clearly specific the Bin Number(marker number)!" <<endl;
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

		   if((arg=="-u")||(arg=="-U"))
	       {
			     cout<<"Welcome to use this program to do Kinship Matrix Calculation" <<endl;
				 cout << "The usuage of input parameter arguments are listed as followings:" <<endl;
				 cout << "-u or -U: Output this help usuage message" <<endl;
				 cout << "-z or -Z: The full name of Omics Z(Additive effect) file"<<endl;
				 cout << "-a or -A: The full name of LUT(Invidual xIndividual) allocation file"<<endl;
				 cout << "-l or -L: The lower index of LUT allocatied for the current cpu node"<<endl;
				 cout << "-h or -H: The higher index of LUT allocatied for the current cpu node"<<endl;
				 cout << "-i or -I: The Individual number" <<endl;
				 cout << "-b or -B: The Bin number" <<endl;
			     return 1; // Parse cmd_line fail
		    }

	    }
	}
    return 0;
}

// This function is used to calculate the Ka
double Ka_Calculation(char *LpZ_Vector_i, char *LpZ_Vector_j, int Bin_Num)
{
    double result=0.0;
    for (int m=0; m<Bin_Num; m++)
    {
        char Z_i = *(LpZ_Vector_i +m);
        char Z_j = *(LpZ_Vector_j +m);
        result += Z_i*Z_j;
    }
    return result;
}

double Ka_Calculation(float *LpZ_Vector_i, float *LpZ_Vector_j, int Bin_Num)
{
    double result=0.0;
    for (int m=0; m<Bin_Num; m++)
    {
        float Z_i = *(LpZ_Vector_i +m);
        float Z_j = *(LpZ_Vector_j +m);
        result += Z_i*Z_j;
    }
    return result;
}

double Kaa_Calculation(char *LpZ_Vector_i, char *LpZ_Vector_j, int Bin_Num)
{
    double result=0.0;
    for (int m=0; m<Bin_Num; m++)
    {
        char Z_i_m =  *(LpZ_Vector_i+m);
        char Z_j_m =  *(LpZ_Vector_j+m);
        for(int n=m+1; n<Bin_Num; n++)
	{
	   char Z_i_n = *(LpZ_Vector_i+n);
           char Z_j_n = *(LpZ_Vector_j+n);
           char aa_i  = Z_i_m*Z_i_n;
           char aa_j  = Z_j_m*Z_j_n;;
           result += aa_i*aa_j;
	}
    }
    return result;
}

double Kaa_Calculation(float *LpZ_Vector_i, float *LpZ_Vector_j, int Bin_Num)
{
    double result=0.0;
    for (int m=0; m<Bin_Num; m++)
    {
        float Z_i_m =  *(LpZ_Vector_i+m);
        float Z_j_m =  *(LpZ_Vector_j+m);
        for(int n=m+1; n<Bin_Num; n++)
	{
	   float  Z_i_n = *(LpZ_Vector_i+n);
           float  Z_j_n = *(LpZ_Vector_j+n);
           double aa_i  = Z_i_m*Z_i_n;
           double aa_j  = Z_j_m*Z_j_n;;
           result += aa_i*aa_j;
	}
    }
    return result;
}

