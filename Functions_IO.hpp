//
//  Functions_IO.hpp
//  PAINTOR_V3.1
//
//  Created by Gleb Kichaev on 10/25/17.
//  Copyright © 2017 Gleb Kichaev. All rights reserved.
//

#ifndef Functions_IO_hpp
#define Functions_IO_hpp

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "Functions_model.hpp"
//#ifndef EIGEN_USE_BLAS
//#define EIGEN_USE_BLAS
//#endif
#include<Eigen>

using namespace Eigen;
using namespace std;

struct input_file_info{
    string file_directory; // directory where all the files are located
    string locus_name; //root where the Locus file is located
    vector<string> zscore_header; //name of Zscore columns
    vector<string> ld_suffix; // ld suffix which is appendend to locus_name
    string annotation_suffix; // annotation suffix which is appendend to locus name
    vector<string> model_annotations; //names of annotations to use in model
    string locus_file_header;
};


//parse command line
void Parse_command_line(int argc, const char * argv[], execution_parameters& parameters );

//split characters based on a given delimiter 'delim'
vector<string> &split(const string &s, char delim, vector<string> &elems) ;

//check to see if file name given exists
vector<string> split(const string &s, char delim);

//mandatory flags that need to be specified are (1) input file list,
//(2) zscore names in locus file header, (3) ld suffix
bool check_file_exists (const string& name);

//Functions to read input locus
void Read_Locus( Locus &  current_locus, input_file_info& info);

//Read in a single LD matrix
void Read_LD(string &input_directory, string& fname , MatrixXf & ld_matrix);


//Read in the annotation matrix for locus
void Read_Annotations(string &input_directory, string& fname , vector<string>& model_annotations, MatrixXf & out_annotations);


//read in all the input for PAINTOR based list of locus files specified in file_list_name
void Get_all_input(string& file_list_name, vector<Locus>& all_loci,  execution_parameters& params);

//display welcome message and command flags
void Welcome_Message();

//parse inititialization of input
void initialize_gammas_input(vector<Locus>& all_loci, execution_parameters & parameters);

void write_all_posteriors(vector<Locus>& all_loci, execution_parameters & parameters);

void write_other_output(VectorXf& gamma_estimates, double log_bf, execution_parameters & parameters);

void write_estimated_prior_variance(vector<Locus>& all_loci, execution_parameters & parameters);

void display_execution_message(execution_parameters & parameters);

#endif /* Functions_IO_hpp */


