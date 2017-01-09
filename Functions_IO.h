//
// Created by Gleb Kichaev on 10/21/15.
//

#ifndef PAINTOR_3_0_FUNCTIONS_H
#define PAINTOR_3_0_FUNCTIONS_H

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen>
//#include "nlopt.hpp"
#include <iterator>
#include <string>
#include <iomanip>
#include <algorithm>

using namespace Eigen;
using namespace std;
//using namespace nlopt;


struct ObjectiveData{
    vector<double> probs;
    vector<vector<double> > Aijs;
};
struct CausalProbs{
    vector<VectorXd> probs_locs;
    vector<double> probs_stacked;
};





vector<string> &split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);
void Read_Locus(string &input_directory, string& fname, vector<string> & zname,  vector<VectorXd>& all_association_stat, vector<string>& snp_info, string & header);
double Regularize_LD(MatrixXd & ld_mat);
void Read_LD(string &input_directory, string& fname , MatrixXd& ld_matrix);
void Get_All_Input(string& file_list_name, string& directory, vector<string> znames, vector<string>& model_annotations, vector<vector<VectorXd>>& all_locus_statistics, vector<vector<MatrixXd>> &all_ld_matrices,  vector<MatrixXd>& all_annotations, vector<string>& LD_suffix, string& annotation_suffix,
                   vector<vector<string>>& all_SNP_info, vector<string>& all_headers);
void Read_Annotations(string &input_directory, string& fname , vector<string>& model_annotations, MatrixXd& out_annotations, vector<string> & annotation_header_split);
void Write_Posterior(string& out_dir, string & out_name, VectorXd& locus_results, vector<string>& locus_info, string& header );
void Write_All_Output(string& input_files, string& out_dir, string& out_suffix, CausalProbs& results, vector<vector<string>> & all_locus_info,
                      VectorXd & gamma_estimates, string& Gname, double log_likeli, string & Lname, vector<string>& all_headers, vector<string>& annot_names);
void Reshape_Input(vector<vector<VectorXd>>& all_locus_statistics, vector<vector<MatrixXd>> &all_ld_matrices,vector<VectorXd>& reshape_stats, vector<MatrixXd> &reshape_ld );
double Get_Gamma_Zero(vector<vector<VectorXd>>& all_locus_statistics);
bool check_file_exists (const std::string& name);
void Check_Mandatory_Flags(int argc, const char * argv[]);
#endif //PAINTOR_3_0_FUNCTIONS_H
