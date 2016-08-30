//
// Created by Gleb Kichaev on 10/21/15.
//

#ifndef PAINTOR_3_0_FUNCTIONS_COREMODEL_H
#define PAINTOR_3_0_FUNCTIONS_COREMODEL_H

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
#include "Functions_IO.h"
#include <unordered_map>
#include <random>
#include <math.h>

using namespace Eigen;
using namespace std;
//using namespace nlopt;

struct Config_Prob
{
    string config;
    double value;
};

void BuildCausalVector(VectorXd& vec2Build , VectorXd& index);
void BuildCausalVector(VectorXd& vec2Build , VectorXd& index, vector<int>& index_as_vector);
int NextCombo(VectorXd& c, int k, int n) ;
inline double Prior_Snp_Probabilty(VectorXd& beta , VectorXd& aij);
double Prior_CausalSet_Probabilty(VectorXd& priorJ,VectorXd& beta, VectorXd& C_vector);
double MVN_Density_ZeroMean(VectorXd& x, MatrixXd& sigma);
VectorXd Zscores2Post(VectorXd& Zs);
double LogSum(double val1, double val2);
vector<double> eigen2vec(VectorXd &vec);
vector<vector<double> > Stack_Matrices(vector<MatrixXd> &mats);
VectorXd vector2eigen(vector<double> &vec);
void string_to_int(string& input_string, vector<int>& output_int);
int Sample_Int(mt19937 & generator ,vector<double>& weights);
int Next_Move(int set_size, double move_probability, mt19937& generator);
double Get_Rand_Unif( mt19937& generator);
int Next_Deletion(VectorXd & zscore_squared, vector<int>& current_set, mt19937& generator);
int Next_Addition(VectorXd & zscore_squared, MatrixXd& LDsquared_complement, vector<int>& current_set, mt19937& generator);
double Calculate_LogBayesFactor(VectorXd& zscore, MatrixXd& ld_matrix, vector<int>& causal_config, double prior_variance);
double Calculate_LogPrior(VectorXd& per_snp_prior,VectorXd& gammas, vector<int> causal_set);
void Compute_SNP_Priors(MatrixXd& annotation_matrix,VectorXd& gammas, VectorXd& per_snp_prior);
void CasualSet_To_String(string& causal_string, vector<int>& causal_set);
void Locus_Sampler(VectorXd& marginal, VectorXd& zscores, VectorXd& gammas, MatrixXd& annotations,  MatrixXd& ld_matrix, double& fullLikeli, double prior_variance, double  move_probability, int num_samples, int sampling_seed);
void Marginalize_Sets(unordered_map<string,double>& causal_posteriors, VectorXd& marginals, double& running_sum);
double Estep(vector<VectorXd> &Zscores, VectorXd &betas, vector<MatrixXd> &Aijs, vector<MatrixXd>& ld_matrices, CausalProbs &E_out, double prior_variance,double  move_probability, int num_samples,int sampling_seed);
double Estep(vector<vector<VectorXd>> &Zscores, VectorXd &betas, vector<MatrixXd> &Aijs, vector<vector<MatrixXd>>& ld_matrices, CausalProbs &E_out, double prior_variance, int num_samples, int sampling_seed);
double Estep(vector<vector<VectorXd>> &Zscores, VectorXd &betas, vector<MatrixXd> &Aijs, vector<vector<MatrixXd>>& ld_matrices, CausalProbs &E_out, double prior_variance,int  max_causals);

double EM_Run(CausalProbs &probabilites, int iter_max, vector<vector<VectorXd>> &Zscores,  VectorXd &beta_int, vector<MatrixXd> &Aijs,vector<vector<MatrixXd>> &ld_matrix , double prior_variance, int num_samples, int sampling_seed);
double EM_Run(CausalProbs &probabilites, int iter_max, vector<vector<VectorXd>> &Zscores,  VectorXd &beta_int, vector<MatrixXd> &Aijs,vector<vector<MatrixXd>> &ld_matrix , double prior_variance,int max_causals);
void Stack_EigenMatrices(vector<MatrixXd> &mats, MatrixXd& stacked_matrices);
void Eigen2Vector(vector<double> & vector_in, VectorXd& eigen_out);
void Gradient_Ascent(VectorXd& current_gammas, VectorXd& new_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, double gradient_tolerance, int max_iterations, VectorXd& return_values);

void Sample_Causal_Vector(VectorXd& causal_prior, mt19937& generator, vector<int>& causal_set);
double Calc_Importance_logWeight(VectorXd& causal_prior, vector<int>& causal_set);

void Locus_Importance_Sampler(VectorXd& marginal, vector<VectorXd>& zscores, VectorXd& gammas, MatrixXd& annotations, vector<MatrixXd>& ld_matrix, double& fullLikeli, double prior_variance, int num_samples, int sampling_seed);

void PrintVector(vector<int>& in_vec);
void Marginalize_Sets_Importance(unordered_map<string,double>& causal_posteriors, VectorXd& marginals, double& running_sum, double total_log_weights);
void Locus_Sampler_Multi(VectorXd& marginal, vector<VectorXd>& zscores, VectorXd& gammas, MatrixXd& annotations,  vector<MatrixXd>& ld_matrix, double& fullLikeli, double prior_variance, double  move_probability, int num_samples, int sampling_seed);
double Importance_Expectation( unordered_map<string,double>& causal_posteriors, double total_log_weights);
double Importance_Expectation( vector<Config_Prob>& causal_posteriors, double total_log_weights);
void Marginalize_Sets(vector<Config_Prob>& causal_posteriors, VectorXd& marginals, double& running_sum);
void Marginalize_Sets_Importance(vector<Config_Prob>& causal_posteriors, VectorXd& marginals, double total_log_weights);
void Marginalize_Importance_Sets(vector<vector<int>>& sampled_causal_sets,  vector<double>& log_posterior_weight,  double log_normalizer, VectorXd& marginal_probs);
double Calc_Importance_logWeight(VectorXd& per_snp_prior, VectorXd& proposal_prior, vector<int>& causal_set);
double PreCompute_Enrichment(int iter_max, vector<vector<VectorXd>> &Zscores,  VectorXd &gamma_intitial, vector<MatrixXd> &annotations,vector<vector<MatrixXd>> &ld_matrix , double prior_variance, int max_causals);
double PAINTOR_Importance_Sampling(vector<vector<VectorXd>> &Zscores, VectorXd &betas, vector<MatrixXd> &Aijs, vector<vector<MatrixXd>>& ld_matrices, CausalProbs &E_out, double prior_variance, int num_samples, int sampling_seed);
double log_normal_density_zero_mean(double x, double sigma_sq);
double determine_optimal_prior_variance(VectorXd& zscore_vector);


#endif //PAINTOR_3_0_FUNCTIONS_COREMODEL_H
