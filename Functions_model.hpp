//
//  Functions_model.hpp
//  PAINTOR_V3.1
//
//  Created by Gleb Kichaev on 10/25/17.
//  Copyright © 2017 Gleb Kichaev. All rights reserved.
//

#ifndef Functions_model_hpp
#define Functions_model_hpp

#include <cstdio>
#include <Eigen>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen>
#include <iterator>
#include <string>
#include <iomanip>
#include <algorithm>
#include <unordered_map>
#include <random>
#include <cmath>


using namespace Eigen;
using namespace std;


struct optimization_input{
    vector<double> current_posterior;
    vector<vector<double> > all_annotations;
};


//default executation parameters
struct execution_parameters{
    string in_dir = "./";
    string input_files = "input.files";
    string out_dir = "./";
    vector<string> model_annotations;
    int max_causal = 2;
    string gammaName = "Enrichment.Values";
    string likeli_name= "Log.BayesFactor";
    int maxIter= 15;
    string annotation_suffix = "annotations";
    string results_suffix = "results";
    vector <string> ld_suffix;
    vector<string> zscore_header;
    string single_post_flag = "False";
    int sampling_seed = time(NULL);
    string gamma_initial;
    int enumerate_flag = 0;
    int mcmc_flag = 0;
    int burn_in = 5000;
    int max_samples = 50000;
    int initialized_gammas=0;
    float prop_ld_eigenvalues = 0.95;
    int num_chains = 5;

};

//Class definition of a locus
//Each locus has its requires, Z-scores, LD, annotations
//Initialized with the gamma parameters (enrichments) along with the prior variance
class Locus{
public:
    Locus(){
        vector<VectorXf> z_scores;
        vector<MatrixXf> ld_matrix;
        MatrixXf annotations;
        vector<double> prior_variance;
        vector<int> config;
        VectorXf posterior_probability;
        VectorXf gamma_parameters;
        vector<string> snp_info;
        default_random_engine generator;
        size_t num_snps;
        size_t num_pops;
        double expected_log_likeli;
    };

    
    //initialzing functions
    void set_annotations(MatrixXf & update);
    void set_Zscore(VectorXf& update);
    void set_LD(MatrixXf& update);
    void add_info(string& curr_snp);
    void estimate_prior_variance(float proportion_eigenvalues);
    void set_num_snps(unsigned long number_of_snps);
    void set_total_pops(unsigned long number_of_pops);
    void set_locus_name(string& name);

    
    //acesss private member variables
    size_t get_num_snps();
    void set_gamma(VectorXf& gamma_update);
    MatrixXf get_annotations();
    VectorXd get_posterior_prob();
    double get_prior_variance(int index);
    VectorXf get_gamma();
    
    //model functions
    double enumerate_posterior(int max_number_of_causals);
    double  mcmc_posterior(size_t burn_in, size_t max_samples);
    double compute_log_prior();
    double compute_log_prior(int flip_bit, bool already_causal);
    void initialize_causal_config(int max_causals);
    void initialize_causal_config();
    void causal_config_to_vector(int max_causals);
    void update_log_posterior(double current_evaluation);
    void normalize_log_posterior(double normalizer);
    void reset_posterior();
    void write_results(string& fname);
    double mcmc_posterior(size_t burn_in, size_t max_samples, int num_chains);
    void next_configuration(int next_bit, bool done_burning);

    //public member variables
    uniform_real_distribution<double> uniform_double_generator;
    uniform_int_distribution<int> uniform_int_generator;
    default_random_engine generator;
    size_t num_pops;
    vector<int> ld_pcs;

private:
    string locus_name;
    vector<VectorXf> z_scores;
    vector<MatrixXf> ld_matrix;
    size_t num_snps;
    MatrixXf annotations;
    vector<double> prior_variance;
    vector<int> causal_index_vector;
    VectorXf causal_index;
    VectorXf causal_config_bit_vector;
    VectorXf annotation_eval;
    VectorXd posterior_probability;
    VectorXf gamma_parameters;
    vector<string> snp_info;
    double expected_log_likeli;

};



//use fixed effects estimator with truncated SVD to estimate the prior variance
// on the non-centrality parameter
void all_prior_variance_esimation(vector<Locus> & all_loci, float prop_eigenvalues);

float compute_prior_variance(VectorXf & z_scores, MatrixXf & ld_matrix);

int next_combo(VectorXf& causal_config, int k, int n);

void initialize_gamma_zero(vector<Locus>& all_loci, int number_of_annotations);

double compute_all_posteriors_enumerate(vector<Locus>& all_loci, int max_number_of_cauals);

double compute_all_posteriors_mcmc(vector<Locus>& all_loci, size_t burn_in, size_t max_samples, int num_chains);

void create_optimization_input(vector<Locus> & all_loci, optimization_input& opt_in);

vector<double> eigen_to_vector(VectorXf& eigen);

void update_all_gammas(vector<Locus>& all_loci, vector<double>& updated_gammas);

double run_EM(vector<Locus> & all_loci, execution_parameters& parameters );

int rand_int(size_t min, size_t max);

double rand_double_uniform(size_t min, size_t max) ;

void initialize_generators(vector<Locus>& all_loci, int seed_value);



#endif /* Functions_model_hpp */
