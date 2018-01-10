//
//  Functions_model.cpp
//  PAINTOR_V3.1
//
//  Created by Gleb Kichaev on 10/25/17.
//  Copyright © 2017 Gleb Kichaev. All rights reserved.
//

#include "Functions_model.hpp"
#include "Functions_optimize.hpp"


//given the current causal configuration index (causal_config) of size  = max # of cauals, give the next configuration
// for up to k causals.
int next_combo(VectorXf& causal_config, int k, int n) {
    for (int i= k; --i >= 0;) {
        if (++causal_config[i] <= n-(k-i)) {
            while (++i < k)
                causal_config[i]= causal_config[i-1]+1;
            return 1;
        }
    }
    return 0;
}

//convert a vector stored in Eigen vector class to a STL vector
vector<double> eigen_to_vector(VectorXf& eigen){
    vector<double> out_vector;
    for(int i = 0; i < eigen.size(); i++){
        out_vector.push_back(eigen(i));
    }
    return out_vector;
}


///Log sum exp trick for adding two numbers in log-space
double log_sum(double val1, double val2){
    double logsum = 0;
    if(val1 > val2){
        logsum = log(1 + exp(val2-val1)) + val1;
    }
    else{
        logsum = log(1 + exp(val1-val2)) + val2;
    }
    
    return(logsum);
}


//return density log MVN with mean zero
double logMVN_zero_mean(VectorXd& x, MatrixXd& sigma){
    MatrixXd b;
    b.setIdentity(sigma.rows(),sigma.cols());
    MatrixXd sigma_inv = sigma.ldlt().solve(b);
    double exponential_term = x.transpose()*(sigma_inv*x);
    double log_det = -.5*log(sigma.determinant());
    return log_det + -.5*exponential_term;
}

double calc_log_bayes_factor(VectorXf& zscore, MatrixXf& ld_matrix, vector<int>& causal_config, double prior_variance){
    //extract ld and zscores for current causal configuration
    size_t set_size = causal_config.size();
    if(set_size == 0){
        return 0; // corresponds to bayes factor of nul model
    }
    MatrixXd ld_c(set_size, set_size);
    VectorXd zscores_c(set_size);
    for(unsigned int i =0; i < set_size; i++){
        zscores_c(i)=zscore(causal_config[i]);
        for(unsigned int j = 0; j < set_size; j++){
            ld_c(i,j)= ld_matrix(causal_config[i], causal_config[j]);
        }
    }
    //evaluate log bayes factor
    MatrixXd causal_variance = ld_c+(prior_variance/set_size)*(ld_c*ld_c);
    double bf_num = logMVN_zero_mean(zscores_c, causal_variance);
    double bf_denom = logMVN_zero_mean(zscores_c, ld_c);
    return bf_num-bf_denom;
}

//compute the log prior for current causal configuration and gamam estimates
double Locus::compute_log_prior(){
    double log_prior = 0;

    for(int i = 0; i < num_snps; i++){
        if(causal_config_bit_vector[i] == 1){
            log_prior += -1*log(1 + exp(annotation_eval[i])); //evaluate log logistic function
        }
        else{
            log_prior += -1*log(1 + exp(-1*annotation_eval[i])); //since 1-logit(x) = logit(-x)
        }
    }
    return log_prior;
}

//compute the log prior for current causal configuration and gamma estimates with a bit flipped
double Locus::compute_log_prior(int flip_bit, bool already_causal){
    double log_prior = 0;
    for(int i = 0; i < num_snps; i++){
        if(causal_config_bit_vector[i] == 1){
            log_prior += -1*log(1 + exp(annotation_eval[i])); //evaluate log logistic function
        }
        else{
            log_prior += -1*log(1 + exp(-1*annotation_eval[i])); //since 1-logit(x) = logit(-x)
        }
    }

    //flip the remove contribution of flip_bit and replace with alternate contribution
    if(already_causal){
        log_prior -= -1*log(1 + exp(annotation_eval[flip_bit]));
        log_prior += -1*log(1 + exp(-1*annotation_eval[flip_bit]));
    }
    else{
        log_prior -= -1*log(1 + exp(-1*annotation_eval[flip_bit]));
        log_prior += -1*log(1 + exp(annotation_eval[flip_bit]));
    }

    return log_prior;
}


//enumerate overall all possible configuration up to max_number_of_causals
//return the sum of (log bayes factors + log prior)
double Locus::enumerate_posterior(int max_number_of_causals){
    //initialize all the variables
    annotation_eval = annotations*gamma_parameters; // compute inner product with with SNP annotation values and current gamma values
    double current_prior;
    initialize_causal_config(max_number_of_causals);
    reset_posterior();
    //iterate over all causal configurations up to max_number_of_causals
    int counter = 1;
    double log_bayes_factor = 0;
    double running_sum = -1e150;
    double current_evaluation = 0;
    while(next_combo(causal_index, max_number_of_causals, num_snps+1)){
        counter ++;
        causal_config_to_vector(max_number_of_causals);
        current_prior = compute_log_prior();
        log_bayes_factor = 0;
        for(int i = 0; i< num_pops; i++){
            log_bayes_factor += calc_log_bayes_factor(z_scores[i], ld_matrix[i], causal_index_vector, prior_variance[i]);
        }
        if(isnan(log_bayes_factor) || isinf(log_bayes_factor)){
            log_bayes_factor = -1e150;
        }
        
        current_evaluation = current_prior + log_bayes_factor;
        running_sum = log_sum(running_sum, current_prior + log_bayes_factor);
        update_log_posterior(current_evaluation);
    }
    normalize_log_posterior(running_sum);
    return(running_sum);
}

double Locus::mcmc_posterior(size_t burn_in, size_t max_samples) {

    //initalizer sampler parameters
   // reset_posterior();
    posterior_probability.setZero();
    initialize_causal_config();
    int next_bit=-1;
    VectorXf configuration_tracker(num_snps);
    configuration_tracker.setZero();
    annotation_eval = annotations*gamma_parameters; //only need to evaluate annotatios once
    expected_log_likeli = -1e150;

    ///run gibb sampler, discard burn-in
    for(int i = 0; i < burn_in+max_samples; i++){
        next_bit = uniform_int_generator(generator);
        if (i > burn_in) {
            next_configuration(next_bit, true);
            configuration_tracker = configuration_tracker + causal_config_bit_vector;
        }
        else{
            next_configuration(next_bit, false);
        }
    }

    //compute posterior probability as proportion of times causal configuration was on
    configuration_tracker = configuration_tracker/max_samples;
    posterior_probability = configuration_tracker.cast<double> ();


    return expected_log_likeli - log(max_samples);
}


double Locus::mcmc_posterior(size_t burn_in, size_t max_samples, int num_chains) {
    //initalizer sampler parameters and reset posterior
    posterior_probability.setZero();
    annotation_eval = annotations*gamma_parameters; //only need to evaluate annotatinos once
    VectorXf avg_config_tracker(num_snps);
    avg_config_tracker.setZero();
    double avg_expect_log_likeli = -1e150;

    for(int j = 0; j < num_chains; ++j) {
        initialize_causal_config(); //restart the chain at null model
        int next_bit = -1;
        VectorXf configuration_tracker(num_snps);
        configuration_tracker.setZero();
        expected_log_likeli = -1e150;
        ///run gibb sampler, discard burn-in
        for (int i = 0; i < burn_in + max_samples; i++) {
            next_bit = uniform_int_generator(generator);
            if (i > burn_in) {
                next_configuration(next_bit, true);
                configuration_tracker = configuration_tracker + causal_config_bit_vector;
            }
            else{
                next_configuration(next_bit, false);
            }
        }

        //compute posterior probability as proportion of times causal configuration was on
        configuration_tracker = configuration_tracker / max_samples;
        avg_config_tracker = avg_config_tracker+configuration_tracker;
        avg_expect_log_likeli = log_sum(avg_expect_log_likeli, expected_log_likeli);
    }

    //average over configurations
    avg_config_tracker = avg_config_tracker/num_chains;
    posterior_probability = avg_config_tracker.cast<double> ();

    return avg_expect_log_likeli - log(num_chains);
}

void Locus::next_configuration(int next_bit, bool done_burning){
    bool currently_causal;
    vector<int> alt_configuration = causal_index_vector;
    if(causal_config_bit_vector(next_bit) == 1){
        currently_causal = true;
        alt_configuration.erase(remove(alt_configuration.begin(), alt_configuration.end(), next_bit), alt_configuration.end());
    } else{
        currently_causal = false;
        alt_configuration.push_back(next_bit);
    }

    //evaluate the two probabilities for the current and alternative configurations
    double log_prior_current = compute_log_prior();
    double log_prior_alt = compute_log_prior(next_bit, currently_causal);
    double log_bayes_factor_current = 0;
    for(int i = 0; i< num_pops; i++){
        log_bayes_factor_current += calc_log_bayes_factor(z_scores[i], ld_matrix[i], causal_index_vector, prior_variance[i]);
    }
    if(isnan(log_bayes_factor_current) || isinf(log_bayes_factor_current)){
        log_bayes_factor_current = -1e150;
    }

    //evaluate alternative configuration
    double log_bayes_factor_alt = 0;
    for(int i = 0; i< num_pops; i++){
        log_bayes_factor_alt += calc_log_bayes_factor(z_scores[i], ld_matrix[i], alt_configuration, prior_variance[i]);
    }
    if(isnan(log_bayes_factor_alt) || isinf(log_bayes_factor_alt)){
        log_bayes_factor_alt = -1e150;
    }

    //determine
    double log_prob_curr_config = log_prior_current + log_bayes_factor_current;
    double log_prob_alt_config = log_prior_alt + log_bayes_factor_alt;
    double log_total = log_sum(log_prob_curr_config,log_prob_alt_config);
    double log_flip_bit_prob = log_prob_alt_config - log_total;

    double random_log_unif = log(uniform_double_generator(generator));

    //accept the bit flip and update causal configuration if the probability is greater than a random uniform
    if(log_flip_bit_prob > random_log_unif){
        if(done_burning){
            expected_log_likeli = log_sum(expected_log_likeli, log_prob_alt_config);
        }

        causal_index_vector = alt_configuration;
        if(currently_causal){
            causal_config_bit_vector(next_bit) = 0;
        }
        else{
            causal_config_bit_vector(next_bit) = 1;
        }
    }
    else{
        if(done_burning) {
            expected_log_likeli = log_sum(expected_log_likeli, log_prob_curr_config);
        }
    }
}

void initialize_generators(vector<Locus>& all_loci, int seed_value){

    size_t num_loci = all_loci.size();
    uniform_real_distribution<double> uniform_double_init(0, 1);
    for(int i = 0; i < num_loci; i++){
        all_loci[i].generator.seed(seed_value);
        int number_of_snps = int(all_loci[i].get_num_snps());
        uniform_int_distribution<int> uniform_int_init(0, number_of_snps - 1);
        all_loci[i].uniform_int_generator = uniform_int_init;
        all_loci[i].uniform_double_generator = uniform_double_init;
    }
}

//update the marginal probabilities at the curretn causal configuration with current_evaluation
void Locus::update_log_posterior(double current_evaluation){
    for(int i = 0; i < causal_index_vector.size(); i++){
        posterior_probability(causal_index_vector[i]) = log_sum(posterior_probability(causal_index_vector[i]), current_evaluation);
    }
}

//normalize the log posteriors and exponentiate to return normal probabilties
void Locus::normalize_log_posterior(double normalizer){
    for(int i = 0; i < num_snps; i++){
        posterior_probability(i) = exp(posterior_probability(i)-normalizer);
    }
}

//reset the posterior probability for the current iteration
void Locus::reset_posterior(){
    VectorXd reset_vector(num_snps);
    reset_vector.setZero();
    posterior_probability = reset_vector;
    for(int i = 0; i < num_snps; i++){
        posterior_probability(i) = -1e150;
    }
}

//convert causal configuration to a vector and reset bit vector
void Locus::causal_config_to_vector(int max_causals){
    causal_index_vector.clear();
    causal_config_bit_vector.setZero();
    for(int i = 0; i < max_causals; i ++){
        if(causal_index[i] > 0){
            causal_index_vector.push_back(causal_index[i]-1);
            causal_config_bit_vector(causal_index[i]-1) = 1;
        }
    }
}

//reset the causal configuration for new round of sampling
void Locus::initialize_causal_config(int max_causals){
    VectorXf temp_vector(max_causals);
    temp_vector.setZero();
    causal_index = temp_vector;
    VectorXf temp_bit_vector(num_snps);
    temp_bit_vector.setZero();
    causal_config_bit_vector = temp_bit_vector;
}

void Locus::initialize_causal_config(){
    VectorXf temp_bit_vector(num_snps);
    temp_bit_vector.setZero();
    causal_config_bit_vector = temp_bit_vector;
    vector<int> no_caus;
    causal_index_vector  = no_caus;
}


//accessor and mutator functions for member variables
void Locus::set_Zscore(VectorXf& update){
    z_scores.push_back(update);
}

void Locus::set_LD(MatrixXf& update){
    ld_matrix.push_back(update);
}

void Locus::add_info(string &curr_snp){
    snp_info.push_back(curr_snp);
}

void Locus::set_annotations(MatrixXf & update){
    annotations = update;
}

void Locus::set_num_snps(unsigned long number_of_snps){
    num_snps = number_of_snps;
}

void Locus::set_total_pops(unsigned long number_of_pops){
    num_pops = number_of_pops;
}

void Locus::set_locus_name(string &name){
    locus_name = name;
}

VectorXf Locus::get_gamma(){
    return gamma_parameters;
}

VectorXd Locus::get_posterior_prob(){
    return posterior_probability ;
}

double Locus::get_prior_variance(int index){
    if(index < prior_variance.size()){
        return prior_variance[index];
    }
    else{
        cout << "Error, index out of bounds" << endl;
        return(-1);
    }
}

MatrixXf Locus::get_annotations(){
    return annotations;
}

void Locus::estimate_prior_variance(float proportion_eigenvalues){
    unsigned long num_sets = ld_matrix.size();
    int max_var = 500;
    for(int i=0; i < num_sets; i++){
        //run eigen decomposition (eigen values are sorted lowest to smallest)
        MatrixXf curr_ld = ld_matrix[i];
        VectorXf curr_z = z_scores[i];
        SelfAdjointEigenSolver<MatrixXf> eigen_decomp(curr_ld);
        MatrixXf D = eigen_decomp.eigenvalues().asDiagonal();
        MatrixXf V = eigen_decomp.eigenvectors();
        
        //determine where to threshold on K based on proportion of variance explained
        //in the LD matrix. Using the fact that Tr(X)  = number of snps =sum(eigen values)
        float total_eigenvalues=0;
        float max_k = 1;
        
        for(int k = num_snps-1; k >= 0; k--){
            total_eigenvalues += D(k,k);
            if(total_eigenvalues/num_snps >= proportion_eigenvalues){
                D(k,k) = 0; //truncate any eigenvectors beyond max_k
            }
            else{
                max_k = num_snps-k;
                D(k,k) = 1/D(k,k); //invert eigenvalue to get a truncated approximation of ld_inv
            }
        }
        MatrixXf trunc_ld_inv = V*(D*V.transpose());
        float curr_prior_variance = curr_z.transpose()*(trunc_ld_inv*curr_z)- max_k  ;
        if(curr_prior_variance > max_var){
            cout << "Warning! The estimated N*h2_g,local for locus " <<  locus_name << " is: " << curr_prior_variance << endl;
            cout << "This may potentially indicate mismatch/error in the LD-matrix" << endl;
        }
        prior_variance.push_back(curr_prior_variance);
        ld_pcs.push_back(max_k);
    }
}

size_t Locus::get_num_snps(){
    return num_snps;
}

void Locus::set_gamma(VectorXf& gamma_update){
    gamma_parameters = gamma_update;
}

//initialize the gammas such that the prior probability to becausal is 1/#average number of SNPs
void initialize_gamma_zero(vector<Locus>& all_loci, int number_of_annotations){
    int total_snps = 0;
    int num_loci = all_loci.size();
    for(int i = 0; i < num_loci; i++){
        total_snps += all_loci[i].get_num_snps();
    }
    float ave_snps = total_snps/num_loci;
    float gamma_zero = log(ave_snps-1);
    VectorXf gamma_initial(number_of_annotations);
    gamma_initial.setZero();
    gamma_initial(0) = gamma_zero;
    for(int i = 0; i < num_loci; i++){
        all_loci[i].set_gamma(gamma_initial);
    }
}

//iterate over all loci and estimate the prior standardized effect size variance
void all_prior_variance_esimation(vector<Locus> & all_loci, float prop_eigenvalues){
    unsigned long total_loci = all_loci.size();
    for(int i = 0; i < total_loci ; i ++){
        all_loci[i].estimate_prior_variance(prop_eigenvalues);
    }
}

//iterate over all loci and estimate the posterior probablity via enumeration of max_number_of_cauals
double compute_all_posteriors_enumerate(vector<Locus>& all_loci, int max_number_of_cauals){
    size_t total_loci = all_loci.size();
    double total_log_bayes_factors = 0;
    for(int i = 0; i < total_loci; i++){
        total_log_bayes_factors += all_loci[i].enumerate_posterior(max_number_of_cauals);
    }

    return total_log_bayes_factors;
}

double compute_all_posteriors_mcmc(vector<Locus>& all_loci, size_t burn_in, size_t max_samples, int num_chains){
    size_t total_loci = all_loci.size();
    double total_log_bayes_factors = 0;
    double curr_log_bayes;
    for(int i = 0; i < total_loci; i++){
       // curr_log_bayes = all_loci[i].mcmc_posterior(burn_in, max_samples);
        curr_log_bayes = all_loci[i].mcmc_posterior(burn_in, max_samples, num_chains);
        total_log_bayes_factors += curr_log_bayes;
    }

    return total_log_bayes_factors;
}

//restructure the data and to create apporpriate input for optimization routine
void create_optimization_input(vector<Locus> & all_loci, optimization_input& opt_in){
    unsigned long total_loci = all_loci.size();
    VectorXd current_posterior;
    MatrixXf current_annotation_matrix;
    unsigned long current_num_snps;
    for(int i = 0; i < total_loci; i ++){
        current_posterior = all_loci[i].get_posterior_prob();
        current_annotation_matrix = all_loci[i].get_annotations();
        current_num_snps = all_loci[i].get_num_snps();
        unsigned long num_annotations = current_annotation_matrix.cols();
        for(int j = 0; j < current_num_snps; j++){
            opt_in.current_posterior.push_back(current_posterior(j));
            vector<double> snp_annotation;
            for(int k = 0; k < num_annotations; k++){
                snp_annotation.push_back(current_annotation_matrix(j,k));
            }
            opt_in.all_annotations.push_back(snp_annotation);
        }
    }
}

void update_all_gammas(vector<Locus>& all_loci, vector<double>& updated_gammas){
    unsigned long total_loci = all_loci.size();
    unsigned long total_annotations = updated_gammas.size();
    VectorXf gammas(total_annotations);
    for(int i = 0; i < total_annotations; i++){
        gammas(i) = updated_gammas[i];
    }
    for(int i = 0; i < total_loci; i++){
        all_loci[i].set_gamma(gammas);
    }
}

double run_EM(vector<Locus> & all_loci, execution_parameters& parameters ){

    int total_iterations = 1;
    double total_log_bayes_factor;
    double new_log_bayes_factor = 0;
    VectorXf current_gammas;
    float epsilon = 0.001;
    if(parameters.mcmc_flag == 1){
        initialize_generators(all_loci, parameters.sampling_seed);
    }

    while(total_iterations <= parameters.maxIter){
        total_iterations++;
        current_gammas = all_loci[0].get_gamma();
        cout << "Enrichment estimates at iteration " << total_iterations-1 << " is :" << endl << current_gammas << endl;
        //E Step,compute posterior probabilities for each SNP to be causal
        if(parameters.mcmc_flag == 1){
            epsilon = 0.1;
            total_log_bayes_factor = compute_all_posteriors_mcmc(all_loci, parameters.burn_in, parameters.max_samples, parameters.num_chains);
        }
        else{
            total_log_bayes_factor = compute_all_posteriors_enumerate(all_loci, parameters.max_causal);
        }

        cout << "Average log bayes factor: " << total_log_bayes_factor << endl << endl;

        //M step, run optimization for current iterate
        optimization_input opt_in;
        create_optimization_input(all_loci, opt_in);
        vector<double> current_gamma_vector = eigen_to_vector(current_gammas);
        void *optimization_input_ptr = &opt_in;

        //sometimes nlopt will throw error, safely catch and exit at current iterate
        try {
            optimize_complete_data_log_liklei(current_gamma_vector,-20,20, optimization_input_ptr);
            update_all_gammas(all_loci, current_gamma_vector);
            if(abs(total_log_bayes_factor - new_log_bayes_factor) < epsilon){
                break;
            }
            else{
                new_log_bayes_factor = total_log_bayes_factor;
            }

        }catch (exception& e) {
            //if an optimization error occurs, catch exception and exit program, modifying suffix to let user know
            //that error has occured but still write output from previous iteration.
            cout << "Optimization error encountered. Program exiting and writing output for current iteration " << endl;
            parameters.results_suffix = parameters.results_suffix + ".failed";
            parameters.likeli_name = parameters.results_suffix + ".failed";
            parameters.gammaName = parameters.gammaName + ".failed";
            break;
        }

    }
    return total_log_bayes_factor;
}
