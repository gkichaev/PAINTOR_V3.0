//
// Created by Gleb Kichaev on 10/21/15.
//


#include <unordered_map>
#include "Functions_CoreModel.h"
#include "Functions_Optimize.h"
#include "Functions_IO.h"

void BuildCausalVector(VectorXd& vec2Build , VectorXd& index){
    vec2Build.setZero();
    for(int i = 0; i <index.size(); i++){
        if(index[i] >0){
            vec2Build[index[i]-1] = 1;
        }
    }
}

void BuildCausalVector(VectorXd& vec2Build , VectorXd& index, vector<int>& index_as_vector){
    vec2Build.setZero();
    for(int i = 0; i <index.size(); i++){
        if(index[i] >0){
            vec2Build[index[i]-1] = 1;
            index_as_vector.push_back(index[i]-1);
        }
    }
}

inline double LogSum(double val1, double val2){
    double logsum = 0;
    if(val1 > val2){
        logsum = log(1 + exp(val2-val1)) + val1;
    }
    else{
        logsum = log(1 + exp(val1-val2)) + val2;
    }

    return(logsum);
}

int NextCombo(VectorXd& c, int k, int n) {
    for (int i= k; --i >= 0;) {
        if (++c[i] <= n-(k-i)) {
            while (++i < k)
                c[i]= c[i-1]+1;
            return 1;
        }
    }
    return 0;
}

double MVN_Density_ZeroMean(VectorXd& x, MatrixXd& sigma){
    MatrixXd b;
    b.setIdentity(sigma.rows(),sigma.cols());
    MatrixXd sigma_inv = sigma.ldlt().solve(b);
    double exponential_term = x.transpose()*(sigma_inv*x);
    double log_det = -.5*log(sigma.determinant());
    return log_det + -.5*exponential_term;
}

double MVN_Density_ZeroMean_print (VectorXd& x, MatrixXd& sigma){
    MatrixXd b;
    b.setIdentity(sigma.rows(),sigma.cols());
    MatrixXd sigma_inv = sigma.ldlt().solve(b);
    double exponential_term = x.transpose()*(sigma_inv*x);
    double log_det = -.5*log(sigma.determinant());
    cout << "det, exponential" << endl;
    cout << log_det << " " << -.5*exponential_term << endl;
    cout << "sigma " << endl;
    cout << sigma << endl;
    return log_det + -.5*exponential_term;
}
/*
double CalculateLogBayesFactor(VectorXd& zscore, MatrixXd& ld_matrix, VectorXd& causal_config, double prior_variance){
    MatrixXd ld_c(causal_config.size(), causal_config.size());
    VectorXd zscores_c(causal_config.size());
    for(int i =0; i < causal_config.size(); i++){
        zscores_c(i)=zscore(causal_config[i]);
        for(int j = 0; j < causal_config.size(); j++){
            ld_c(i,j)= ld_matrix(causal_config[i], causal_config[j]);
        }
    }

    MatrixXd offset = ld_c+prior_variance*(ld_c*ld_c);

    double bf_num = MVN_Density_ZeroMean(zscores_c, offset);
    double bf_denom = MVN_Density_ZeroMean(zscores_c, ld_c);

    return bf_num-bf_denom;
}
*/
int Sample_Int(mt19937 & generator ,vector<double>& weights){
    /* Function that is analagous to sample.int function in R */
    discrete_distribution<int> distribution(weights.begin(), weights.end());
    int sampled_int = distribution(generator);
    return sampled_int;
}

double Get_Rand_Unif( mt19937& generator){
    uniform_real_distribution<double> distribution(0.0,1.0);
    double rand_uniform = distribution(generator);
    return rand_uniform;
}

int Next_Move(int set_size, double move_probability, mt19937& generator){
// #return 0 if swap,  1 if add, 2 if delete
    double rand_uniform= Get_Rand_Unif(generator);
    int next_move;
    int power = set_size-1;
    if(set_size ==2){
        if(rand_uniform < move_probability){
            next_move=1;
        }
        else {
            next_move = 0;
        }
    }
    else{
        if(rand_uniform < pow(move_probability,power)){
            next_move = 1;
        }
   //     else if(rand_uniform < pow((1-move_probability),power)){
     //       next_move = 2;
      //  }
        else{
            next_move = 2;
           // next_move = 2;
        }
    }
    return(next_move);
}

int Next_Deletion(VectorXd & zscore_squared, vector<int>& current_set, mt19937& generator){
    //drop a snp from the causal set proportional to 1/zscore^2
    vector<double> probs(current_set.size());
    for(unsigned int i=0; i < current_set.size(); i++){
        probs[i] = 1/zscore_squared[current_set[i]];
    }
    int drop_from_current_set = Sample_Int(generator,probs);
    return drop_from_current_set;
}

int Next_Addition(VectorXd & zscore_squared, MatrixXd& LDsquared_complement, vector<int>& current_set, mt19937& generator){
    //get the next addition based on ld and z_score
    unsigned int N = zscore_squared.size();
    VectorXd ld_weight(N);
    ld_weight.setZero();
    //get max ld to the causal set
    for(unsigned int i =0; i < N; i++){
        for(unsigned int j=0 ; j<current_set.size(); j++){
            if(ld_weight[i]< LDsquared_complement(current_set[j], i)){
                ld_weight[i] = LDsquared_complement(current_set[j], i);
            }
        }
    }
    for(unsigned int j=0 ; j<current_set.size(); j++){
        ld_weight[current_set[j]]=0;
    }
    vector<double> probs(N);
    for(int i=0; i < N; i++){
        probs[i]=zscore_squared[i]*ld_weight[i];
    }
    int add_to_current_set =  Sample_Int(generator, probs);
    return add_to_current_set;
}

double Calculate_LogBayesFactor(VectorXd& zscore, MatrixXd& ld_matrix, vector<int>& causal_config, double prior_variance){
    unsigned int set_size = causal_config.size();
    MatrixXd ld_c(set_size, set_size);
    VectorXd zscores_c(set_size);
    for(unsigned int i =0; i < set_size; i++){
        zscores_c(i)=zscore(causal_config[i]);
        for(unsigned int j = 0; j < set_size; j++){
            ld_c(i,j)= ld_matrix(causal_config[i], causal_config[j]);
        }
    }
    MatrixXd offset = ld_c+prior_variance*(ld_c*ld_c);
    double bf_num = MVN_Density_ZeroMean(zscores_c, offset);
    double bf_denom = MVN_Density_ZeroMean(zscores_c, ld_c);
/*
    if(isnan(bf_num-bf_denom)){
      //  cout << "The offset matrix is:" << endl;
      //  cout <<offset << endl;
        cout << "num, denom" << endl;
        cout << bf_num << " " << bf_denom << endl;
        cout << "numerator" << endl;
        MVN_Density_ZeroMean_print(zscores_c, offset);
        cout << "denom" << endl;
        MVN_Density_ZeroMean_print(zscores_c, ld_c);
    }
*/
    return bf_num-bf_denom;
}

double Calculate_LogPrior(VectorXd& per_snp_prior,VectorXd& gammas, vector<int> causal_set){
    int N = per_snp_prior.size();
    VectorXd causal_indicator(N);
    causal_indicator.setZero();
    for(unsigned int i=0; i < causal_set.size(); i++){
        causal_indicator[causal_set[i]] = 1;
    }
    double probs = 0;
    for(unsigned int i = 0; i < per_snp_prior.size(); i++){
        if(causal_indicator[i] == 1){
            probs += log(per_snp_prior[i]);
        }
        else{
            probs += log(1-per_snp_prior[i]);
        }
    }
    return(probs);
}

inline double Prior_Snp_Probabilty(VectorXd& beta , VectorXd& aij){
    double dotprod = beta.dot(aij);
    double prob = 1/(1+exp(dotprod));
    return(prob);
}

void Compute_SNP_Priors(MatrixXd& annotation_matrix,VectorXd& gammas, VectorXd& per_snp_prior){
    for(int i= 0; i < annotation_matrix.rows(); i ++){
        VectorXd snp_annotation = annotation_matrix.row(i);
        per_snp_prior[i] = Prior_Snp_Probabilty(gammas, snp_annotation);
    }
}

double Prior_CausalSet_Probabilty(VectorXd& priorJ,VectorXd& beta, VectorXd& C_vector){
    double probs = 0;
    for(int i = 0; i < C_vector.size(); i++){
        if(C_vector[i] == 1){
            probs += log(priorJ[i]);
        }
        else{
            probs += log(1-priorJ[i]);
        }
    }
    return(probs);
}

void CasualSet_To_String(string& causal_string, vector<int>& causal_set){
    causal_string = "";
    for(unsigned int i =0; i < causal_set.size(); i++){
        causal_string += to_string(causal_set[i]);
        causal_string += ",";
    }
    causal_string.pop_back();
}

void Locus_Sampler(VectorXd& marginal, VectorXd& zscores, VectorXd& gammas, MatrixXd& annotations,  MatrixXd& ld_matrix, double& fullLikeli, double prior_variance, double  move_probability, int num_samples){
    unsigned int num_snps = zscores.size();
    double log_runsum = -1e150;
    VectorXd per_snp_priors(num_snps);
    random_device rd;
    mt19937 generator(rd());
    generator.seed(time(NULL));

    Compute_SNP_Priors(annotations, gammas, per_snp_priors);

    unordered_map <string, double> causal_set_posteriors;
    causal_set_posteriors.reserve(num_samples);
    vector<int> causal_set;
    causal_set.push_back(0);
    double log_bayes_factor;
    double log_prior;
    //compute single causal sets
    for(unsigned int i = 0; i < num_snps; i++){
        causal_set[0] = i;
        log_bayes_factor = Calculate_LogBayesFactor(zscores, ld_matrix, causal_set, prior_variance);
        log_prior = Calculate_LogPrior(per_snp_priors, gammas, causal_set);
        marginal[i] = log_bayes_factor+log_prior;
        log_runsum = LogSum(log_runsum, log_bayes_factor+log_prior);
    }

    causal_set[0] = 0;
    causal_set.push_back(num_snps/2);
    MatrixXd ld_sq_complement(num_snps,num_snps);
    VectorXd Zscore_sq(num_snps);
    ld_sq_complement = 1-ld_matrix.array().square();
    Zscore_sq = zscores.array().square();
    int counter = 0;
    //sample multiple causal sets
    int move;
    int index_add;
    int index_delete;
    string causal_string;
    unordered_map<string,double>::const_iterator found;
    while(counter <= num_samples){
        counter++;
        //sample next causal set
        move = Next_Move(causal_set.size(), move_probability,generator);
        if(move == 0){
            index_delete = Next_Deletion(Zscore_sq,causal_set,generator);
            causal_set.erase(causal_set.begin() + index_delete);
            index_add = Next_Addition(Zscore_sq,ld_sq_complement, causal_set,generator);
            causal_set.push_back(index_add);
        }
        else if(move ==1){
            index_add = Next_Addition(Zscore_sq,ld_sq_complement, causal_set,generator);
            causal_set.push_back(index_add);
        }
        else{
            index_delete = Next_Deletion(Zscore_sq,causal_set,generator);
            causal_set.erase(causal_set.begin()+index_delete);
        }
        //check if causal set has already been evaluated and store in hash table if it has not been
        CasualSet_To_String(causal_string, causal_set);
        found = causal_set_posteriors.find (causal_string);
        if ( found == causal_set_posteriors.end()) {
            log_bayes_factor = Calculate_LogBayesFactor(zscores, ld_matrix, causal_set, prior_variance);
            log_prior = Calculate_LogPrior(per_snp_priors, gammas, causal_set);
            causal_set_posteriors[causal_string] = log_bayes_factor+log_prior;
        }
    }
    Marginalize_Sets(causal_set_posteriors, marginal, log_runsum);

    for(int i = 0; i < num_snps; i ++){
        marginal[i] = marginal[i]-log_runsum;
    }

    fullLikeli = fullLikeli+log_runsum;
}

void Locus_Sampler_Multi(VectorXd& marginal, vector<VectorXd>& zscores, VectorXd& gammas, MatrixXd& annotations,  vector<MatrixXd>& ld_matrix, double& fullLikeli, double prior_variance, double  move_probability, int num_samples){
    unsigned int num_snps = zscores[0].size();
    unsigned int num_sets = zscores.size();
    double log_runsum = -1e150;
    VectorXd per_snp_priors(num_snps);
    random_device rd;
    mt19937 generator(rd());
    generator.seed(time(NULL));

    Compute_SNP_Priors(annotations, gammas, per_snp_priors);
    unordered_map <string, double> causal_set_posteriors;
    causal_set_posteriors.reserve(num_samples);
    vector<int> causal_set;
    causal_set.push_back(0);
    double log_bayes_factor;
    double log_prior;
    //compute single causal sets
    for(unsigned int i = 0; i < num_snps; i++){
        causal_set[0] = i;
        log_prior = Calculate_LogPrior(per_snp_priors, gammas, causal_set);
        log_bayes_factor=0;
        for(unsigned int j=0; j < num_sets; j++) {
            log_bayes_factor += Calculate_LogBayesFactor(zscores[j], ld_matrix[j], causal_set, prior_variance);
        }
        marginal[i] = log_bayes_factor + log_prior;
        log_runsum = LogSum(log_runsum, log_bayes_factor + log_prior);
    }

    causal_set[0] = 0;
    causal_set.push_back(num_snps/2);
    MatrixXd ld_sq_complement(num_snps,num_snps);
    VectorXd Zscore_sq(num_snps);
    VectorXd divisor(num_snps);
    divisor.setZero();
    divisor = divisor.array()+(num_sets);
    Zscore_sq.setZero();
    ld_sq_complement = 1-ld_matrix[0].array().square();
    for(unsigned int j=0; j < num_sets; j++) {
        Zscore_sq = Zscore_sq.array() + zscores[j].array().square();
    }
    Zscore_sq = Zscore_sq.array()/divisor.array();
    //cout << Zscore_sq << endl;
    int counter = 0;
    //sample multiple causal sets
    int move;
    int index_add;
    int index_delete;
    string causal_string;
    unordered_map<string,double>::const_iterator found;
    while(counter <= num_samples){
        counter ++;
        //sample next causal set
        move = Next_Move(causal_set.size(), move_probability,generator);
        if(move == 0){
            index_delete = Next_Deletion(Zscore_sq,causal_set,generator);
            causal_set.erase(causal_set.begin() + index_delete);
            index_add = Next_Addition(Zscore_sq,ld_sq_complement, causal_set,generator);
            causal_set.push_back(index_add);
        }
        else if(move ==1){
            index_add = Next_Addition(Zscore_sq,ld_sq_complement, causal_set,generator);
            causal_set.push_back(index_add);
        }
        else{
            index_delete = Next_Deletion(Zscore_sq,causal_set,generator);
            causal_set.erase(causal_set.begin()+index_delete);
        }
        //check if causal set has already been evaluated and store in hash table if it has not been
        CasualSet_To_String(causal_string, causal_set);
        found = causal_set_posteriors.find (causal_string);
        if (found == causal_set_posteriors.end()) {
            log_bayes_factor=0;
            for(unsigned int j=0; j < num_sets; j++) {
                log_bayes_factor += Calculate_LogBayesFactor(zscores[j], ld_matrix[j], causal_set, prior_variance);
            }
            log_prior = Calculate_LogPrior(per_snp_priors, gammas, causal_set);
            causal_set_posteriors[causal_string] = log_bayes_factor+log_prior;
        }
    }
    Marginalize_Sets(causal_set_posteriors, marginal, log_runsum);
    for(int i = 0; i < num_snps; i ++){
        marginal[i] = marginal[i]-log_runsum;
    }
    fullLikeli = fullLikeli+log_runsum;
}

void string_to_int(string& input_string, vector<int>& output_int){
    vector<string> string_split = split(input_string, ',');
    for(int i = 0; i < string_split.size(); i++){
        output_int.push_back(stoi(string_split[i]));
    }
}

void Marginalize_Sets(unordered_map<string,double>& causal_posteriors, VectorXd& marginals, double& running_sum){
    string string_set;
    double posterior;
    int counter =0;
    vector<int> num_k_sampled;
    for(int k = 0; k < 6; k++){
        num_k_sampled.push_back(0);
    }
    for( auto it = causal_posteriors.begin(); it != causal_posteriors.end(); ++it){
        counter++;
        string_set =  it->first;
        posterior =  it->second;
        vector<int> causal_set;
        string_to_int(string_set, causal_set);

        int causal_set_size = causal_set.size();
        for(int j =0; j<causal_set_size; j++){
            marginals[causal_set[j]] = LogSum(marginals[causal_set[j]], posterior);
        }
        if(causal_set_size> 5){
            num_k_sampled[5]++;
        }
        else{
            num_k_sampled[causal_set_size]++;
        }
        running_sum = LogSum(running_sum, posterior);
    }

    for(int k = 2; k < 5; k++){
        cout << "Sampled " << num_k_sampled[k] << " unique multi-causal sets of size: " << k << endl;
    }
    cout << "Sampled " << num_k_sampled[5] << " unique multi-causal sets of size 5 or greater " << endl;
    cout << "Sampled " << counter << " unique multi-causal sets in total" << endl << endl;
}

void Enumerate_Posterior(VectorXd& Marginal, vector<VectorXd>& zscores, VectorXd& beta, MatrixXd& annotations,  vector<MatrixXd>& ld_matrix, int NC, double& fullLikeli, double prior_variance){
    int numsnps = zscores[0].size();
    int num_pops = zscores.size();
    double runsum = -1e150;
    for(int i =0 ; i < Marginal.size(); i++){
        Marginal(i) = -1e150;
    }
    double c_prob = 0;
    double sum = 0;
    VectorXd causConfig(numsnps);
    causConfig.setZero();
    VectorXd causal_index(NC);
    causal_index.setZero();
    VectorXd per_snp_prior(numsnps);
    Compute_SNP_Priors(annotations, beta, per_snp_prior);

    double log_bayes_factor = 0;
    double log_bayes_factor_pop;
    //runsum = sum;
    int counter = 1;
    while(NextCombo(causal_index, NC, numsnps+1) == 1){
        vector<int> causal_index_as_vector;
        BuildCausalVector(causConfig, causal_index,causal_index_as_vector);
        c_prob = Prior_CausalSet_Probabilty(per_snp_prior, beta, causConfig);
        log_bayes_factor = 0;
        for (int i = 0;  i<num_pops ; i++) {
            log_bayes_factor_pop=Calculate_LogBayesFactor(zscores[i], ld_matrix[i], causal_index_as_vector, prior_variance);
            if(std::isnan(log_bayes_factor_pop)){
                log_bayes_factor_pop = -1e150;
            }
            log_bayes_factor += log_bayes_factor_pop;
        }

        sum = c_prob+log_bayes_factor;
        runsum = LogSum(runsum, sum);

        for(int j = 0; j < causal_index.size(); j++){
            if(causal_index[j] >0){
                Marginal[causal_index[j]-1] = LogSum(Marginal[causal_index[j]-1], sum);
            }
        }
        counter ++;
    }
    for(int f = 0 ; f < Marginal.size(); f++){
        Marginal[f] = Marginal[f]- runsum;
    }
    cout << "Enumerated " << counter << " causal sets" << endl;
    //cout <<runsum << endl;
    fullLikeli = fullLikeli+runsum;
}

//E step with full enumeration
double Estep(vector<vector<VectorXd>> &Zscores, VectorXd &betas, vector<MatrixXd> &Aijs, vector<vector<MatrixXd>>& ld_matrices, CausalProbs &E_out, double prior_variance,int  max_causals){
    vector<VectorXd> marginal_i;
    VectorXd temp;
    VectorXd exp_temp;
    vector<double> stacker;
    vector<double> stack_temp;
    double fullLikeli = 0;
    for(unsigned int i = 0; i < Zscores.size(); i ++){
        VectorXd locus_log_marginals(Zscores[i][0].size());
        Enumerate_Posterior(locus_log_marginals,Zscores[i], betas, Aijs[i], ld_matrices[i], max_causals, fullLikeli, prior_variance);
        exp_temp = locus_log_marginals.array().exp();
        marginal_i.push_back(exp_temp);
        stack_temp =  eigen2vec(exp_temp);
        stacker.insert(stacker.end(), stack_temp.begin(), stack_temp.end());
    }
    E_out.probs_locs = marginal_i;
    E_out.probs_stacked = stacker;
    return(fullLikeli);
}


double Estep(vector<vector<VectorXd>> &Zscores, VectorXd &betas, vector<MatrixXd> &Aijs, vector<vector<MatrixXd>>& ld_matrices, CausalProbs &E_out, double prior_variance,double  move_probability, int num_samples){
    vector<VectorXd> marginal_i;
    VectorXd temp;
    VectorXd exp_temp;
    vector<double> stacker;
    vector<double> stack_temp;
    double fullLikeli = 0;
    for(unsigned int i = 0; i < Zscores.size(); i ++){
        VectorXd locus_log_marginals(Zscores[i][0].size());
        Locus_Sampler_Multi(locus_log_marginals, Zscores[i], betas, Aijs[i], ld_matrices[i], fullLikeli, prior_variance, move_probability, num_samples);
        exp_temp = locus_log_marginals.array().exp();
        marginal_i.push_back(exp_temp);
        stack_temp =  eigen2vec(exp_temp);
        stacker.insert(stacker.end(), stack_temp.begin(), stack_temp.end());
    }
    E_out.probs_locs = marginal_i;
    E_out.probs_stacked = stacker;

    return(fullLikeli);
}

double Estep(vector<VectorXd> &Zscores, VectorXd &betas, vector<MatrixXd> &Aijs, vector<MatrixXd>& ld_matrices, CausalProbs &E_out, double prior_variance,double  move_probability, int num_samples){
    vector<VectorXd> marginal_i;
    VectorXd temp;
    VectorXd exp_temp;
    vector<double> stacker;
    vector<double> stack_temp;
    double fullLikeli = 0;
    for(unsigned int i = 0; i < Zscores.size(); i ++){
        VectorXd locus_log_marginals(Zscores[i].size());
        Locus_Sampler(locus_log_marginals, Zscores[i], betas, Aijs[i], ld_matrices[i], fullLikeli, prior_variance, move_probability, num_samples);
        exp_temp = locus_log_marginals.array().exp();
        marginal_i.push_back(exp_temp);
        stack_temp =  eigen2vec(exp_temp);
        stacker.insert(stacker.end(), stack_temp.begin(), stack_temp.end());
    }
    E_out.probs_locs = marginal_i;
    E_out.probs_stacked = stacker;

    return(fullLikeli);
}

vector<double> eigen2vec(VectorXd &vec){
    vector<double> outvec;
    for(unsigned i = 0; i < vec.size(); i++){
        outvec.push_back(vec[i]);
    }
    return(outvec);
}

vector<vector<double> > Stack_Matrices(vector<MatrixXd> &mats){
    vector<vector<double> > out_stack;

    for(unsigned i = 0; i < mats.size(); i++){
        MatrixXd tempMat = mats[i];
        for(unsigned j = 0; j < tempMat.rows(); j++){
            VectorXd temp = tempMat.row(j);
            out_stack.push_back(eigen2vec(temp));
        }
    }
    return(out_stack);
}

void Stack_EigenMatrices(vector<MatrixXd> &mats, MatrixXd& stacked_matrices){
    int place_holder = 0;
    for(unsigned int i = 0; i < mats.size(); i++){
        int num_snps = mats[i].rows();
        stacked_matrices.middleRows(place_holder, num_snps) = mats[i];
        place_holder += num_snps;
    }
}

void Copy_CausalProbs(CausalProbs& source, CausalProbs& destination ){
    unsigned int num_loci = source.probs_locs.size();
    if(destination.probs_locs.size()==0){
        for(unsigned int i = 0; i < num_loci ; i++){
            destination.probs_locs.push_back(source.probs_locs[i]);
        }
        destination.probs_stacked= source.probs_stacked;
    }
    else {
        for (unsigned int i = 0; i < num_loci; i++) {
            destination.probs_locs[i] = source.probs_locs[i];
        }
        destination.probs_stacked = source.probs_stacked;
    }
}

double EM_Run(CausalProbs &probabilites, int iter_max, vector<vector<VectorXd>> &Zscores,  VectorXd &beta_int, vector<MatrixXd> &Aijs,vector<vector<MatrixXd>> &ld_matrix , double prior_variance,double  move_probability, int num_samples){
    //Run EM with sampling
    vector<double> beta_run = eigen2vec(beta_int);
    ObjectiveData Opt_in;
    Opt_in.Aijs = Stack_Matrices(Aijs);
    int iterations = 0;
    VectorXd beta_update;
    double likeliOld;
    double likeli =1;
    double gradient_tolerance = 1e-5;
    int max_descents = 200;
    int succesful_opt;
    CausalProbs curr_probs;
    while(iterations < iter_max) {
        likeli = Estep(Zscores, beta_int, Aijs, ld_matrix, probabilites, prior_variance, move_probability, num_samples);
        VectorXd stacked_probabilities = vector2eigen(probabilites.probs_stacked);
        MatrixXd stacked_annotations(stacked_probabilities.size(), beta_int.size());
        Stack_EigenMatrices(Aijs, stacked_annotations);
        if(likeliOld < likeli || iterations==0){
            Copy_CausalProbs(probabilites, curr_probs);
            VectorXd new_gammas(beta_int.size());
            succesful_opt= Gradient_Ascent(beta_int,  new_gammas, stacked_probabilities, stacked_annotations, gradient_tolerance, max_descents);
            likeliOld = likeli;
            if(succesful_opt == 0){
                beta_int = new_gammas;
                cout << std::setprecision(9) << "Sum of log Bayes Factors at iteration: " << iterations+1 << ": " << likeli << endl;
                cout << "Parameter Estimates:" << endl << new_gammas << endl << endl;
            }
            else{
                cout << "Warning: Optimization unstable at iteration: " << iterations+1 << endl;
                cout << "Maximum liklelihood estimate for enrichments potentially not reached" << endl;
                cout << "Program is exiting and returning values corresponding to estimates from iteration: " << iterations << endl;
                break;
            }
        }
        else{
            cout << "Last Iteration did not improve model fit. Resampling at prevoius estimates. " << endl;
        }
        iterations++;
    }
    Copy_CausalProbs(curr_probs,probabilites);
    return(likeliOld);
}

double EM_Run(CausalProbs &probabilites, int iter_max, vector<vector<VectorXd>> &Zscores,  VectorXd &beta_int, vector<MatrixXd> &Aijs,vector<vector<MatrixXd>> &ld_matrix , double prior_variance,int max_causals){
    //Run EM with full enumeration
    vector<double> beta_run = eigen2vec(beta_int);
    ObjectiveData Opt_in;
    Opt_in.Aijs = Stack_Matrices(Aijs);
    int iterations = 0;
    VectorXd beta_update;
    double likeliOld;
    double likeli =1;
    double gradient_tolerance = 1e-5;
    int max_descents = 200;
    int succesful_opt;
    CausalProbs curr_probs;
    while(iterations < iter_max) {
        likeli = Estep(Zscores, beta_int, Aijs, ld_matrix, probabilites, prior_variance, max_causals);
        VectorXd stacked_probabilities = vector2eigen(probabilites.probs_stacked);
        MatrixXd stacked_annotations(stacked_probabilities.size(), beta_int.size());
        Stack_EigenMatrices(Aijs, stacked_annotations);
        if(abs(likeliOld - likeli)> 0.01){
            Copy_CausalProbs(probabilites, curr_probs);
            VectorXd new_gammas(beta_int.size());
            succesful_opt= Gradient_Ascent(beta_int,  new_gammas, stacked_probabilities, stacked_annotations, gradient_tolerance, max_descents);
            likeliOld = likeli;
            if(succesful_opt == 0){
                beta_int = new_gammas;
                cout << std::setprecision(9) << "Sum of log Bayes Factors at iteration: " << iterations+1 << ": " << likeli << endl;
                cout << "Parameter Estimates:" << endl << new_gammas << endl << endl;
            }
            else{
                cout << "Warning: Optimization unstable at iteration: " << iterations+1 << endl;
                cout << "Maximum liklelihood estimate for enrichments potentially not reached" << endl;
                cout << "Program is exiting and returning values corresponding to estimates from iteration: " << iterations << endl;
                break;
            }
        }
        else {
            break;
        }
        iterations++;
    }
    Copy_CausalProbs(curr_probs,probabilites);
    return(likeliOld);
}

VectorXd Zscores2Post(VectorXd& Zs){
    VectorXd post(Zs.size());
    VectorXd Zsq  = Zs.array().square();
    for(int i = 0; i < Zsq.size(); i ++){
        VectorXd Ztemp = (Zsq.array() - Zsq[i])/2;
        VectorXd Zexp = Ztemp.array().exp();
        post[i] = 1/Zexp.sum();
    }
    return(post);
}

VectorXd vector2eigen(vector<double> &vec){
    VectorXd outVec(vec.size());
    for(unsigned i =0; i <vec.size(); i++){
        outVec[i] = vec[i];
    }
    return(outVec);
}


