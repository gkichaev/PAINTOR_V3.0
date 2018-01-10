//
//  Functions_optimize.cpp
//  PAINTOR_V3.1
//
//  Created by Gleb Kichaev on 11/21/17.
//  Copyright © 2017 Gleb Kichaev. All rights reserved.
//

#include "Functions_optimize.hpp"

//return the dot prodcut of two vectors
double dot_prod(vector<double>& vec1, vector<double>& vec2){
    double runsum = 0;
    for(unsigned int i = 0; i < vec1.size(); i++){
        runsum += vec1[i]*vec2[i];
    }
    return(runsum);
}


//return scalar vector product
vector<double> vector_scalar_product(const double& scalar, vector<double>& vec){
    vector<double> scalvec(vec.size(),0);
    for(unsigned i = 0; i < vec.size(); i++){
        scalvec[i] = vec[i]*scalar;
    }
    return(scalvec);
}


//return standard vector sum
vector<double> vector_sum(vector<double>& vec1, vector<double>& vec2){
    vector<double> outsum;
    for(unsigned i = 0; i< vec2.size(); i++){
        outsum.push_back(vec1[i] +vec2[i]);
    }
    return(outsum);
}


vector<double> compute_gradient(vector<double>& gamma, optimization_input *in_data){
    
    unsigned long numsnps = in_data -> current_posterior.size();
    double cij1 = 0;
    double cij0 = 0;
    double dp  = 0;
    vector<double> aij1(numsnps, 0);
    vector<double> aij0(numsnps, 0);
    vector<double> aij1_0(numsnps,0);
    vector<double> aij_out(numsnps,0);
    
    for(int i = 0; i < numsnps; i ++){
        dp = dot_prod(gamma, in_data->all_annotations[i]);
        cij1 = in_data->current_posterior[i]*1/(1+exp(-1*dp));
        cij0 = (1-in_data->current_posterior[i])*1/(1+exp(dp));
        aij1= vector_scalar_product(-1*cij1, in_data->all_annotations[i]);
        aij0 = vector_scalar_product(cij0, in_data->all_annotations[i]);
        aij1_0 = vector_sum(aij1, aij0);
        aij_out = vector_sum(aij_out, aij1_0);
    }
    
    return(aij_out);
}



double complete_data_log_liklei(const vector<double> &x, vector<double> &grad, void *data){
    
    optimization_input *in_data = static_cast<optimization_input*>(data);
    int numsnps = in_data -> current_posterior.size();
    double marginal_i;
    double runsum = 0;
    double cij1 = 0;
    double cij0 = 0;
    double dp  = 0;
    vector<double> temp;
    vector<double> gammas = x;
    grad = compute_gradient(gammas, in_data);
    
    for(int i = 0; i < numsnps; i++){
        temp = in_data->all_annotations[i];
        marginal_i = in_data -> current_posterior[i];
        dp = dot_prod(gammas, temp);
        cij1 = marginal_i*log(1+exp(dp));
        cij0 = (1-marginal_i)*log(1+exp(-dp));
        runsum = runsum - cij1 - cij0;
    }
    
    return(runsum);
}


void optimize_complete_data_log_liklei(vector<double>& gamma_initial, double lower, double upper, void* in_data){
    
    opt opt(LD_LBFGS, gamma_initial.size());
    vector<double> lb;
    vector<double> ub;
    
    for(unsigned i = 0 ; i < gamma_initial.size(); i++){
        lb.push_back(lower);
        ub.push_back(upper);
    }
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_max_objective(complete_data_log_liklei, in_data);
    double minB;
    opt.optimize(gamma_initial, minB);
}

