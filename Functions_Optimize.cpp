//
// Created by Gleb Kichaev on 10/21/15.
//

#include "Functions_Optimize.h"

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

/*
double dot_prod(vector<double>& vec1, vector<double>& vec2){
    double runsum = 0;
    for(unsigned i = 0; i < vec1.size(); i++){
        runsum += vec1[i]*vec2[i];
    }
    return(runsum);
}

vector<double> scalar_product(const double& scalar, vector<double>& vec){
    vector<double> scalvec(vec.size(),0);
    for(unsigned i = 0; i < vec.size(); i++){
        scalvec[i] = vec[i]*scalar;
    }
    return(scalvec);
}

vector<double> Vector_Sum(vector<double>& vec1, vector<double>& vec2){
    vector<double> outsum;
    for(unsigned i = 0; i< vec2.size(); i++){
        outsum.push_back(vec1[i] +vec2[i]);
    }
    return(outsum);
}

vector<double> GradientFxn(vector<double>& betas, ObjectiveData *in_data){
    ObjectiveData f_data;
    f_data.probs = in_data->probs;
    f_data.Aijs = in_data->Aijs;
    int numsnps = f_data.probs.size();
    double cij1 = 0;
    double cij0 = 0;
    double dp  = 0;
    vector<double> aij1(numsnps, 0);
    vector<double> aij0(numsnps, 0);
    vector<double> aij1_0(numsnps,0);
    vector<double> aij_out(numsnps,0);

    for(int i = 0; i < numsnps; i ++){
        dp = dot_prod(betas, f_data.Aijs[i]);
        cij1 = f_data.probs[i]*1/(1+exp(-1*dp));
        cij0 = (1-f_data.probs[i])*1/(1+exp(dp));
        aij1= scalar_product(-1*cij1, f_data.Aijs[i]);
        aij0 = scalar_product(cij0, f_data.Aijs[i]);
        aij1_0 = Vector_Sum(aij1, aij0);
        aij_out = Vector_Sum(aij_out, aij1_0);
    }

    return(aij_out);
}

double ObjectiveFxn(const vector<double> &x, vector<double> &grad, void *data){

    ObjectiveData *f_data = static_cast<ObjectiveData*>(data);
    int numsnps = f_data -> probs.size();
    double marginal_i;
    double runsum = 0;
    double cij1 = 0;
    double cij0 = 0;
    double dp  = 0;
    vector<double> temp;
    vector<double> betas = x;
    grad = GradientFxn(betas, f_data);

    for(int i = 0; i < numsnps; i++){
        temp = f_data -> Aijs[i];
        marginal_i = f_data -> probs[i];
        dp = dot_prod(betas, temp);
        cij1 = marginal_i*log(1+exp(dp));
        cij0 = (1-marginal_i)*log(1+exp(-dp));
        runsum = runsum - cij1 - cij0;

    }

    return(runsum);
}

void Optimize_Nlopt(vector<double>& x, double lower, double upper, double betaZero,  void* in_data){

    opt opt(LD_LBFGS, x.size());
    vector<double> lb;
    vector<double> ub;

    for(unsigned i = 0 ; i < x.size(); i++){
        lb.push_back(lower);
        ub.push_back(upper);
    }
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_max_objective(ObjectiveFxn, in_data);
    double minB;
    opt.optimize(x, minB);
}
*/

double Compute_Objective(VectorXd& current_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations){
    unsigned int total_snps = stacked_probabilites.size();
    double exponential_term;
    double objective_value=0;
    double value_c1 = 0;
    double value_c0 = 0;
    for(int i =0; i < total_snps; i++){
        exponential_term = current_gammas.dot(stacked_annotations.row(i));
        value_c1 = stacked_probabilites[i]*log(1+exp(exponential_term));
        value_c0 = (1-stacked_probabilites[i])*log(1+exp(-1*exponential_term));
        objective_value = objective_value - value_c1 -value_c0;
    }
    return objective_value;
}

void Compute_Gradient(VectorXd& current_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, VectorXd& gradient){
    unsigned int total_snps = stacked_probabilites.size();
    double exponential_term;
    double scaling_factor_c1;
    double scaling_factor_c0;
    gradient.setZero();
    for(int i =0; i < total_snps; i++){
        exponential_term = current_gammas.dot(stacked_annotations.row(i));
        scaling_factor_c1 = stacked_probabilites[i]*(1/(1+exp(-1*exponential_term)));
        scaling_factor_c0 = (1-stacked_probabilites[i])*(1/(1+exp(exponential_term)));
        gradient = gradient - scaling_factor_c1*stacked_annotations.row(i).transpose()+ scaling_factor_c0*stacked_annotations.row(i).transpose();
    }
}

void Compute_Hessian(VectorXd& current_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, MatrixXd& hessian){
    unsigned int total_snps = stacked_probabilites.size();
    double exponential_term;
    double scaling_factor_c1;
    double scaling_factor_c0;
    hessian.setZero();
    MatrixXd outer_product;
    for(int i = 0; i < total_snps; i++){
        outer_product = stacked_annotations.row(i).transpose()*stacked_annotations.row(i);
        exponential_term = current_gammas.dot(stacked_annotations.row(i));
        scaling_factor_c1 =  stacked_probabilites[i]*(exp(-1*exponential_term))/((1+exp(-1*exponential_term))*(1+exp(-1*exponential_term)));
        scaling_factor_c0 = (1-stacked_probabilites[i])*(exp(exponential_term))/((1+exp(exponential_term))*(1+exp(exponential_term)));
        hessian = hessian - scaling_factor_c1*outer_product - scaling_factor_c0*outer_product;
    }
}

int Gradient_Ascent(VectorXd& current_gammas, VectorXd& new_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, double gradient_tolerance, int max_iterations){
    unsigned int num_parameters = current_gammas.size();
    VectorXd gradient(num_parameters);
    VectorXd gamma_iterate(num_parameters);
    int max_line_search = 500;
    int search_counter;
    double tuner;
    double new_objective;
    double current_objective = Compute_Objective(current_gammas, stacked_probabilites, stacked_annotations);
    for(int i = 0; i < max_iterations;i++){
        tuner = 0.9;
        Compute_Gradient(current_gammas, stacked_probabilites,stacked_annotations, gradient);
        new_gammas = current_gammas + tuner * gradient;
        new_objective = Compute_Objective(new_gammas, stacked_probabilites, stacked_annotations);
        if(gradient.norm()>gradient_tolerance){
            search_counter=0;
            while (new_objective - current_objective <= 0){
                search_counter++;
                if(search_counter > max_line_search){
                    return -9;
                }
                else {
                    tuner = tuner * .9;
                    new_gammas = current_gammas + tuner * gradient;
                    new_objective = Compute_Objective(new_gammas, stacked_probabilites, stacked_annotations);
                }
            }
            current_objective = new_objective;
            current_gammas = new_gammas;
        }
        else{
            break;
        }
    }
    return 0;
}

