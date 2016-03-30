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

void Gradient_Ascent(VectorXd& current_gammas, VectorXd& new_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, double gradient_tolerance, int max_iterations, VectorXd& return_values){
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
                    return_values[0] = -9;
                    return;
                }
                else {
                    tuner = tuner * .9;
                    new_gammas = current_gammas + tuner * gradient;
                    new_objective = Compute_Objective(new_gammas, stacked_probabilites, stacked_annotations);
                }
            }
            current_objective = new_objective;
            return_values[1] = new_objective;
            current_gammas = new_gammas;
        }
        else{
            break;
        }
    }
    return_values[0] = 0;
}
