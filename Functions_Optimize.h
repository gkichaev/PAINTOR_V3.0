//
// Created by Gleb Kichaev on 10/21/15.
//

#ifndef PAINTOR_3_0_FUNCTIONS_OPTIMIZE_H
#define PAINTOR_3_0_FUNCTIONS_OPTIMIZE_H

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
#include <random>
#include "Functions_IO.h"
using namespace Eigen;
using namespace std;
//using namespace nlopt;




vector<double> GradientFxn(vector<double>& betas, ObjectiveData *in_data);
double ObjectiveFxn(const vector<double> &x, vector<double> &grad, void *data);
void Optimize_Nlopt(vector<double>& x, double lower, double upper, double betaZero,  void* in_data);
double dot_prod(vector<double>& vec1, vector<double>& vec2);
vector<double> scalar_product(const double& scalar, vector<double>& vec);
vector<double> Vector_Sum(vector<double>& vec1, vector<double>& vec2);

double Compute_Objective(VectorXd& current_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations);
void Compute_Gradient(VectorXd& current_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, VectorXd& gradient);
void Compute_Hessian(VectorXd& current_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, MatrixXd& hessian);
int Gradient_Ascent(VectorXd& current_gammas, VectorXd& new_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, double gradient_tolerance, int max_iterations);

#endif //PAINTOR_3_0_FUNCTIONS_OPTIMIZE_H
