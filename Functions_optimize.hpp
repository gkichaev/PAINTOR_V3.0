//
//  Functions_optimize.hpp
//  PAINTOR_V3.1
//
//  Created by Gleb Kichaev on 11/21/17.
//  Copyright © 2017 Gleb Kichaev. All rights reserved.
//

#ifndef Functions_optimize_hpp
#define Functions_optimize_hpp

#include <stdio.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen>

#include "nlopt.hpp"
#include "Functions_model.hpp"



using namespace Eigen;
using namespace std;
using namespace nlopt;


/*
struct optimization_input{
    vector<double> current_posterior;
    vector<vector<double> > all_annotations;
};
*/


double dot_prod(vector<double>& vec1, vector<double>& vec2);
vector<double> vector_scalar_product(const double& scalar, vector<double>& vec);
vector<double> vector_sum(vector<double>& vec1, vector<double>& vec2);
double complete_data_log_liklei(const vector<double> &x, vector<double> &grad, void *data);
vector<double> compute_gradient(vector<double>& gamma, optimization_input *in_data);
void optimize_complete_data_log_liklei(vector<double>& gamma_initial, double lower, double upper, void* in_data);

#endif /* Functions_optimize_hpp */
