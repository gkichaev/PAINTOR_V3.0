//
//  main.cpp
//  PAINTOR_V3.1
//
//  Created by Gleb Kichaev on 10/25/17.
//  Copyright © 2017 Gleb Kichaev. All rights reserved.
//

#include <iostream>
#include "Functions_IO.hpp"

int main(int argc, const char * argv[]) {
    
    
    if(argc < 2){
        Welcome_Message();
        return 0;
    }
    
    else{
        
        execution_parameters parameters;
        Parse_command_line(argc, argv, parameters);
        vector<Locus> all_loci;
        display_execution_message(parameters);
        Get_all_input(parameters.input_files, all_loci, parameters);
        all_prior_variance_esimation(all_loci, parameters.prop_ld_eigenvalues);
        size_t num_annotations = parameters.model_annotations.size()+1;
        if(parameters.gamma_initial.size() == 0){
            initialize_gamma_zero(all_loci, num_annotations);
        }
        else{
            initialize_gammas_input(all_loci, parameters);
        }

        double total_log_bayes_factor = run_EM(all_loci, parameters);
        write_all_posteriors(all_loci, parameters);
        write_estimated_prior_variance(all_loci, parameters);
        VectorXf final_gamma = all_loci[0].get_gamma();
        write_other_output(final_gamma, total_log_bayes_factor, parameters);
    }

    return 0;
}
