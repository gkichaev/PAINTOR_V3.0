//
//  main.cpp
//  TestTarget
//
//  Created by Gleb Kichaev on 4/16/15.
//  Copyright (c) 2015 Gleb Kichaev. All rights reserved.
//

#include "Functions_Optimize.h"
#include "Functions_CoreModel.h"
#include "Functions_IO.h"

using namespace Eigen;
using namespace std;
//using namespace nlopt;

void Welcome_Message(){
    cout << "Welcome to PAINTOR v3.0. Copyright Gleb Kichaev 2015. This software is free for academic use. Please e-mail gkichaev@ucla.edu for questions." << endl;
    cout << "For questions or bug reports please contact: gkichaev@ucla.edu \n \n "<< endl;
    cout << "For required files and format specifications see User Manual \n \n" << endl;
    cout << "Usage: PAINTOR -input [input_filename] -in [input_directory] -out [output_directory] -Zhead [Zscore_header(s)] -LDname [LD_suffix(es)]  -annotations [annot_name1,annot_name2,...]  <other options> \n"<< endl;
    cout << "OPTIONS: -flag \t Description [default setting]  \n" << endl;
    cout << "-input \t (required) Filename of the input file containing the list of the fine-mapping loci [default: input.files]" << endl;
    cout << "-c \t The number of causal variants to consider per locus [default: 2]" << endl;
    cout << "-Zhead \t (required) The name(s) of the Zscores in the header of the locus file (comma separated) [default: N/A]" << endl;
    cout << "-LDname \t (required) Suffix(es) for LD files. Must match the order of Z-scores in locus file (comma separated) [Default:N/A]" << endl;
    cout << "-annotations \t The names of the annotations to include in model (comma separted) [default: N/A]" << endl;
    cout << "-in \t Input directory with all run files [default: ./ ]" << endl;
    cout << "-out \t Output directory where output will be written [default: ./ ]" << endl;
    cout << "-Gname \t Output Filename for enrichment estimates [default: Enrichment.Estimate]" << endl;
    cout << "-Lname \t Output Filename for log likelihood [default: Log.Likelihood]" << endl;
    cout << "-RESname \t Suffix for ouput files of results [Default: results] " << endl;
    cout << "-ANname \t Suffix for annotation files [Default: annotations]" << endl;
    cout << "-MI \t Maximum iterations for algorithm to run [Default: 10]" << endl;
    cout << "-post1CV \t fast conversion of Z-scores to posterior probabilites assuming a single casual variant and no annotations [Default: False]" << endl;
    cout << "-GAMinitial \t inititalize the enrichment parameters to a pre-specified value (comma separated) [Default: 0,...,0]" << endl;
    cout << "-variance \t specify prior variance on effect sizes scaled by sample size [Default: 25]" << endl ;
    cout << "-num_samples  \t specify number of samples to draw for each locus [Default: 50000]" << endl ;
    cout << "-prob_add \t specify geometric probablity to add a causal [Default: 0.25]" << endl ;
    cout << "-enumerate\t specify this flag if you want to enumerate all possible configurations followed by the max number of causal SNPs (eg. -enumerate 3 considers up to 3 causals at each locus) [Default: not specified]" << endl;
    cout << "-set_seed\t specify an integer as a seed for random number generator [default: clock time at execution]" << endl;
    cout << endl << endl ;
}

void Check_Flags(){

}

int main(int argc, const char * argv[])
{

    string in_dir = "./";
    string input_files = "input.files";
    string out_dir = "./";
    vector<string> annot_names;
    int max_causal = 2;
    string gammaName = "Enrichment.Values";
    string likeli_name= "Log.Likelihood";
    string single_post_flag;
    int maxIter= 10;
    string LD_suffix = "ld";
    string annot_suffix = "annotations";
    string results_suffix = "results";
    vector <string> LD_all_names;
    vector<string> z_headers;
    string ncp_flag= "default";
    int num_samples = 50000;
    double prob_add = .25;
    single_post_flag = "False";
    int sampling_seed = time(NULL);
    vector<VectorXd> z_score_loc;
    vector<string> snp_info;
    MatrixXd chol_factor;
    MatrixXd chol_factor_inv;
    MatrixXd out_annotations;
    vector<string> annot_header;
    string file_list_name = "input.files";
    string gamma_initial;
    double prior_variance = 25;
    int enumerate_flag = 0;
    if(argc < 2){
        Welcome_Message();
        return 0;
    }
    for(int i = 1; i < argc; i++){
        string argComp = argv[i];


        if(argComp.compare("-input")== 0){
            input_files = argv[i+1];
        }

        else if(argComp.compare("-in") == 0){
            in_dir = argv[i+1];
            string last_char = string(&in_dir[in_dir.size()-1]);
            if(last_char.compare("/")!=0){
                in_dir = in_dir + "/";
            }
        }

        else if(argComp.compare("-out") == 0){
            out_dir = argv[i+1];
            string last_char = string(&out_dir[out_dir.size()-1]);
            if(last_char.compare("/")!=0){
                out_dir = out_dir + "/";
            }
        }

        else if(argComp.compare("-c") == 0){
            max_causal = stoi(argv[i+1]);
        }


        else if(argComp.compare("-LDname") == 0){
            string temp_ldnames = argv[i+1];
            size_t n = count(temp_ldnames.begin(), temp_ldnames.end(), ',');
            if(n >0) {
                LD_all_names = split(temp_ldnames, ',');
            }
            else{
                LD_all_names = {temp_ldnames};
            }
        }

        else if(argComp.compare("-Zhead") == 0){
            string header_temp = argv[i+1];
            size_t n = count(header_temp.begin(), header_temp.end(), ',');
            if(n>0){
                z_headers = split(header_temp, ',');
            }
            else{
                z_headers = {header_temp};
            }
        }

        else if(argComp.compare("-annotations") == 0){
            string input_annotations = argv[i+1];
            size_t n = count(input_annotations.begin(), input_annotations.end(), ',');
            if(n >0){
                annot_names = split(input_annotations, ',');
            }
            else{
                annot_names= {input_annotations};
            }

        }

        else if(argComp.compare("-Gname") == 0){
            gammaName = argv[i+1];
        }

        else if(argComp.compare("-Lname") == 0){
            likeli_name = argv[i+1];
        }

        else if(argComp.compare("-post1CV") == 0){
            single_post_flag = argv[i+1];
        }

        else if(argComp.compare("-MI") == 0){
            maxIter = stoi(argv[i+1]);
        }

        else if(argComp.compare("-ANname") == 0){
            annot_suffix = argv[i+1];
        }

        else if(argComp.compare("-RESname") == 0){
            results_suffix = argv[i+1];
        }

        else if(argComp.compare("-GAMinitial") == 0){
            gamma_initial = argv[i+1];
        }

        else if(argComp.compare("-variance")==0){
            prior_variance = stod(argv[i+1]);

        }
        else if(argComp.compare("-num_samples")==0){
            num_samples = stoi(argv[i+1]);

        }
        else if(argComp.compare("-prob_add")==0){
            prob_add = stod(argv[i+1]);
        }
        else if(argComp.compare("-enumerate")==0){
            max_causal = stoi(argv[i+1]);
            enumerate_flag = 1;
        }
        else if(argComp.compare("-set_seed")==0){
            sampling_seed = stoi(argv[i+1]);
        }
    }

    /* initialize PAINTOR model parameters */

    vector<vector<VectorXd>> all_transformed_statistics;
    vector<vector<MatrixXd>> all_sigmas;
    vector<vector<string>> all_snp_info;
    vector<MatrixXd> all_annotations;
    vector<string> all_headers;
    Get_All_Input(input_files, in_dir, z_headers, annot_names, all_transformed_statistics, all_sigmas, all_annotations, LD_all_names, annot_suffix, all_snp_info, all_headers);

    vector<VectorXd> run_zscores;
    vector<MatrixXd> run_ld;
   // Reshape_Input(all_transformed_statistics, all_sigmas, run_zscores, run_ld);

    CausalProbs runProbs;
    VectorXd gamma_estimates(all_annotations[0].cols());
    gamma_estimates.setZero();
    gamma_estimates[0] = Get_Gamma_Zero(all_transformed_statistics);
    if(gamma_initial.size() > 0){
        vector<string> gamma_initial_split = split(gamma_initial, ',');
        if(gamma_initial_split.size() != annot_names.size()+1){
            cout << "Warning: Incorrect number of Enrichment parameters specified. Pre-setting all paramters to zero" << endl;
        }
        else{
            for(unsigned int i =0; i < gamma_initial_split.size(); i++){
                gamma_estimates[i]=stod(gamma_initial_split[i]);
            }
        }
    }

    /* run PAINTOR */

    if (single_post_flag.compare("True")==0) {
        if (z_headers.size() > 1) {
            cout << "Single Posterior Conversion not valid for more than 1 set of Z-scores per locus" << endl;
            return 0;
        }
        else{
            vector<VectorXd> all_single_post;
            VectorXd single_post;
            for(unsigned int i =0; i < all_transformed_statistics.size(); i++){
                single_post = Zscores2Post(all_transformed_statistics[i][0]);
                all_single_post.push_back(single_post);
            }
            Write_All_Output(input_files, out_dir, results_suffix, runProbs, all_snp_info, gamma_estimates,gammaName, 0, likeli_name, all_headers, annot_names);
        }
    }

    else{
        double Final_loglikeli;
        if(enumerate_flag == 1) {
            cout << "Running PAINTOR with full enumeration of a maximum of " << max_causal << " causals per locus" << endl;
            Final_loglikeli = EM_Run(runProbs, maxIter , all_transformed_statistics, gamma_estimates, all_annotations, all_sigmas, prior_variance, max_causal);

        }
        else{
            cout << "Running PAINTOR in sampling mode with seed value: "<< sampling_seed << endl;
            Final_loglikeli = EM_Run(runProbs, maxIter , all_transformed_statistics, gamma_estimates, all_annotations, all_sigmas, prior_variance, prob_add, num_samples, sampling_seed );
        }

        Write_All_Output(input_files, out_dir, results_suffix, runProbs, all_snp_info, gamma_estimates,gammaName, Final_loglikeli, likeli_name, all_headers, annot_names);
    }

    return 0;
}
