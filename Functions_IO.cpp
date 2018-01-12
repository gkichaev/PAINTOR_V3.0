//
//  Functions_IO.cpp
//  PAINTOR_V3.1
//
//  Created by Gleb Kichaev on 10/25/17.
//  Copyright © 2017 Gleb Kichaev. All rights reserved.
//

#include "Functions_IO.hpp"

#include <sys/stat.h>




//split characters based on a given delimiter 'delim'
vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

//check to see if file name given exists
bool check_file_exists (const string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

//mandatory flags that need to be specified are (1) input file list,
//(2) zscore names in locus file header, (3) ld suffix
void Check_Mandatory_Flags(int argc, const char * argv[]){
    bool input_missing = true;
    bool in_dir_missing = true;
    bool ld_name_missing = true;
    bool zscore_flag_missing = true;
    for(int i = 0; i < argc; i++){
        string argComp = argv[i];
        if(argComp.compare("-input")== 0){
            input_missing = false;
        }
        else if(argComp.compare("-in") == 0){
            in_dir_missing = false;
        }
        
        else if(argComp.compare("-LDname") == 0){
            ld_name_missing = false;
        }
        else if(argComp.compare("-Zhead") == 0){
            zscore_flag_missing = false;
        }
    }
    
    if(input_missing){
        cout << "Error: Please specify -input flag" << endl;
        exit(0);
    }
    
    if(in_dir_missing){
        cout << "Error: Please specify -in flag" << endl;
        exit(0);
    }
    if(ld_name_missing){
        cout << "Error: Please specify -LDname flag" << endl;
        exit(0);
    }
    if(zscore_flag_missing){
        cout << "Error: Please specify -Zhead flag" << endl;
        exit(0);
    }
    
}


//Functions to read input locus
void Read_Locus( Locus &  current_locus, string & locus_name, execution_parameters & info){
    string full_path = info.in_dir+ locus_name;
    if(!check_file_exists(full_path)){
        cout << "Error: Locus file " << locus_name << " not found in directory : " << info.in_dir << endl;
        cout << "Terminating Execution" << endl;
        exit(0);
    }
    
    //determine the index of the Zscore column(s)
    string header;
    ifstream locus;
    locus.open(full_path);
    getline(locus, header);
    header.erase(remove(header.begin(), header.end(), '\r'), header.end()); //remove line terminator
    char delimiter = ' ';
    current_locus.add_info(header);
    vector<string> split_header = split(header, delimiter);
    vector<int> z_index;
    for(unsigned int i=0; i < split_header.size(); i++){
        for(unsigned int j=0; j < info.zscore_header.size(); j++){
            if(info.zscore_header[j].compare(split_header[i])==0){
                z_index.push_back(i);
            }
        }
    }
    //Terminate if header that is supplied with -Zhead flag is incorrect
    if(z_index.size() != info.zscore_header.size()){
        cout << "Error! Specified Z-score headers not found in locus file" << endl;
        cout << "Terminating Execution" << endl;
        exit(0);
    }
    
    //read the rest of the locus file
    string input_line;
    vector<string> split_line;
    vector<vector<string>> snp_info_split;
    while(getline(locus,input_line)){
        split_line = split(input_line, ' ');
        current_locus.add_info(input_line);
        snp_info_split.push_back(split_line);
    }
    locus.close();
    unsigned long num_snps = snp_info_split.size();
    VectorXf association_stat;
    for(unsigned int i =0; i < z_index.size();i++){
        association_stat.resize(num_snps,1);
        for(unsigned int j=0; j < num_snps; j++){
            if(snp_info_split[j][z_index[i]].compare("NA") == 0){
                cout <<  "Missing values not supported in PAINTOR 3.XX. Please remove SNP or use PAINTOR 2.XX" << endl;
                cout << "Terminating Execution. Check Locus file: " << locus_name << " at index (" << i << "," << j << ")" <<endl;
                exit(0);
            }
            else{
                association_stat[j] = stod(snp_info_split[j][z_index[i]]);
            }
        }
        current_locus.set_Zscore(association_stat);
        current_locus.set_num_snps(association_stat.size());
    }
}

//Read in a single LD matrix 
void Read_LD(string &input_directory, string& fname , MatrixXf & ld_matrix){
    
    //check if file exists
    if(!check_file_exists(input_directory+fname)){
        cout << "Error: LD file " << fname << " not found in directory : " << input_directory << endl;
        cout << "Terminating Execution" << endl;
        exit(0);
    }
    
    //read in ld file and store into matrix
    ifstream ld_file;
    ld_file.open(input_directory+fname);
    string ld_line;
    vector<string> ld_line_split;
    vector<vector<string>> all_ld_line_split;
    while(getline(ld_file, ld_line)){
        ld_line_split = split(ld_line, ' ');
        all_ld_line_split.push_back(ld_line_split);
    }
    ld_file.close();
    unsigned long num_snps = all_ld_line_split.size();
    string val_str;
    int check_ind;
    MatrixXf out_matrix(num_snps,num_snps);
    float epsilon = 1e-5;
    for(unsigned int i = 0; i < num_snps; i++){
        for(unsigned int j = 0; j < num_snps; j++){
            //check to make sure last digit is a number, otherwise throw error
            val_str = all_ld_line_split[i][j];
            check_ind = val_str.length()-1;
            if(isdigit(val_str[check_ind])) {
                out_matrix(i, j) = stod(val_str);
            }
            else{
                cout <<  "Error: Infinite/non-numeric value in LD matrix detected!!" << endl;
                cout << "Terminating Execution. Check ld matrix file: " << fname << " at index (" << i << "," << j << ")" <<endl;
                exit(0);
            }
        }
        out_matrix(i,i) = 1+epsilon;
    }
    ld_matrix = out_matrix;
}

//Read in the annotation matrix for locus
void Read_Annotations(string &input_directory, string& fname , vector<string>& model_annotations, MatrixXf & out_annotations){
    
    
    if(!check_file_exists(input_directory+fname)){
        cout << "Error: Annotation file " << fname << " not found in directory : " << input_directory << endl;
        cout << "Software requires annotation file to execute" <<endl;
        cout << "Terminating Execution" << endl;
        exit(0);
    }
    
    //load in all annotations
    ifstream annotation_file;
    annotation_file.open(input_directory+fname);
    string annotation_header;
    getline(annotation_file,annotation_header);
    vector<string> annotation_header_split;
    annotation_header_split = split(annotation_header, ' ');
    string annotation_line;
    vector<string> annotation_line_split;
    vector<vector<string>> all_annotation_line_split;
    while(getline(annotation_file, annotation_line)){
        annotation_line_split = split(annotation_line, ' ');
        all_annotation_line_split.push_back(annotation_line_split);
    }
    unsigned long num_annotations = annotation_header_split.size();
    unsigned long num_snps = all_annotation_line_split.size();
    MatrixXf annotation_matrix(num_snps, num_annotations);
    for(unsigned int i = 0; i < num_snps; i++){
        for(unsigned int j= 0; j < num_annotations; j++){
            annotation_matrix(i,j)= stod(all_annotation_line_split[i][j]);
        }
    }
    
    //extract annotations specified in model
    MatrixXf output;
    if(model_annotations.size()>0){
        vector<int> annotation_index;
        bool annotation_found;
        for(unsigned int i = 0; i < model_annotations.size();  i++){
            annotation_found=false;
            for(unsigned int j=0; j < annotation_header_split.size(); j++){
                if(annotation_header_split[j].compare(model_annotations[i])==0){
                    annotation_index.push_back(j);
                    annotation_found=true;
                }
            }
            if(!annotation_found){
                cout << "Error: Annotation " << model_annotations[i] << " not found in annotation file" << endl;
            }
        }
        output.resize(num_snps, model_annotations.size()+1);
        VectorXf A0(num_snps);
        A0.setOnes();
        output.col(0) = A0;
        for(unsigned int i = 0; i < model_annotations.size(); i++){
            output.col(i+1) = annotation_matrix.col(annotation_index[i]);
        }
    }
    else{
        output.resize(num_snps, 1);
        VectorXf A0(num_snps);
        A0.setOnes();
        output.col(0) = A0;
    }
    out_annotations=output;
}

//read in all the input for PAINTOR based list of locus files specified in file_list_name
void Get_all_input(string& file_list_name, vector<Locus>& all_loci,  execution_parameters& params){
    
    ifstream input_list;
    input_list.open(file_list_name);
    string locus_name;
    string directory = params.in_dir;
    
    while(getline(input_list,locus_name)){
        Locus current_locus;
        current_locus.set_locus_name(locus_name);
        //Read in locus info and parse out z-scores
        cout << "Reading in files for: " << locus_name << endl;
        vector<VectorXf> locus_statistics;
        vector<string> snp_info;
        vector<MatrixXf > locus_ld;
        
        Read_Locus(current_locus, locus_name, params);
        current_locus.set_total_pops(params.zscore_header.size());
        
        
        //Read in LD for locus
        vector<string> LD_suffix = params.ld_suffix;
        vector<MatrixXf > locus_ld_matrices;
        for(unsigned int i = 0; i < LD_suffix.size(); i++){
            MatrixXf ld_matrix;
            string ld_name = locus_name+"."+LD_suffix[i];
            VectorXf transformed_z;
            Read_LD(directory, ld_name, ld_matrix);
            current_locus.set_LD(ld_matrix);
        }
        
        //Read in annotations for the locus
        string annot_name_filename = locus_name+ "." + params.annotation_suffix;
        vector<string> model_annotations = params.model_annotations;
        MatrixXf locus_annotations;
        Read_Annotations(directory, annot_name_filename, model_annotations, locus_annotations);
        current_locus.set_annotations(locus_annotations);
        all_loci.push_back(current_locus);
    }
    cout << "**********" <<  endl;
}

void Parse_command_line(int argc, const char * argv[], execution_parameters& parameters ){
    
    for(int i = 1; i < argc; i++){
        string argComp = argv[i];
        
        
        if(argComp.compare("-input")== 0){
            parameters.input_files = argv[i+1];
        }
        
        else if(argComp.compare("-in") == 0){
            string in_dir = argv[i+1];
            string last_char = string(&in_dir[in_dir.size()-1]);
            if(last_char.compare("/")!=0){
                in_dir = in_dir + "/";
            }
            parameters.in_dir = in_dir;
        }
        
        else if(argComp.compare("-out") == 0){
            string out_dir = argv[i+1];
            string last_char = string(&out_dir[out_dir.size()-1]);
            if(last_char.compare("/")!=0){
                out_dir = out_dir + "/";
            }
            parameters.out_dir = out_dir;
        }
        
        else if(argComp.compare("-LDname") == 0){
            string temp_ldnames = argv[i+1];
            size_t n = count(temp_ldnames.begin(), temp_ldnames.end(), ',');
            if(n >0) {
                parameters.ld_suffix = split(temp_ldnames, ',');
            }
            else{
                parameters.ld_suffix = {temp_ldnames};
            }
        }
        
        else if(argComp.compare("-Zhead") == 0){
            string header_temp = argv[i+1];
            size_t n = count(header_temp.begin(), header_temp.end(), ',');
            if(n>0){
                parameters.zscore_header = split(header_temp, ',');
            }
            else{
                parameters.zscore_header = {header_temp};
            }
        }
        
        else if(argComp.compare("-annotations") == 0){
            string input_annotations = argv[i+1];
            size_t n = count(input_annotations.begin(), input_annotations.end(), ',');
            if(n >0){
                parameters.model_annotations = split(input_annotations, ',');
            }
            else{
                parameters.model_annotations= {input_annotations};
            }
            
        }
        
        else if(argComp.compare("-Gname") == 0){
            parameters.gammaName = argv[i+1];
        }
        
        else if(argComp.compare("-Lname") == 0){
            parameters.likeli_name = argv[i+1];
        }
        
        else if(argComp.compare("-post1CV") == 0){
            parameters.single_post_flag = argv[i+1];
        }
        
        else if(argComp.compare("-MI") == 0){
            parameters.maxIter = stoi(argv[i+1]);
        }
        
        else if(argComp.compare("-ANname") == 0){
            parameters.annotation_suffix = argv[i+1];
        }
        
        else if(argComp.compare("-RESname") == 0){
            parameters.results_suffix = argv[i+1];
        }
        
        else if(argComp.compare("-gamma_initial") == 0){
            parameters.gamma_initial = argv[i+1];
            parameters.initialized_gammas = 1;
        }
        else if(argComp.compare("-enumerate")==0){
            parameters.max_causal = stoi(argv[i+1]);
            parameters.enumerate_flag = 1;
        }
        else if(argComp.compare("-set_seed")==0){
            parameters.sampling_seed = stoi(argv[i+1]);
        }
        else if(argComp.compare("-max_causal")==0){
            parameters.max_causal = stoi(argv[i+1]);
        }
        else if(argComp.compare("-prop_ld_eigenvalues")==0){
            parameters.prop_ld_eigenvalues = stof(argv[i+1]);
        }
        else if(argComp.compare("-mcmc")==0){
            parameters.mcmc_flag = 1;
        }
        else if(argComp.compare("-burn_in")==0){
            parameters.burn_in = stof(argv[i+1]);
        }
        else if(argComp.compare("-max_samples")==0) {
            parameters.max_samples = stof(argv[i + 1]);
        }
        else if(argComp.compare("-num_chains")==0) {
            parameters.num_chains = stoi(argv[i + 1]);
        }
    }
}

void display_execution_message(execution_parameters & parameters){
    if(parameters.mcmc_flag == 0 && parameters.mcmc_flag == 0){
        cout << "Error! Please specify either the -mcmc or -enumerate flags" << endl;
        exit(0);
    }
    if(parameters.mcmc_flag == 1){
        cout << "Running PAINTOR in MCMC mode" << endl;
        cout << "Number of chains: " << parameters.num_chains << endl;
        cout << "Samples per chain: " << parameters.max_samples << endl;
        cout << "Burn in: " << parameters.burn_in << endl;

    }
    else{
        cout << "Running PAINTOR with full enumeration" << endl;
        cout << "Maximum number of causals per locus: " << parameters.max_causal << endl;

    }
    cout << "Proportion of LD variance kept when performing truncated SVD for estimating N*h2g: " << parameters.prop_ld_eigenvalues << endl;
    cout << "Model annotations: ";
    for (auto i: parameters.model_annotations)
        cout << i << ' ';
    cout << endl;
    cout << "**********" <<  endl;
}

//initiailzie gamma such that the prior probability at each locus has 1 causal variant.
void initialize_gammas_input(vector<Locus>& all_loci, execution_parameters & parameters){
    vector<string> gamma_initial_split = split(parameters.gamma_initial, ',');
    VectorXf gamma_initial(gamma_initial_split.size());
    if(gamma_initial_split.size() != parameters.model_annotations.size()+1){
        cout << "Warning: Incorrect number of Enrichment parameters specified. Ignoring --gamma_initial flag" << endl;
    }
    else{
        for(unsigned int i =0; i < gamma_initial_split.size(); i++){
            gamma_initial[i] = stof(gamma_initial_split[i]);
        }
    }
    size_t num_loci = all_loci.size();
    for(int i = 0; i < num_loci; i++){
        all_loci[i].set_gamma(gamma_initial);
    }
}

void Welcome_Message(){
    cout << "Welcome to PAINTOR v3.0. Copyright Gleb Kichaev 2015. This software is free for academic use. Please e-mail gkichaev@ucla.edu for questions or reporting bugs" << endl;
    cout << "For required files and format specifications see User Manual \n \n" << endl;
    cout << "Usage: PAINTOR -input [input_filename] -in [input_directory] -out [output_directory] -Zhead [Zscore_header(s)] -LDname [LD_suffix(es)]  -annotations [annot_name1,annot_name2,...]  <other options> \n"<< endl;
    cout << "OPTIONS: -flag \t Description [default setting]  \n" << endl;
    cout << "-input \t (required) Filename of the input file containing the list of the fine-mapping loci [default: input.files]" << endl;
    cout << "-Zhead \t (required) The name(s) of the Zscores in the header of the locus file (comma separated) [default: N/A]" << endl;
    cout << "-LDname \t (required) Suffix(es) for LD files. Must match the order of Z-scores in locus file (comma separated) [Default:N/A]" << endl;
    cout << "-annotations \t The names of the annotations to include in model (comma separted) [default: N/A]" << endl;
    cout << "-in \t Input directory with all run files [default: ./ ]" << endl;
    cout << "-out \t Output directory where output will be written [default: ./ ]" << endl;
    cout << "-Gname \t Output Filename for enrichment estimates [default: Enrichment.Estimate]" << endl;
    cout << "-Lname \t Output Filename for log likelihood [default: Log.BayesFactor]" << endl;
    cout << "-RESname \t Suffix for ouput files of results [Default: results] " << endl;
    cout << "-ANname \t Suffix for annotation files [Default: annotations]" << endl;
    cout << "-MI \t Maximum iterations for algorithm to run [Default: 10]" << endl;
    cout << "-gamma_initial \t inititalize the enrichment parameters to a pre-specified value (comma separated) [Default: 0,...,0]" << endl;
    cout << "-num_samples  \t specify number of samples to draw for each locus [Default: 50000]" << endl ;
    cout << "-enumerate\t specify this flag if you want to enumerate all possible configurations followed by the max number of causal SNPs (eg. -enumerate 3 considers up to 3 causals at each locus) [Default: not specified]" << endl;
    cout << "-set_seed\t specify an integer as a seed for random number generator [default: clock time at execution]" << endl;
    cout << "-prop_ld_eigenvalues\t specify the proprotion of eigenvalues of LD matrix to keep when estimating the prior variance for each locus [default: 0.95]" << endl;
    cout << "-mcmc\t should the algorithm be run with MCMC? [Default: not specified]" << endl;
    cout << "-burn_in\t specify how many samples to discard during burn-in period [default: 50000]" << endl;
    cout << "-max_samples\t specify the number of samples to keep [default: 10000]" << endl;
    cout << "-num_chains\t specify the number of chains to run [default: 5]" << endl;
    cout << endl << endl ;
}

//append posterior probability vector to locus file and write output
void Locus::write_results(string& fname){
    ofstream file;
    file.open(fname);
    file << snp_info[0] + " Posterior_Prob\n";
    string out_line;
    for(int i = 0; i < num_snps ; i++){
        out_line = snp_info[i+1];
        out_line.erase(remove(out_line.begin(),out_line.end(), '\r'), out_line.end());
        file << out_line + " " + to_string(posterior_probability(i)) + "\n";
    }
    file.close();
}

void write_all_posteriors(vector<Locus>& all_loci, execution_parameters & parameters){
    ifstream input_list;
    input_list.open(parameters.input_files);
    string locus_name;
    size_t j=0;
    string out_name;
    while(getline(input_list,locus_name)){
        out_name = parameters.out_dir + locus_name + "."+ parameters.results_suffix;
        all_loci[j].write_results(out_name);
        j++;
    }
}

void write_other_output(VectorXf& gamma_estimates, double log_bf, execution_parameters & parameters){

    //write gamma estimates
    ofstream gamma_out;
    gamma_out.open( parameters.out_dir+parameters.gammaName);
    gamma_out << "Baseline";
    for(unsigned int i =0; i < parameters.model_annotations.size(); i++){
        gamma_out << " ";
        gamma_out << parameters.model_annotations[i];
    }
    gamma_out << endl;

    gamma_out << gamma_estimates[0];
    for(unsigned int i=1; i < gamma_estimates.size(); i ++){
        gamma_out << " ";
        gamma_out << gamma_estimates[i];
    }
    gamma_out << endl;
    gamma_out.close();

    //output log bayes factors
    ofstream likeli_out;
    likeli_out.open( parameters.out_dir+parameters.likeli_name);
    likeli_out<< std::setprecision(10) << log_bf;
    likeli_out << "\n";
    likeli_out.close();
}

void write_estimated_prior_variance(vector<Locus>& all_loci, execution_parameters & parameters){
    ifstream input_list;
    input_list.open(parameters.input_files);
    string locus_name;
    size_t j=0;
    string out_name;
    size_t num_pops = all_loci[0].num_pops;
    string output_log_fn = parameters.out_dir + "LogFile."+ parameters.results_suffix;
    ofstream file;
    int ld_pcs;
    file.open(output_log_fn);
    file << "locus_name prior_variance ld_pcs\n";
    double prior_var;
    while(getline(input_list,locus_name)){
        file << locus_name;
        file << " ";
        for(int i = 0; i < num_pops; i++){
            prior_var = all_loci[j].get_prior_variance(i);
            file << prior_var;
            if(i+1 <num_pops){
                file << ",";
            }
        }
        file << " ";
        for(int i = 0; i < num_pops; i++){
            ld_pcs = all_loci[j].ld_pcs[i];
            file << ld_pcs;
            if(i+1 <num_pops){
                file << ",";
            }
        }
        file << "\n";
        j++;
    }
}

///


