//
// Created by Gleb Kichaev on 10/21/15.
//

#include "Functions_IO.h"


//Functions to read input

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

void Read_Locus(string &input_directory, string& fname, vector<string> & zname,  vector<VectorXd>& all_association_stat, vector<string>& snp_info, string & header){

    ifstream locus;
    locus.open(input_directory+fname);
    getline(locus, header);
    char delimiter = ' ';
    vector<string> split_header = split(header, delimiter);

    //vector<string> fields = {"SNP_ID", "CHR", "POS", "EFFECT_ALLELE"};
    vector<int> z_index;
    for(unsigned int i=0; i < split_header.size(); i++){
        for(unsigned int j=0; j < zname.size(); j++){
            if(zname[j].compare(split_header[i])==0){
                z_index.push_back(i);
            }
        }
    }
    string input_line;
    vector<string> split_line;
    vector<vector<string>> snp_info_split;
    string stat_holder;

    while(getline(locus,input_line)){
        split_line = split(input_line, ' ');
        snp_info.push_back(input_line);
        snp_info_split.push_back(split_line);
    }
    locus.close();
    unsigned long num_snps = snp_info.size();
    VectorXd association_stat;
    for(unsigned int i =0; i < z_index.size();i++){
        association_stat.resize(num_snps,1);
        for(unsigned int j=0; j < num_snps; j++){
            if(snp_info_split[j][z_index[i]].compare("NA") == 0){
                association_stat[j] = 0;
            }
            else{
                association_stat[j] = stod(snp_info_split[j][z_index[i]]);
            }
        }
        all_association_stat.push_back(association_stat);
    }

}

double Regularize_LD(MatrixXd & ld_mat){
    double reg_factor = 0.001;
    int reg_flag = 0;
    double det;
    VectorXd X(ld_mat.rows());
    MatrixXd Y = X.asDiagonal();
    while(reg_flag == 1){
        X.fill(reg_factor);
        Y = X.asDiagonal();
        det = (ld_mat+Y).determinant();
        if(det <= 1e-100){
            reg_factor = reg_factor*1.25;
        }
        else{
            reg_flag = 1;
        }

    }
    X.fill(reg_factor);
    Y = X.asDiagonal();
    ld_mat= ld_mat +Y;
    return reg_factor;
}

void Read_LD(string &input_directory, string& fname , MatrixXd& ld_matrix){
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
    MatrixXd numeric_ld(num_snps, num_snps);
    double val;
    for(unsigned int i = 0; i < num_snps; i++){
        for(unsigned int j = 0; j < num_snps; j++){
            val =  stod(all_ld_line_split[i][j]);
            if(isfinite(val)) {
                numeric_ld(i, j) = stod(all_ld_line_split[i][j]);
            }
            else{
                cout <<  "Error: Infinite value in LD matrix detected!!" << endl;
            }
        }
    }
    double reg_factor;
    reg_factor = Regularize_LD(numeric_ld);
    ld_matrix = numeric_ld;
}

void Get_All_Input(string& file_list_name, string& directory, vector<string> znames, vector<string>& model_annotations, vector<vector<VectorXd>>& all_locus_statistics, vector<vector<MatrixXd>> &all_ld_matrices,  vector<MatrixXd>& all_annotations, vector<string>& LD_suffix, string& annotation_suffix,
                   vector<vector<string>>& all_SNP_info, vector<string>& all_headers) {

    ifstream input_list;
    input_list.open(file_list_name);
    string locus_name;
    string header;

    while(getline(input_list,locus_name)){
        //Read in locus info and parse out z-scores
        cout << "Reading in files for: " << locus_name << endl;
        vector<VectorXd> locus_statistics;
        vector<string> snp_info;
        vector<MatrixXd> locus_ld;
        Read_Locus(directory, locus_name, znames, locus_statistics, snp_info, header);
        all_locus_statistics.push_back(locus_statistics);
        all_headers.push_back(header);
        all_SNP_info.push_back(snp_info);

        //Read in LD for locus
        vector<MatrixXd> locus_ld_matrices;
        for(unsigned int i = 0; i < LD_suffix.size(); i++){
            MatrixXd ld_matrix;
            string ld_name = locus_name+"."+LD_suffix[i];
            VectorXd transformed_z;
            Read_LD(directory, ld_name, ld_matrix);
            locus_ld_matrices.push_back(ld_matrix);
        }
        all_ld_matrices.push_back(locus_ld_matrices);

        //Read in annotations for the locus
        string annot_name = locus_name+"."+annotation_suffix;
        vector<string> annotation_header;
        MatrixXd locus_annotations;
        Read_Annotations(directory, annot_name, model_annotations, locus_annotations, annotation_header);
        all_annotations.push_back(locus_annotations);
    }

}


void Read_Annotations(string &input_directory, string& fname , vector<string>& model_annotations, MatrixXd& out_annotations, vector<string> & annotation_header_split){
    ifstream annotation_file;
    annotation_file.open(input_directory+fname);
    string annotation_header;
    getline(annotation_file,annotation_header);
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
    MatrixXd annotation_matrix(num_snps, num_annotations);
    for(unsigned int i = 0; i < num_snps; i++){
        for(unsigned int j= 0; j < num_annotations; j++){
            annotation_matrix(i,j)= stod(all_annotation_line_split[i][j]);
        }
    }
    MatrixXd output;
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
        VectorXd A0(num_snps);
        A0.setOnes();
        output.col(0) = A0;
        for(unsigned int i = 0; i < model_annotations.size(); i++){
            output.col(i+1) = annotation_matrix.col(annotation_index[i]);
        }
    }
    else{
        output.resize(num_snps, 1);
        VectorXd A0(num_snps);
        A0.setOnes();
        output.col(0) = A0;
    }
    out_annotations=output;
}

//Functions to output results


void Write_Posterior(string& out_dir, string & out_name, VectorXd& locus_results, vector<string>& locus_info, string& header ){

    vector<double> probs_as_strings;
    ofstream myfile;
    string fname = out_dir+out_name;
    myfile.open(fname);
    for(int i =0; i<locus_results.size(); i++){
        probs_as_strings.push_back((locus_results[i]));
    }
    myfile << header + " Posterior_Prob" << endl;
    for(int i=0 ;i<locus_results.size(); i++){
        myfile << locus_info[i] + " ";
        myfile << probs_as_strings[i];
        myfile << "\n";
    }

    myfile.close();


}

void Write_All_Output(string& input_files, string& out_dir, string& out_suffix, CausalProbs& results, vector<vector<string>> & all_locus_info,
                      VectorXd & gamma_estimates, string& Gname, double log_likeli, string & Lname, vector<string>& all_headers, vector<string>& annot_names){
    ifstream input_list;
    input_list.open(input_files);
    string locus_name;
    unsigned int j=0;
    while(getline(input_list,locus_name)){
        string out_fname = locus_name + "." + out_suffix;
        Write_Posterior(out_dir, out_fname, results.probs_locs[j], all_locus_info[j], all_headers[j]);
        j++;
        if(j >= all_locus_info.size()) {
            break;
        }
    }
    input_list.close();

    ofstream gamma_out;
    string G_open = out_dir+Gname;
    gamma_out.open(G_open);

    gamma_out << "Baseline";
    for(unsigned int i =0; i < annot_names.size(); i++){
        gamma_out << " ";
        gamma_out << annot_names[i];
    }
    gamma_out << endl;

    gamma_out << gamma_estimates[0];
    for(unsigned int i=1; i < gamma_estimates.size(); i ++){
        gamma_out << " ";
        gamma_out << gamma_estimates[i];
    }
    gamma_out << endl;
    gamma_out.close();

    string L_open = out_dir+Lname;
    ofstream likeli_out;
    likeli_out.open(L_open);
    likeli_out<< std::setprecision(10) << log_likeli;
    likeli_out << "\n";
    likeli_out.close();

}

void Reshape_Input(vector<vector<VectorXd>>& all_locus_statistics, vector<vector<MatrixXd>> &all_ld_matrices, vector<VectorXd>& reshape_stats, vector<MatrixXd> &reshape_ld ){
    for(unsigned int i=0; i< all_locus_statistics.size(); i++){
        reshape_stats.push_back(all_locus_statistics[i][0]);
        reshape_ld.push_back(all_ld_matrices[i][0]);
    }
}



///


