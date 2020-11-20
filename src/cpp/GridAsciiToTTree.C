// -----------------------------------------------------------------------------
// Brute forcing this a bit, gonna replace these using bash scripts
static const std::string process = 
  "PROCESS_MARKER"
;

static const std::string output_dir = 
  "OUTPUT_DIR_MARKER"
;
// -----------------------------------------------------------------------------


#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitter.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TMath.h"
#include "TMatrix.h"
#include "TMatrixT.h"
#include "TMinuit.h"
#include "TPave.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVectorT.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string> 
#include <time.h> 
#include <vector>

const int color_beamer_display[] = {kBlack,kBlue,kRed,kGreen+2,kViolet+1,kOrange+2,kMagenta,kAzure+1,kYellow+1,kCyan+1 };
double legendposition_header[4] = {0.17,0.93,0.95,0.995};


string prints(const char* format, ...){
	
	char buffer[255];
	va_list args;
	va_start (args, format);
	vsprintf (buffer,format, args);
	string output(buffer);
	va_end (args);
	
	return output;
	
	
}

vector<string> split(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}
	

void GridAsciiToTTree(){
	
	const unsigned int numberchiralities = 4;
	const string TGCdeviation = "0.0001";
	const string EnergyCMS = "250GeV";
	ifstream ifs;
	
	vector < string > file_description;
  map < string, unsigned int > n_process_variables;
  map < string, vector < unsigned int > > MAXFileIndicesOfProcesses;
  
  if (process == "ww_sl0muq") {
    file_description.push_back("grid_ww_sl0muq_leptonic");
    file_description.push_back("grid_ww_sl0muq_hadronic");
    
    n_process_variables["grid_ww_sl0muq_leptonic"] = 5;
    n_process_variables["grid_ww_sl0muq_hadronic"] = 5;
    
    MAXFileIndicesOfProcesses["grid_ww_sl0muq_leptonic"].push_back(20);
    MAXFileIndicesOfProcesses["grid_ww_sl0muq_leptonic"].push_back(10);
    MAXFileIndicesOfProcesses["grid_ww_sl0muq_leptonic"].push_back(10);
    MAXFileIndicesOfProcesses["grid_ww_sl0muq_leptonic"].push_back(1);
    MAXFileIndicesOfProcesses["grid_ww_sl0muq_hadronic"].push_back(20);
    MAXFileIndicesOfProcesses["grid_ww_sl0muq_hadronic"].push_back(10);
    MAXFileIndicesOfProcesses["grid_ww_sl0muq_hadronic"].push_back(10);
    MAXFileIndicesOfProcesses["grid_ww_sl0muq_hadronic"].push_back(1);
  } else if (process == "sw_sl0qq") {
    file_description.push_back("grid_sw_sl0qq_minus");
    file_description.push_back("grid_sw_sl0qq_plus");
    
    n_process_variables["grid_sw_sl0qq_minus"] = 6;
    n_process_variables["grid_sw_sl0qq_plus"] = 6;
    
    MAXFileIndicesOfProcesses["grid_sw_sl0qq_minus"].push_back(20);
    MAXFileIndicesOfProcesses["grid_sw_sl0qq_minus"].push_back(10);
    MAXFileIndicesOfProcesses["grid_sw_sl0qq_minus"].push_back(20);
    MAXFileIndicesOfProcesses["grid_sw_sl0qq_plus"].push_back(20);
    MAXFileIndicesOfProcesses["grid_sw_sl0qq_plus"].push_back(10);
    MAXFileIndicesOfProcesses["grid_sw_sl0qq_plus"].push_back(20);
  } else {
    throw std::invalid_argument("Unknown process " + process);
  }
  
	//file_description.push_back("grid_ww_sl0eq_leptonic");
	//file_description.push_back("grid_ww_sl0eq_hadronic");
	//file_description.push_back("grid_zz_sl0mu_down");
	
	//n_process_variables["grid_ww_sl0eq_leptonic"] = 5;
	//n_process_variables["grid_ww_sl0eq_hadronic"] = 5;
	//n_process_variables["grid_zz_sl0mu_down"] = 2;
	
	//MAXFileIndicesOfProcesses["grid_ww_sl0eq_leptonic"];
	//MAXFileIndicesOfProcesses["grid_ww_sl0eq_hadronic"];
	//MAXFileIndicesOfProcesses["grid_zz_sl0mu_down"].push_back(20);
	//MAXFileIndicesOfProcesses["grid_zz_sl0mu_down"].push_back(10);
	//MAXFileIndicesOfProcesses["grid_zz_sl0mu_down"].push_back(10);
	//MAXFileIndicesOfProcesses["grid_zz_sl0mu_down"].push_back(1);
	
	const string input_path = output_dir + "/omega/";
	const string output_path = output_dir + "/grids_root/";
	
	
	for(unsigned int n = 0; n < file_description.size(); n++){
		
		const unsigned int numbervariables = n_process_variables[file_description[n]];
		
		cout << file_description[n] << "\t" << numbervariables << endl;
		time_t start_file_progressing_time = time(0);
		
		unsigned int MaxNumberRelements = 0;
		vector < vector < double > > gird_variables;
		vector < vector < double > > grid_SM_elements;
		vector < vector < vector < double > > > grid_R_elements;
		
		unsigned int max_file_number[4];
		
		double total_file_number = 1;
		for(unsigned int i = 0; i < 4; i++){
			max_file_number[i] = MAXFileIndicesOfProcesses[file_description[n]][i];
			total_file_number *= (double)max_file_number[i];
		}
		
		double total_file_count = 0;
		
		for(unsigned int i = 0; i < max_file_number[0]; i++){
		for(unsigned int j = 0; j < max_file_number[1]; j++){
		for(unsigned int k = 0; k < max_file_number[2]; k++){
		for(unsigned int l = 0; l < max_file_number[3]; l++){
			
			unsigned int file_index[4] = {i+1,j+1,k+1,l+1};
			string currentfile = prints("%s%s_%s",input_path.c_str(),file_description[n].c_str(),EnergyCMS.c_str());
			for(unsigned int m = 0; m < 4; m++){
				if(max_file_number[m]>1){
					currentfile += prints("_%i",file_index[m]);
				}
			}
			currentfile += ".txt";
			
			if( (int)(total_file_count)%(int)(0.05*total_file_number) == 0 ){
				cout << 100 * total_file_count / total_file_number << "% of the files done after " << time(0) - start_file_progressing_time << "s" << endl;
			}
			
			
			ifs.open( currentfile.c_str(), std::ifstream::in);
			if(ifs){
				bool print_out_file_name = false;
				while (ifs.good()) {
					string current_line;
					getline(ifs,current_line);
					vector<string> stringentries = split(current_line," ");
					vector < double > currentvalues;
					for(unsigned int m = 0; m < stringentries.size(); m++){
						currentvalues.push_back( stod( stringentries[m] ) );
					}
					if(currentvalues.size() > 1){
						vector < double >::iterator current_R_counter = currentvalues.begin();
						gird_variables.push_back(vector<double>(current_R_counter,current_R_counter+numbervariables));
						current_R_counter += numbervariables;
						grid_SM_elements.push_back(vector<double>(current_R_counter,current_R_counter+numberchiralities));
						current_R_counter += numberchiralities;
						vector < vector < double > > current_grid_R_elements;
						while( current_R_counter+numberchiralities <= currentvalues.end() ){
							current_grid_R_elements.push_back(vector<double>(current_R_counter,current_R_counter+numberchiralities));
							current_R_counter+=numberchiralities;
						}
						
						if(current_grid_R_elements.size()>0)
							grid_R_elements.push_back(current_grid_R_elements);
						
						if(current_grid_R_elements.size()>MaxNumberRelements)
							MaxNumberRelements = current_grid_R_elements.size();
						
						for(unsigned n = 0; n < grid_SM_elements[grid_SM_elements.size()-1].size(); n++){
							if(grid_SM_elements[grid_SM_elements.size()-1][n] > 10){
								print_out_file_name = true;
							}
						}
							
					}
				}
				
				if(print_out_file_name){
					cout << currentfile << endl;
				}
				
			}else{
				cout << "file " << currentfile << " couldn't be opend" << endl;
			}
			ifs.close();
			total_file_count++;
			
		}
		}
		}	
		}
		
		cout << gird_variables.size()<< "\t" << grid_SM_elements.size() << "\t" << grid_R_elements.size() << endl;
		
		if(gird_variables.size() == grid_SM_elements.size() && gird_variables.size() == grid_R_elements.size()){
		
			double variables[numbervariables];
			double SM_elements[numberchiralities];
			double R_elements[MaxNumberRelements][numberchiralities];
			
			string variables_description = "par";
			string SM_elements_description = "SM";
			string R_elements_description = "R";
			
			TFile* outputfile = new TFile( ( output_path + file_description[n] + "_" + EnergyCMS + "_" + TGCdeviation + ".root").c_str() , "RECREATE" );
			TTree* matrix_elements_distribution = new TTree( "MatrixElements", "" );
			
			matrix_elements_distribution->Branch(variables_description.c_str(),		variables,		prints("%s[%i]/D",variables_description.c_str(),numbervariables).c_str());
			matrix_elements_distribution->Branch(SM_elements_description.c_str(),	SM_elements,	prints("%s[%i]/D",SM_elements_description.c_str(),numberchiralities).c_str());
			matrix_elements_distribution->Branch(R_elements_description.c_str(),	R_elements,		prints("%s[%i][%i]/D",R_elements_description.c_str(),MaxNumberRelements,numberchiralities).c_str());
		
			for(unsigned int i = 0; i < gird_variables.size(); i++){
				
				for(unsigned int j = 0; j < numbervariables; j++){
					variables[j] = gird_variables[i][j];
				}
				for(unsigned int j = 0; j < numberchiralities; j++){
					SM_elements[j] = grid_SM_elements[i][j];
				}
				for(unsigned int j = 0; j < grid_R_elements[i].size(); j++){
					for(unsigned int k = 0; k < numberchiralities; k++){
						R_elements[j][k] = grid_R_elements[i][j][k];
					}
				}
				
				matrix_elements_distribution->Fill();
					
			}
			matrix_elements_distribution->Write();
			outputfile->Write();
			outputfile->Close();
			
		}else if(gird_variables.size() == grid_SM_elements.size()){
			cerr << "Array size for R not matching -> ignoring R" << endl;
			
			double variables[numbervariables];
			double SM_elements[numberchiralities];
			
			string variables_description = "par";
			string SM_elements_description = "SM";
			
			TFile* outputfile = new TFile( ( output_path + file_description[n] + "_" + EnergyCMS + "_" + TGCdeviation + ".root").c_str() , "RECREATE" );
			TTree* angular_distribution = new TTree( "MatrixElements", "" );
			
			angular_distribution->Branch(variables_description.c_str(),		variables,		prints("%s[%i]/D",variables_description.c_str(),numbervariables).c_str());
			angular_distribution->Branch(SM_elements_description.c_str(),	SM_elements,	prints("%s[%i]/D",SM_elements_description.c_str(),numberchiralities).c_str());
		
			for(unsigned int i = 0; i < gird_variables.size(); i++){
				
				for(unsigned int j = 0; j < numbervariables; j++){
					variables[j] = gird_variables[i][j];
				}
				for(unsigned int j = 0; j < numberchiralities; j++){
					SM_elements[j] = grid_SM_elements[i][j];
				}
				
				angular_distribution->Fill();
					
			}
			angular_distribution->Write();
			outputfile->Write();
			outputfile->Close();
			
		}else{
			cerr << "Array size not matching" << endl;
		}
		
		cout << "File finished after " << time(0) - start_file_progressing_time << "s" << endl;
		
	}
	
	
}
