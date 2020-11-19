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
const string path_angular_distributions = output_dir + "/distributions/angular/";
const string path_PNPC_distributions = output_dir + "/distributions/TGC/";
const string path_plots = output_dir + "/distributions/plots/";

unsigned int global_hist_count = 0;
ofstream logfile;

double GetMaximum( vector < double > values){
	
	double max = -TMath::Infinity();
	for(unsigned int i = 0; i < values.size(); i++){
		if( max < values[i] )
			max = values[i];
	}
	
	return max;
}

double GetMinimum( vector < double > values){
	
	double min = TMath::Infinity();
	for(unsigned int i = 0; i < values.size(); i++){
		if( min > values[i] )
			min = values[i];
	}
	
	return min;
}

double GetN0Minimum( vector < double > values){
	
	double min = TMath::Infinity();
	for(unsigned int i = 0; i < values.size(); i++){
		if( values[i] != 0 && min > values[i] )
			min = values[i];
	}
	
	return min;
}


double GetL0Minimum( vector < double > values){
	
	double min = TMath::Infinity();
	for(unsigned int i = 0; i < values.size(); i++){
		if( values[i] > 0 && min > values[i] )
			min = values[i];
	}
	
	return min;
}



void TH2D_RightMarginCorrection(){
	gPad->SetRightMargin(0.15);
}

string prints(const char* format, ...){
	
	char buffer[255];
	va_list args;
	va_start (args, format);
	vsprintf (buffer,format, args);
	string output(buffer);
	va_end (args);
	
	return output;
	
	
}

void read_grid_files(
vector < vector < double > > &grid_variables,
vector < vector < double > > &grid_SM_elements,
vector < vector < vector < double > > > &grid_R_elements,
const unsigned int numbervariables,
const unsigned int numberchiralities,
const unsigned int MaxNumberRelements,
vector < string > file_description,
vector < vector < double > > chiral_scaling,
string EnergyCMS,
string TGCdeviation,
bool filebased_nomalization = false
){
	

	double variables[numbervariables];
	double SM_elements[numberchiralities];
	double R_elements[MaxNumberRelements][numberchiralities];
	

	string variables_description = "par";
	string SM_elements_description = "SM";
	string R_elements_description = "R";
	
	vector < vector < vector < double > > > grid_SM_elements_file_separated(file_description.size());
	vector < vector < vector < vector < double > > > > grid_R_elements_file_separated(file_description.size());
	
	for(unsigned int n = 0; n < file_description.size(); n++){

		TFile* inputfile = new TFile( ( file_description[n] + "_" + EnergyCMS + "_" + TGCdeviation + ".root").c_str() , "READ" );
		if(inputfile->IsOpen()){
			TTree* matrix_elements_distribution = (TTree*)inputfile->Get( "MatrixElements" );

			if(numbervariables>0)
				matrix_elements_distribution->SetBranchAddress(variables_description.c_str(),	variables	);
			if(numberchiralities>0)
				matrix_elements_distribution->SetBranchAddress(SM_elements_description.c_str(),	SM_elements	);
			if(MaxNumberRelements>0)
				matrix_elements_distribution->SetBranchAddress(R_elements_description.c_str(),	R_elements	);

			if(numbervariables>0 || numberchiralities>0 || MaxNumberRelements>0){
				for(unsigned int i = 0; i < (unsigned int)matrix_elements_distribution->GetEntries(); i++){
					
					matrix_elements_distribution->GetEntry(i);
					
					if(numbervariables>0){
						vector < double > current_grid_variables;
						for(unsigned int j = 0; j < numbervariables; j++){
							current_grid_variables.push_back(variables[j]);

						}
						grid_variables.push_back(current_grid_variables);
					}
					
					if(numberchiralities>0){
						vector < double > current_grid_SM_elements;
						for(unsigned int j = 0; j < numberchiralities; j++){
							current_grid_SM_elements.push_back(chiral_scaling[n][j]*SM_elements[j]);

						}
						//grid_SM_elements.push_back(current_grid_SM_elements);
						grid_SM_elements_file_separated[n].push_back(current_grid_SM_elements);
					}
					
					if(MaxNumberRelements>0){
						vector < vector < double > > current_grid_R_elements(MaxNumberRelements);					
						for(unsigned int j = 0; j < MaxNumberRelements; j++){
							for(unsigned int k = 0; k < numberchiralities; k++){
								current_grid_R_elements[j].push_back(chiral_scaling[n][k]*R_elements[j][k]);
							}
						}
						//grid_R_elements.push_back(current_grid_R_elements);
						grid_R_elements_file_separated[n].push_back(current_grid_R_elements);
					}
						
				}
			}
		}
		else{
			cout << "File " << (file_description[n] + "_" + EnergyCMS + "_" + TGCdeviation + ".root") << " could not be opend" << endl;
		}
		inputfile->Close();
	}

	if(filebased_nomalization){
		for(unsigned int n = 0; n < file_description.size(); n++){
			double sm_nborms[numberchiralities];
			for(unsigned int chiral = 0; chiral < numberchiralities; chiral++){
				sm_nborms[chiral]=0;
			}
			for(unsigned int i = 0; i < grid_SM_elements_file_separated[n].size(); i++ ){
				for(unsigned int j = 0; j < grid_SM_elements_file_separated[n][i].size(); j++){
					sm_nborms[j] += grid_SM_elements_file_separated[n][i][j];
				}
			}
			for(unsigned int i = 0; i < grid_SM_elements_file_separated[n].size(); i++ ){
				for(unsigned int j = 0; j < grid_SM_elements_file_separated[n][i].size(); j++){
					if(sm_nborms[j] > 0){
						grid_SM_elements_file_separated[n][i][j] /= sm_nborms[j];
					}
				}
			}
			for(unsigned int i = 0; i < grid_R_elements_file_separated[n].size(); i++){
				for(unsigned int j = 0; j < grid_R_elements_file_separated[n][i].size(); j++){
					for(unsigned int k = 0; k < grid_R_elements_file_separated[n][i][j].size(); k++){
						if(sm_nborms[k] > 0){
							grid_R_elements_file_separated[n][i][j][k] /= sm_nborms[k];
						}
					}
				}
			}
			
		}
	}
	
	for(unsigned int n = 0; n < file_description.size(); n++){
		for(unsigned int i = 0; i < grid_SM_elements_file_separated[n].size(); i++ ){
			grid_SM_elements.push_back(grid_SM_elements_file_separated[n][i]);
		}
		for(unsigned int i = 0; i < grid_R_elements_file_separated[n].size(); i++){
			grid_R_elements.push_back(grid_R_elements_file_separated[n][i]);
		}
	}

}

void write_angular_file( 
	string filename,
	string TTreeName,
	vector < vector < double > > Distribution,
	vector < vector < double > > center,
	vector < vector < double > > width,
	vector < string > chiral_description
){
	
	if( Distribution.size() == center.size() && Distribution.size() == width.size() ){
		map< string, unsigned int > chiral_assigment;
		
		chiral_assigment["eLpR"] = 0;
		chiral_assigment["eRpL"] = 1;
		chiral_assigment["eRpR"] = 2;
		chiral_assigment["eLpL"] = 3;
			
		double chiral_angular_value[4] = {0,0,0,0};
		int angle_dim = 0;
		double angle_center[256];
		double angle_width[256];
		
		for( unsigned int i = 0; i < Distribution.size(); i++){
			angle_dim = TMath::Max( angle_dim, (int)center[i].size() );
			angle_dim = TMath::Max( angle_dim, (int)width[i].size() );
		}
		
		
		TFile* outputfile = new TFile( filename.c_str() , "RECREATE" );
		TTree* angular_distribution = new TTree( TTreeName.c_str(), "" );
		
		angular_distribution->Branch("angle_center",		angle_center,				prints("angle_center[%i]/D",angle_dim).c_str() );
		angular_distribution->Branch("angle_width",			angle_width,				prints("angle_width[%i]/D",angle_dim).c_str() );
		angular_distribution->Branch("relative_sigma_LR",	&chiral_angular_value[0],	"relative_sigma_LR/D");
		angular_distribution->Branch("relative_sigma_RL",	&chiral_angular_value[1],	"relative_sigma_RL/D");
		angular_distribution->Branch("relative_sigma_RR",	&chiral_angular_value[2],	"relative_sigma_RR/D");
		angular_distribution->Branch("relative_sigma_LL",	&chiral_angular_value[3],	"relative_sigma_LL/D");
		
		for(unsigned int i = 0; i < Distribution.size(); i++){
			
			for( int j = 0; j < angle_dim; j++){
				angle_center[j] = 0;
				angle_width[j] = 0;
			}
			
			for( unsigned int j = 0; j < center[i].size(); j++){
				angle_center[j] = center[i][j];
			}
			for( unsigned int j = 0; j < width[i].size(); j++){
				angle_width[j] = width[i][j];
			}
			
			for(unsigned int j = 0; j < Distribution[i].size(); j++){
				chiral_angular_value[chiral_assigment[chiral_description[j]]] = Distribution[i][j];
			}
			
			angular_distribution->Fill();
			
		}
		angular_distribution->Write();
		outputfile->Write();
		outputfile->Close();
		
	}else{
		cerr << " array lenghts are not matching, " << filename << " has not been created!!! " << endl;
	}
}

void write_PNPC_file( 
	string filename,
	string TTreeName,
	vector < vector < vector < double > > > coefficients,
	vector < vector < double > > center,
	vector < vector < double > > width,
	vector < string > PNPC_description,
	vector < string > chiral_description
){
	
	if( coefficients.size() == center.size() && coefficients.size() == width.size() ){
		
		map< string, unsigned int > chiral_assigment;		
		chiral_assigment["eLpR"] = 0;
		chiral_assigment["eRpL"] = 1;
		chiral_assigment["eRpR"] = 2;
		chiral_assigment["eLpL"] = 3;
		
		int angle_dim = 0;
		int coefficiont_dim = 0;
		double chiral_angular_value[4][256];
		double angle_center[256];
		double angle_width[256];
		
		string coefficiont_label = "";
		
		for( unsigned int i = 0; i < coefficients.size(); i++){
			angle_dim = TMath::Max( angle_dim, (int)center[i].size() );
			angle_dim = TMath::Max( angle_dim, (int)width[i].size() );
			coefficiont_dim = TMath::Max( coefficiont_dim, (int)coefficients[i].size() );
		}
		
		
		TFile* outputfile = new TFile( filename.c_str() , "RECREATE" );
		TTree* angular_distribution = new TTree( TTreeName.c_str(), "" );
		
		angular_distribution->Branch("PNPC_Label",			&coefficiont_label);
		angular_distribution->Branch("angle_center",		angle_center,				prints("angle_center[%i]/D",angle_dim).c_str() );
		angular_distribution->Branch("angle_width",			angle_width,				prints("angle_width[%i]/D",angle_dim).c_str() );
		angular_distribution->Branch("relative_PNPC_LR",	chiral_angular_value[0],	prints("relative_PNPC_LR[%i]/D",coefficiont_dim).c_str() );
		angular_distribution->Branch("relative_PNPC_RL",	chiral_angular_value[1],	prints("relative_PNPC_RL[%i]/D",coefficiont_dim).c_str() );
		angular_distribution->Branch("relative_PNPC_RR",	chiral_angular_value[2],	prints("relative_PNPC_RR[%i]/D",coefficiont_dim).c_str() );
		angular_distribution->Branch("relative_PNPC_LL",	chiral_angular_value[3],	prints("relative_PNPC_LL[%i]/D",coefficiont_dim).c_str() );
		
		
		for(unsigned int i = 0; i < coefficients.size(); i++){
			
			for( int j = 0; j < angle_dim; j++){
				angle_center[j] = 0;
				angle_width[j] = 0;
			}
			
			for( unsigned int j = 0; j < center[i].size(); j++){
				angle_center[j] = center[i][j];
			}
			for( unsigned int j = 0; j < width[i].size(); j++){
				angle_width[j] = width[i][j];
			}
			
			for( int j = 0; j < coefficiont_dim; j++){
				for(unsigned int k = 0; k < 4; k++){
					chiral_angular_value[k][j] = 0;
				}
			}
			
			coefficiont_label = "";
			for(unsigned int j = 0; j < PNPC_description.size(); j++){
				coefficiont_label += PNPC_description[j] + ";";
			}
			
			for(unsigned int j = 0; j < coefficients[i].size(); j++){
				for(unsigned int k = 0; k < coefficients[i][j].size(); k++){
					chiral_angular_value[chiral_assigment[chiral_description[k]]][j] = coefficients[i][j][k];
				}
			}
			angular_distribution->Fill();
			
		}
		angular_distribution->Write();
		outputfile->Write();
		outputfile->Close();
		
	}else{
		cerr << " array lenghts are not matching, " << filename << " has not been created!!! " << endl;
	}
	
	
}

template< typename Hist >
vector < vector< Hist* > > generate_TGC_coefficients( vector < double > TGC_scale, vector < vector< Hist* > > Hist_Rs, vector < string > chiral_label, string process_label ){
	
	logfile << Hist_Rs.size() << " coefficients for " << process_label << endl;
	vector < vector< Hist* > > TGC_coefficient(Hist_Rs.size());

	if(TGC_scale.size() == chiral_label.size()){
		for(unsigned int i = 0; i < TGC_coefficient.size(); i++){
			TGC_coefficient[i].resize(chiral_label.size());
		}
		
		
		string coefficient_labeling[] = {"A","B","C","D","E","F","G","H","I"};
	
		
		for(unsigned int i = 0; i < chiral_label.size(); i++){
			
			logfile << "Start chirality " << chiral_label[i] << endl;
			
			TGC_coefficient[0][i] = (Hist*)Hist_Rs[0][i]->Clone( (process_label + "_" + chiral_label[i] + "_0").c_str() );
			logfile << "Integral " << TGC_coefficient[0][i]->Integral() << endl;
			
			for(unsigned int j = 0; j < 3; j++){
				
				logfile << "Proceed coefficient " << coefficient_labeling[j] << endl;						
				TGC_coefficient[j+1][i] = (Hist*)Hist_Rs[j+1][i]->Clone( (process_label + "_" + chiral_label[i] + "_" + coefficient_labeling[j] ).c_str() );			
				TGC_coefficient[j+1][i]->Add(Hist_Rs[j+4][i],-1);	
				TGC_coefficient[j+1][i]->Scale(0.5  / TGC_scale[i]);		

				logfile << "Proceed coefficient " << coefficient_labeling[j+3] << endl;
				TGC_coefficient[j+4][i] = (Hist*)Hist_Rs[j+1][i]->Clone( (process_label + "_" + chiral_label[i] + "_" + coefficient_labeling[j+3] ).c_str() );
				TGC_coefficient[j+4][i]->Add(Hist_Rs[j+4][i],1);
				TGC_coefficient[j+4][i]->Add(Hist_Rs[0][i],-2);
				TGC_coefficient[j+4][i]->Scale(0.5 * TMath::Power(TGC_scale[i],-2));
				
			}
			logfile << "Proceed coefficient G" << endl;
			TGC_coefficient[7][i] = (Hist*)Hist_Rs[7][i]->Clone( (process_label + "_" + chiral_label[i] + "_G").c_str() );
			TGC_coefficient[7][i]->Add(Hist_Rs[1][i],-1);
			TGC_coefficient[7][i]->Add(Hist_Rs[2][i],-1);
			TGC_coefficient[7][i]->Add(Hist_Rs[0][i],1);
			TGC_coefficient[7][i]->Scale( TMath::Power(TGC_scale[i],-2) );
			logfile << "Proceed coefficient H" << endl;
			TGC_coefficient[8][i] = (Hist*)Hist_Rs[9][i]->Clone( (process_label + "_" + chiral_label[i] + "_H").c_str() );
			TGC_coefficient[8][i]->Add(Hist_Rs[1][i],-1);
			TGC_coefficient[8][i]->Add(Hist_Rs[3][i],-1);
			TGC_coefficient[8][i]->Add(Hist_Rs[0][i],1);
			TGC_coefficient[8][i]->Scale( TMath::Power(TGC_scale[i],-2) );
			logfile << "Proceed coefficient I" << endl;
			TGC_coefficient[9][i] = (Hist*)Hist_Rs[8][i]->Clone( (process_label + "_" + chiral_label[i] + "_I").c_str() );
			TGC_coefficient[9][i]->Add(Hist_Rs[2][i],-1);
			TGC_coefficient[9][i]->Add(Hist_Rs[3][i],-1);
			TGC_coefficient[9][i]->Add(Hist_Rs[0][i],1);
			TGC_coefficient[9][i]->Scale( TMath::Power(TGC_scale[i],-2) );
			
			logfile << "End chirality " << chiral_label[i] << endl;
			
			
		}
		
		logfile << "begin Normalizing" << endl;
		/*
		for(unsigned int i = 0; i < chiral_label.size(); i++){
			for(unsigned int m = 0; m < 5; m++){
				cout <<TGC_coefficient[0][i]->GetBinContent(m+1) << "\t";
			}
		cout << endl;
		}
		*/
		for(unsigned int i = 1; i < TGC_coefficient.size(); i++){
			for( unsigned int j = 0; j < TGC_coefficient[i].size(); j++){
				/*
				for(unsigned int m = 0; m < 5; m++){
					cout << TGC_coefficient[i][j]->GetBinContent(m+1) << "\t";
				}
				cout << endl;
				*/
				if(TGC_coefficient[0][j]->Integral() != 0){
					TGC_coefficient[i][j]->Divide(TGC_coefficient[0][j]);
				}else{
					logfile << "Integration for " << chiral_label[j] << " yields 0 no division" << endl;
				}
			}
		}
	}else{
		logfile << "ERROR: TGC_scale and chirlaity vector do NOT have the same size\n No output will be created\n" << endl;
	}
	logfile << "return\n" << endl;
	return TGC_coefficient;
	
	
}

void GeneratedHistConfigMultiAngularDistribution(
vector < vector < double > > grid_variables,
vector<unsigned int> variable_index,
vector<unsigned int> & Genereated_binning,
vector<double> & Genereated_range_min,
vector<double> & Genereated_range_max
){
	vector < vector < double > > bin_centers(3);
	vector < map < double, bool > > bin_center_used(3);
	
	for(unsigned int i = 0; i < grid_variables.size(); i++){
		for(unsigned int j = 0; j < 3; j++){
			if(variable_index[j] < grid_variables[i].size()){
				if(!bin_center_used[j][grid_variables[i][variable_index[j]]]){
					bin_centers[j].push_back(grid_variables[i][variable_index[j]]);
					bin_center_used[j][grid_variables[i][variable_index[j]]] = true;
				}
			}else{
				cout << "index" << variable_index[j] << " exseeds vector range!" << endl;
			}
		}
	}
	
	for(unsigned int i = 0; i < 3; i++){
		Genereated_binning[i] = bin_centers[i].size();
		Genereated_range_min[i] = GetMinimum( bin_centers[i] );
		Genereated_range_max[i] = GetMaximum( bin_centers[i] );
		if(Genereated_binning[i] > 1){
			double bin_width = TMath::Abs( ( Genereated_range_max[i] - Genereated_range_min[i] ) / ((double)(Genereated_binning[i]-1)) );
			Genereated_range_min[i] -= 0.5*bin_width;
			Genereated_range_max[i] += 0.5*bin_width;
		}
	}
}

void CalculateMultiAngularDistribution(
vector < vector < double > > &Distribution,
vector < vector < double > > &center,
vector < vector < double > > &width,
vector < vector < double > > grid_variables,
vector < vector < double > > grid_SM_elements,
vector<unsigned int> variable_index,
vector<unsigned int> binning,
vector<double> range_min,
vector<double> range_max
){
	
	Distribution.clear();
	center.clear();
	width.clear();
	
	unsigned int max_chiralities = 1;
	for(unsigned int i = 0; i < grid_SM_elements.size(); i++){
		max_chiralities = TMath::Max( (int)max_chiralities, (int)grid_SM_elements[i].size());
	}
	
	TH3D* h_relative_angular_distribution[max_chiralities];
	TH3D* h_relative_angular_distribution_bin_norm[max_chiralities];
	for(unsigned int i = 0; i < max_chiralities; i++){
		h_relative_angular_distribution[i] = new TH3D( prints("h_relative_angular_distribution_%i",global_hist_count).c_str(),"",
		binning[0],range_min[0],range_max[0],
		binning[1],range_min[1],range_max[1],
		binning[2],range_min[2],range_max[2]
		);
		h_relative_angular_distribution_bin_norm[i] = new TH3D( prints("h_relative_angular_distribution_bin_norm_%i",global_hist_count).c_str(),"",
		binning[0],range_min[0],range_max[0],
		binning[1],range_min[1],range_max[1],
		binning[2],range_min[2],range_max[2]
		);
		global_hist_count++;
	}
	
	for(unsigned int i = 0; i < grid_SM_elements.size(); i++){
		for(unsigned int j = 0; j < grid_SM_elements[i].size(); j++){
			h_relative_angular_distribution[j]->Fill( grid_variables[i][variable_index[0]], grid_variables[i][variable_index[1]], grid_variables[i][variable_index[2]], grid_SM_elements[i][j]);
			//h_relative_angular_distribution[j]->Fill(grid_variables[i][variable_index[0]], grid_variables[i][variable_index[1]], grid_variables[i][variable_index[2]]);
			h_relative_angular_distribution_bin_norm[j]->Fill( grid_variables[i][variable_index[0]], grid_variables[i][variable_index[1]], grid_variables[i][variable_index[2]]);
		}
		
	}
	cout << "3D Entries \t";
	for(unsigned int i = 0; i < max_chiralities; i++){
		h_relative_angular_distribution[i]->Divide(h_relative_angular_distribution_bin_norm[i]);
		if(h_relative_angular_distribution[i]->Integral() != 0)
			h_relative_angular_distribution[i]->Scale( 1. / h_relative_angular_distribution[i]->Integral() );
		cout << h_relative_angular_distribution[i]->GetEntries() << "\t";
	}
	cout << endl;
	
	double integral[max_chiralities];
	for(unsigned int i = 0; i < max_chiralities; i++){
		integral[i] = 0;
	}
	
	
	for(unsigned int i = 0; i < binning[0]; i++){
	for(unsigned int j = 0; j < binning[1]; j++){
	for(unsigned int k = 0; k < binning[2]; k++){
		
		vector < double > current_center;
		vector < double > current_width;
		
		current_center.push_back( h_relative_angular_distribution[0]->GetXaxis()->GetBinCenter(i+1) );
		current_center.push_back( h_relative_angular_distribution[0]->GetYaxis()->GetBinCenter(j+1) );
		current_center.push_back( h_relative_angular_distribution[0]->GetZaxis()->GetBinCenter(k+1) );
		
		current_width.push_back( h_relative_angular_distribution[0]->GetXaxis()->GetBinWidth(i+1) );
		current_width.push_back( h_relative_angular_distribution[0]->GetYaxis()->GetBinWidth(j+1) );
		current_width.push_back( h_relative_angular_distribution[0]->GetZaxis()->GetBinWidth(k+1) );
	
		center.push_back(current_center);
		width.push_back(current_width);
		
		vector < double > chirality_content;
		for(unsigned int l = 0; l < max_chiralities; l++){
			chirality_content.push_back(  h_relative_angular_distribution[l]->GetBinContent(i+1,j+1,k+1) );
			integral[l] += h_relative_angular_distribution[l]->GetBinContent(i+1,j+1,k+1);
		}
		Distribution.push_back(chirality_content);
	}	
	}	
	}
	cout << "integral:\t";
	for(unsigned int i = 0; i < max_chiralities; i++){
		cout << integral[i] << "\t";
	}
	cout << endl;
	
}


void CalculateMultiAngularTGCDistribution(
vector < vector < vector < double > > > &coefficients,
vector < vector < double > > &center,
vector < vector < double > > &width,
vector < double > TGC_scale,
vector < string > chiral_label,
string process_label,
vector < vector < double > > grid_variables,
vector < vector < double > > grid_SM_elements,
vector < vector < vector < double > > > grid_R_elements,
vector<unsigned int> variable_index,
vector<unsigned int> binning,
vector<double> range_min,
vector<double> range_max
){
	
	coefficients.clear();
	center.clear();
	width.clear();
	
	unsigned int max_chiralities = TGC_scale.size();
	
	TH3D* h_angular_distribution[max_chiralities][10];
	for(unsigned int i = 0; i < max_chiralities; i++){
		for(unsigned int j = 0; j < 10; j++){
			h_angular_distribution[i][j] = new TH3D( prints("h_angular_distribution_%i",global_hist_count).c_str(),"",
			binning[0],range_min[0],range_max[0],
			binning[1],range_min[1],range_max[1],
			binning[2],range_min[2],range_max[2]
			);
			global_hist_count++;
		}
	}

	for(unsigned int i = 0; i < grid_SM_elements.size(); i++){
		for(unsigned int j = 0; j < grid_SM_elements[i].size(); j++){
			h_angular_distribution[j][0]->Fill( grid_variables[i][variable_index[0]], grid_variables[i][variable_index[1]], grid_variables[i][variable_index[2]], grid_SM_elements[i][j]);
		}
	}

	for(unsigned int i = 0; i < grid_R_elements.size(); i++){
		for(unsigned int j = 0; j < grid_R_elements[i].size(); j++){
			for(unsigned int k = 0; k < grid_R_elements[i][j].size(); k++){
				h_angular_distribution[k][j+1]->Fill( grid_variables[i][variable_index[0]], grid_variables[i][variable_index[1]], grid_variables[i][variable_index[2]], grid_R_elements[i][j][k]);
			}
		}
	}

	vector < vector < TH3D* > > h_relative_angular_distribution(10);
	for(unsigned int i = 0; i < max_chiralities; i++){
		for(unsigned int j = 0; j < 10; j++){
			if(h_angular_distribution[i][j]->Integral() != 0){
				//h_angular_distribution[i][j]->Scale( 1. / h_angular_distribution[i][j]->Integral() );
			}
			h_relative_angular_distribution[j].push_back(h_angular_distribution[i][j]);
		}	
	}

	vector < vector < TH3D* > > h_relative_TGC_distribution = generate_TGC_coefficients( TGC_scale, h_relative_angular_distribution, chiral_label, process_label );

	for(unsigned int i = 0; i < binning[0]; i++){
	for(unsigned int j = 0; j < binning[1]; j++){
	for(unsigned int k = 0; k < binning[2]; k++){
		
		vector < double > current_center;
		vector < double > current_width;
		
		current_center.push_back( h_relative_TGC_distribution[0][0]->GetXaxis()->GetBinCenter(i+1) );
		current_center.push_back( h_relative_TGC_distribution[0][0]->GetYaxis()->GetBinCenter(j+1) );
		current_center.push_back( h_relative_TGC_distribution[0][0]->GetZaxis()->GetBinCenter(k+1) );
		
		current_width.push_back( h_relative_TGC_distribution[0][0]->GetXaxis()->GetBinWidth(i+1) );
		current_width.push_back( h_relative_TGC_distribution[0][0]->GetYaxis()->GetBinWidth(j+1) );
		current_width.push_back( h_relative_TGC_distribution[0][0]->GetZaxis()->GetBinWidth(k+1) );
	
		center.push_back(current_center);
		width.push_back(current_width);
		
		vector < vector < double > > coefficient_content;
		for(unsigned int m = 0; m < 9; m++){
			vector < double > chirality_content;
			for(unsigned int l = 0; l < max_chiralities; l++){
				chirality_content.push_back(  h_relative_TGC_distribution[m+1][l]->GetBinContent(i+1,j+1,k+1) );
			}
			coefficient_content.push_back(chirality_content);
		}	
		coefficients.push_back(coefficient_content);

	}	
	}	
	}
	
}

void CreateMultiProcessMultiAngularDistributions(){
	
  const string EnergyCMS = "250GeV";
  const bool Use_Generated_Hist_binning_range = true;
  string general_File_path = output_dir + "/grids_root/";
  
	TCanvas* c;
	string c_name;
	
	vector < string > output_file_descriptions {};
	vector < vector < string > > input_file_descriptions {};
  vector < string > file_description {};
  vector < unsigned int > numbervariables {};
  vector < unsigned int > numberchiralities {};
  vector < unsigned int > MaxNumberRelements {};
  vector < string > TGCdeviation {};
  vector<vector<unsigned int>> variable_index {};
  vector<vector<unsigned int>> binning {};
  vector<vector<double>> range_min {};
  vector<vector<double>> range_max {};
  
  
  if (process == "ww_sl0muq") {
    output_file_descriptions.push_back("WW_semilep_MuAntiNu");
    output_file_descriptions.push_back("WW_semilep_AntiMuNu");
    
    input_file_descriptions.push_back({"grid_ww_sl0muq_leptonic"});
    input_file_descriptions.push_back({"grid_ww_sl0muq_hadronic"});
    
    numbervariables.push_back(5);
    numbervariables.push_back(5);
    
    numberchiralities.push_back(4);
    numberchiralities.push_back(4);
    
    MaxNumberRelements.push_back(9);
    MaxNumberRelements.push_back(9);
    
    TGCdeviation.push_back("0.0001");
    TGCdeviation.push_back("0.0001");
    
    variable_index.push_back({0,1,2});
    variable_index.push_back({0,1,2});
    
    binning.push_back({20,10,10});
    binning.push_back({20,10,10});
    
    range_min.push_back( { -1+(1./(double)(binning[0][0])), -1+(1./(double)(binning[0][1])), -TMath::Pi() / ( (double)(binning[0][2]) )} );
    range_min.push_back( { -1+(1./(double)(binning[1][0])), -1+(1./(double)(binning[1][1])), -TMath::Pi() / ( (double)(binning[1][2]) )} );
        
    range_max.push_back( { 1+(1./(double)(binning[0][0])), 1+(1./(double)(binning[0][1])), ( 2. * TMath::Pi() ) + range_min[0][2] } );
    range_max.push_back( { 1+(1./(double)(binning[1][0])), 1+(1./(double)(binning[1][1])), ( 2. * TMath::Pi() ) + range_min[1][2] } );
    
  } else if (process == "sw_sl0qq") {
    output_file_descriptions.push_back("sW_semilep_eMinus");
    output_file_descriptions.push_back("sW_semilep_ePlus");
    
    input_file_descriptions.push_back({"grid_sw_sl0qq_plus"});
    input_file_descriptions.push_back({"grid_sw_sl0qq_minus"});
    
    numbervariables.push_back(6);
    numbervariables.push_back(6);

    numberchiralities.push_back(4);
    numberchiralities.push_back(4);
    
    MaxNumberRelements.push_back(9);
    MaxNumberRelements.push_back(9);
    
    TGCdeviation.push_back("0.0001");
    TGCdeviation.push_back("0.0001");
    
    variable_index.push_back({0,1,3}); // Indices: cosThetaW, cosThetae, Menu
    variable_index.push_back({0,1,3}); // Indices: cosThetaW, cosThetae, Menu
    
    binning.push_back({20,10,20});
    binning.push_back({20,10,20});
    
    range_min.push_back( { -1+(1./(double)(binning[0][0])), -1+(1./(double)(binning[0][1])), 0.0 } );
    range_min.push_back( { -1+(1./(double)(binning[1][0])), -1+(1./(double)(binning[1][1])), 0.0 } );
        
    range_max.push_back( { 1+(1./(double)(binning[0][0])), 1+(1./(double)(binning[0][1])), 250.-80.41 } ); // 80.41 = mW in O'Mega, 250 = COM energy
    range_max.push_back( { 1+(1./(double)(binning[1][0])), 1+(1./(double)(binning[1][1])), 250.-80.41 } ); // 80.41 = mW in O'Mega, 250 = COM energy
  } else {
    throw std::invalid_argument("Unknown process " + process);
  }
	
	
	if(output_file_descriptions.size() == input_file_descriptions.size()){
		for(unsigned int i = 0; i < output_file_descriptions.size(); i++){
			
			cout << output_file_descriptions[i] << endl;
			vector < string > chiral_label;
			chiral_label.push_back("eLpL");
			chiral_label.push_back("eLpR");
			chiral_label.push_back("eRpL");
			chiral_label.push_back("eRpR");	
			
			vector < vector < double > > grid_variables;
			vector < vector < double > > grid_SM_elements;
			vector < vector < vector < double > > > grid_R_elements;
			vector < vector < double > > Distribution;
			vector < vector < double > > center;
			vector < vector < double > > width;
			
			vector < vector < double > > scaling(input_file_descriptions[i].size());
			for(unsigned int j = 0; j < scaling.size(); j++){
				for(unsigned int k = 0; k < numberchiralities[i]; k++){
					scaling[j].push_back(1);
				}
			}
			
			vector < string > current_file_names;
			for(unsigned int j = 0; j < input_file_descriptions[i].size(); j++){
				current_file_names.push_back( general_File_path + input_file_descriptions[i][j] );
			}
			
			read_grid_files(grid_variables, grid_SM_elements, grid_R_elements, numbervariables[i], numberchiralities[i], MaxNumberRelements[i], current_file_names, scaling, EnergyCMS, TGCdeviation[i], true );				
			
			vector<unsigned int> Genereated_binning {0,0,0};
			vector<double> Genereated_range_min {0,0,0};
			vector<double> Genereated_range_max {0,0,0};
			
			GeneratedHistConfigMultiAngularDistribution(grid_variables,variable_index[i],Genereated_binning,Genereated_range_min,Genereated_range_max);
			cout << "Created Hist Init:" << endl; 
			for(unsigned int i = 0; i < 3; i++){
				cout << Genereated_binning[i] << "\t";
				cout << Genereated_range_min[i] << "\t";
				cout << Genereated_range_max[i] << "\t";
				cout << endl;
			}
			
			cout << grid_SM_elements.size() << endl;

			vector < vector < double > > grid_SM_elements_transposed;
			vector < double > grid_SM_elements_for_range;

			
			for(unsigned int n = 0; n < grid_SM_elements.size(); n++ ){
				if(grid_SM_elements[n].size() > grid_SM_elements_transposed.size()){
					grid_SM_elements_transposed.resize(grid_SM_elements[n].size());
				}
				for(unsigned int m = 0; m < grid_SM_elements[n].size(); m++){
					grid_SM_elements_transposed[m].push_back(grid_SM_elements[n][m]);
					grid_SM_elements_for_range.push_back(grid_SM_elements[n][m]);
				}
				
			}
			
			unsigned int sm_value_binning = 1000;

			const bool sm_value_distribution_log_mode = true;
			
			double sm_value_min = GetL0Minimum(grid_SM_elements_for_range);
			double sm_value_max = GetMaximum(grid_SM_elements_for_range);
			
			grid_SM_elements_for_range.clear();

			if(sm_value_distribution_log_mode){
				sm_value_min = TMath::Log10(sm_value_min);
				sm_value_max = TMath::Log10(sm_value_max);
			}
			if(sm_value_max <= sm_value_min){
				sm_value_min = 0;
				sm_value_max = 0;
			}
			
			//double sm_value_bin_width = (sm_value_max-sm_value_min)/(double)(2*sm_value_binning);
			//sm_value_min -= 0.1*sm_value_bin_width;
			//sm_value_max += 0.1*sm_value_bin_width;
					
			cout << "[" << sm_value_min << ", " << sm_value_max << "]" << "\t" << sm_value_binning << endl;
					
			
			vector < TH1D* > h_SM_value_distribution(chiral_label.size());
			vector < double > SM_value_distribution_values_max;
			vector < double > SM_value_distribution_values_min;
			for(unsigned int m = 0; m < grid_SM_elements_transposed.size(); m++){
				if( m < chiral_label.size() ){
					
					h_SM_value_distribution[m] = new TH1D(prints("h_SM_value_distribution_%s_%s",output_file_descriptions[i].c_str(),chiral_label[m].c_str()).c_str(),"",sm_value_binning,sm_value_min,sm_value_max);
					h_SM_value_distribution[m]->SetLineColor(color_beamer_display[m]);
					h_SM_value_distribution[m]->SetLineWidth(2);
					h_SM_value_distribution[m]->GetXaxis()->SetTitle("Log_{10}#(){d#sigma / d#Omega}");
					h_SM_value_distribution[m]->GetYaxis()->SetTitle("entries");
					
					for(unsigned int n = 0; n < grid_SM_elements_transposed[m].size(); n++ ){
						double current_sm_value_of_distribution = grid_SM_elements_transposed[m][n];
						if(sm_value_distribution_log_mode)
							current_sm_value_of_distribution = TMath::Log10(current_sm_value_of_distribution);
						h_SM_value_distribution[m]->Fill(current_sm_value_of_distribution);
					}
					
					vector < double > non_zero_bin_entries;
					for(int j = 0; j < h_SM_value_distribution[m]->GetNbinsX(); j++){
						if( h_SM_value_distribution[m]->GetBinContent(j+1) > 0){
							non_zero_bin_entries.push_back(h_SM_value_distribution[m]->GetBinContent(j+1));
						}
					}
					
					SM_value_distribution_values_min.push_back( GetMinimum( non_zero_bin_entries ) );
					SM_value_distribution_values_max.push_back( GetMaximum( non_zero_bin_entries ) );
					
				}else{
					cout << "too many chiralities!!!" << endl;
				}
			}
			
			double SM_value_distribution_max = GetMaximum(SM_value_distribution_values_max);
			double SM_value_distribution_min = GetMinimum(SM_value_distribution_values_min);
			
			SM_value_distribution_max += 0.1*TMath::Abs(SM_value_distribution_max);
			SM_value_distribution_min -= 0.1*TMath::Abs(SM_value_distribution_min);
			
			TLegend* leg_SM_value_distribution = new TLegend(0.17,0.93,0.95,0.99);
			//leg_SM_value_distribution->SetNColumns(chiral_label.size());
			c_name = "SM_Value_Distribution_" + output_file_descriptions[i];
			c = new TCanvas(c_name.c_str(),c_name.c_str(),1000,700);
			c->cd();
			gPad->SetLogy(true);
			
			bool first_hist_of_loop_h_SM_value_distribution = true;
			int h_SM_value_distribution_display_count = 0;
			for(unsigned int j = 0; j < h_SM_value_distribution.size(); j++){
				if(h_SM_value_distribution[j]->GetMaximum() > 0){
					h_SM_value_distribution_display_count++;
					h_SM_value_distribution[j]->SetMinimum(SM_value_distribution_min);
					h_SM_value_distribution[j]->SetMaximum(SM_value_distribution_max);
					leg_SM_value_distribution->AddEntry(h_SM_value_distribution[j],chiral_label[j].c_str(),"L");
					if(first_hist_of_loop_h_SM_value_distribution){
						first_hist_of_loop_h_SM_value_distribution = false;
						h_SM_value_distribution[j]->Draw();
					}else
						h_SM_value_distribution[j]->Draw("SAME");
				}
			
			}
			leg_SM_value_distribution->SetNColumns(h_SM_value_distribution_display_count);
			leg_SM_value_distribution->Draw();
			
			c->SaveAs( (path_plots + c_name + ".pdf").c_str() );
			
			if(Use_Generated_Hist_binning_range){
				CalculateMultiAngularDistribution(Distribution,center,width,grid_variables,grid_SM_elements,variable_index[i],Genereated_binning,Genereated_range_min,Genereated_range_max);
			}else{
				CalculateMultiAngularDistribution(Distribution,center,width,grid_variables,grid_SM_elements,variable_index[i],binning[i],range_min[i],range_max[i]);
			}
			unsigned int size_min = center[0].size();
			unsigned int size_max = center[0].size();
			
			double value_min[3] = {center[0][0],center[0][1],center[0][3]};
			double value_max[3] = {center[0][0],center[0][1],center[0][3]};
			
			vector < TH3D* > h_relative_angular_distribution_current(chiral_label.size());
			
			for(unsigned int l = 0; l < h_relative_angular_distribution_current.size(); l++){
				if(Use_Generated_Hist_binning_range){
					h_relative_angular_distribution_current[l] = new TH3D( prints("h_relative_angular_distribution_%s_%s",output_file_descriptions[i].c_str(),chiral_label[l].c_str()).c_str(),"",
						Genereated_binning[0],Genereated_range_min[0],Genereated_range_max[0],
						Genereated_binning[1],Genereated_range_min[1],Genereated_range_max[1],
						Genereated_binning[2],Genereated_range_min[2],Genereated_range_max[2]
					);
				}else{
					h_relative_angular_distribution_current[l] = new TH3D( prints("h_relative_angular_distribution_%s_%s",output_file_descriptions[i].c_str(),chiral_label[l].c_str()).c_str(),"",
						binning[i][0],range_min[i][0],range_max[i][0],
						binning[i][1],range_min[i][1],range_max[i][1],
						binning[i][2],range_min[i][2],range_max[i][2]
					);
				}
			}
			
			for(unsigned int k = 1; k < center.size(); k++){
				if(size_min > center[k].size()){
					size_min = center[k].size();
				}
				if(size_max < center[k].size()){
					size_max = center[k].size();
				}
				
				for(unsigned int j = 0; j < 3; j++){
					if(value_min[j] > center[k][j]){
						value_min[j] = center[k][j];
					}
					if(value_max[j] < center[k][j]){
						value_max[j] = center[k][j];
					}
				}
				
				if(center[i].size() >= 3){
					for(unsigned int l = 0; l < h_relative_angular_distribution_current.size(); l++){
						int current_bin_X = h_relative_angular_distribution_current[l]->GetXaxis()->FindBin(center[k][0]);
						int current_bin_Y = h_relative_angular_distribution_current[l]->GetYaxis()->FindBin(center[k][1]);
						int current_bin_Z = h_relative_angular_distribution_current[l]->GetZaxis()->FindBin(center[k][2]);
						if(h_relative_angular_distribution_current[l]->GetBinContent(current_bin_X,current_bin_Y,current_bin_Z) == 0.){
							h_relative_angular_distribution_current[l]->SetBinContent(current_bin_X,current_bin_Y,current_bin_Z,Distribution[k][l]);
						}else{
							cout << "Bin " << current_bin_X << "," << current_bin_Y << "," << current_bin_Z << " already used" << endl;
						}
					}
				}
			}	
			for(unsigned int j = 0; j < 3; j++){
				cout << "[" << value_min[j] << ", " << value_max[j] << "];\t";
			}
			cout << endl;
			cout << size_min << "\t" << size_max << endl;
			
			
			vector < vector < TH1D* > > h_relative_angular_distribution_projections(3);
			for(unsigned int l = 0; l < h_relative_angular_distribution_current.size(); l++){
				h_relative_angular_distribution_projections[0].push_back(
					h_relative_angular_distribution_current[l]->ProjectionX( prints("h_relative_angular_distribution_%s_%s_px",output_file_descriptions[i].c_str(),chiral_label[l].c_str()).c_str() )
				);
				h_relative_angular_distribution_projections[1].push_back(
					h_relative_angular_distribution_current[l]->ProjectionY( prints("h_relative_angular_distribution_%s_%s_py",output_file_descriptions[i].c_str(),chiral_label[l].c_str()).c_str() )
				);
				h_relative_angular_distribution_projections[2].push_back(
					h_relative_angular_distribution_current[l]->ProjectionZ( prints("h_relative_angular_distribution_%s_%s_pz",output_file_descriptions[i].c_str(),chiral_label[l].c_str()).c_str() )
				);
			}
			
			for(unsigned int k = 0; k < 3; k++){
				c_name = prints("Relative_Angular_Distribution_%s_%i",output_file_descriptions[i].c_str(),k);
				c = new TCanvas(c_name.c_str(),c_name.c_str(),1000,700);
				c->cd();
				//gPad->SetLogy(true);
				double currrent_hist_min = h_relative_angular_distribution_projections[k][0]->GetMinimum();
				double currrent_hist_max = h_relative_angular_distribution_projections[k][0]->GetMaximum();
				for(unsigned int l = 1; l < h_relative_angular_distribution_projections[k].size(); l++){
					if(h_relative_angular_distribution_projections[k][l]->GetMinimum() > 0 && currrent_hist_min > h_relative_angular_distribution_projections[k][l]->GetMinimum()){
						currrent_hist_min = h_relative_angular_distribution_projections[k][l]->GetMinimum();
					}
					if(h_relative_angular_distribution_projections[k][l]->GetMaximum() > 0 && currrent_hist_max < h_relative_angular_distribution_projections[k][l]->GetMaximum()){
						currrent_hist_max = h_relative_angular_distribution_projections[k][l]->GetMaximum();
					}
				}
				if(currrent_hist_min - 0.1*TMath::Abs(currrent_hist_min) > 0){
					currrent_hist_min -= 0.1*TMath::Abs(currrent_hist_min);
				}
				currrent_hist_max += 0.1*TMath::Abs(currrent_hist_max);
				
				for(unsigned int l = 0; l < h_relative_angular_distribution_projections[k].size(); l++){
					//cout << h_relative_angular_distribution_projections[k][l]->GetEntries() << "\t";
					//cout << h_relative_angular_distribution_projections[k][l]->GetMinimum() << "\t";
					//cout << h_relative_angular_distribution_projections[k][l]->GetMaximum() << endl;					
					h_relative_angular_distribution_projections[k][l]->SetMinimum(currrent_hist_min);
					h_relative_angular_distribution_projections[k][l]->SetMaximum(currrent_hist_max); 
					h_relative_angular_distribution_projections[k][l]->SetLineColor(color_beamer_display[l]);
					if(l==0)
						h_relative_angular_distribution_projections[k][l]->Draw("Hist");
					else
						h_relative_angular_distribution_projections[k][l]->Draw("HistSame");
				}
			}
				
			
			write_angular_file( ( path_angular_distributions + "AngularDistribution_" + EnergyCMS + "_" + output_file_descriptions[i] + ".root" ),
				"AngularDistribution",Distribution,center,width,chiral_label);
						
			if(MaxNumberRelements[i]>0){
				vector < double > TGC_scale;
				vector < vector < vector < double > > > coefficients;
				for(unsigned int j = 0; j < numberchiralities[i]; j++){
					TGC_scale.push_back( stod( TGCdeviation[i] ) ); 
				}
				
				if(Use_Generated_Hist_binning_range){
					CalculateMultiAngularTGCDistribution(
						coefficients,center,width,TGC_scale,chiral_label,output_file_descriptions[i],grid_variables,grid_SM_elements,grid_R_elements,
						variable_index[i],Genereated_binning,Genereated_range_min,Genereated_range_max
					);
				}else{
					CalculateMultiAngularTGCDistribution(
						coefficients,center,width,TGC_scale,chiral_label,output_file_descriptions[i],grid_variables,grid_SM_elements,grid_R_elements,variable_index[i],binning[i],range_min[i],range_max[i]
					);
				}
				string coefficient_labeling[] = {"A","B","C","D","E","F","G","H","I"};
				vector < string > PNPC_description;
				for(unsigned int j = 0; j < 9; j++){
					PNPC_description.push_back( ( "TGC" + coefficient_labeling[j] ) );
				}
				
				write_PNPC_file( ( path_PNPC_distributions + "AngularTGCDistribution_" + EnergyCMS + "_" + output_file_descriptions[i] + ".root" ),
					"AngularDistribution",coefficients,center,width,PNPC_description,chiral_label);
			}
			
		}
	}
}

void CreateAngularDistributions(){
	
	time_t start_creating = time(0);
	cout << start_creating << endl;
	logfile.open("log.txt");
	logfile.precision(16);
	logfile << "" << endl;

	CreateMultiProcessMultiAngularDistributions();
	
	cout << "done in " << time(0)-start_creating << "s" << endl;
}
