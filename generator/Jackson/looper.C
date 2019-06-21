//--------------------LHEtoROOT-----------------------------------------------
//Author: Dan Phan (dan@umail.ucsb.edu)
//Date: December 2014
//Description/Instructions: Go to line 22, change name of file to your LHE file
//Change tree name to your desired root tree name
//Change name of output file to your desired output filename
//----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;

//Parameters
char* dirname = "subjobs";
char* filename  = "EventOutput.Real.00000001.lhe";
char* outputName = "GiBUU_Output";
char* treeName = "Tree";

double radius = 2.0*1.0;
double kF = 0.325*1.0;

int npos = -1;
int pCode = 2212;
int nCode = 2112;

//fill tree for only maxEvents number of events
int numEvents = 0; //initializing counter of events
int maxEvents = 1e9; //if you want to fill all events, make maxEvents huge


int looper(){


  //Declare TTree and TFile
  TFile *file = new TFile(Form("%s.root", outputName), "RECREATE");
  TTree *tree = new TTree("tree", Form("%s",treeName));  //tree called "tree"
  
  //Declare variables that will be stored in tree
  vector <int> pdgID;
  vector <TLorentzVector> four_position;
  vector <TLorentzVector> four_momentum;  //encodes 4-momentum and position

  int nParticles;      //variables in row1, first row beneath "<event>"
  
  //Match up variable with branch
  
  //variables to be filled in second_line_below loop
  tree->Branch("pdgID", &pdgID); 
  tree->Branch("four_position", &four_position);
  tree->Branch("four_momentum", &four_momentum);
  
  //filled separate from first_below_event and second_line_below loops

  //filled in first_below_event loop
  tree->Branch("nParticles",&nParticles);
 
  string line; //line of LHE file, to be looped over in main part of looper()

  bool first_line_below = false; //whether on line right below event
  bool second_line_below = false; //between <event> and </event>, is turned to true one line after "first_line_below", both turn off after filling tree

  //declare vectors used to fill other variables
  vector<float> row;  //row used to fill variables in second_line_below loop
  vector<float> row1; //row used to fill variables in first_line_below part of looper() 

  for (int i = 0; i < 200 ; i++)
    {
      cout << "Doing file " << i << "...\n";
      char istr[5];
      snprintf(istr,5,"%d",i);
      std::string dirstr = dirname;
      std::string filestr = filename;
      std::string s = dirstr + "/" + istr + "/" + filestr;

  //Opening data file
      fstream myfile (s.c_str(), ios_base::in);

  while (getline(myfile,line)){  //loops through file and fills line(the variable) with that line of the file

    char * cline = line;

    if (numEvents > maxEvents) break;  //fills tree only for maxEvents number of events

    //if on 2nd line below "<event>", fill variables(row1 and reweight entries filled in other loops)
    if (second_line_below == true) {   
      do {
        if (strstr(line.c_str(), "#") != NULL) break; //if line contains '#', line does not contain data anymore, break out of loop
        if (strstr(line.c_str(), "</event>") != NULL)
	  { 

      tree->Fill();   //only fill tree if vectors are nonempty(pdgID picked arbitrarily)
      numEvents++;
  
      //clear vectors so we can fill them again for next event
      pdgID.clear();  
      four_position.clear();
      four_momentum.clear();

	    break;
	  }
        //iss contains line
        istringstream iss;   //opening string stream
        iss.str(line);      //copying line into stream

        float val; 
         
     	while (iss >> val) row.push_back(val);  //fill row with data from line

        //Now split row into vectors for each variable
	int pID = row[0];
        TLorentzVector four_position_temp;
        four_position_temp.SetXYZT(row[2],row[3],row[4],row[5]);
        TLorentzVector four_momentum_temp;
        four_momentum_temp.SetPx(row[6]);
        four_momentum_temp.SetPy(row[7]);
        four_momentum_temp.SetPz(row[8]);
        four_momentum_temp.SetE(row[9]);
	double R = four_position_temp.P();
	double P = four_momentum_temp.P();
	if (R > radius)
	  {
	    if (P > kF)
	      {
		if (pID == pCode)
		  {
		    pdgID.push_back(pID); 
		    four_position.push_back(four_position_temp);
		    four_momentum.push_back(four_momentum_temp);
		  }
		if (pID == nCode)
		  {
		    pdgID.push_back(pID); 
		    four_position.push_back(four_position_temp);
		    four_momentum.push_back(four_momentum_temp);
		  }
	      }
	  }
        row.clear();  //clears row so we can fill it again with next line of data

    } while (getline(myfile,line)); 
    //and (strstr(line.c_str(), '</event>') == NULL));

      first_line_below = false;
      second_line_below = false;  // turn first_line_below, second_line_below off so we know we're done with those loops
    }
 
    //Line right below <event>(row1 variables)
    if (first_line_below == true) {
      second_line_below = true; //sets up loop over second_line_below variables

      istringstream iss;
      iss.str(line);

      float val1;

      while (iss >> val1) row1.push_back(val1);  //fill row1 with numbers in line

      nParticles = row1[0];
      row1.clear(); //clears row1 for next event 
      
    }
    
    //Check every line after </event> to see if we've reached the next event 
    if ((line.find("<event>") != -1)) {  //if line contains  "<event>", we turn first_line_below to true(begins filling row1 variables)

      first_line_below = true; //begins process to move onto next event
    }
    
  }
  myfile.close();
  cout << numEvents << "\n";
}
  tree->Fill(); //last event isn't filled, because there's no next "<event>" to trigger it
  numEvents++;
  file->cd();
  tree->Write(); 
  file->Close();

  return 0;
}
