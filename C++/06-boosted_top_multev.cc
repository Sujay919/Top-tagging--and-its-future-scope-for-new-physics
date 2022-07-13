//----------------------------------------------------------------------
/// \file
/// \page Example13 13 - boosted top tagging
///
/// fastjet example program, illustration of carrying out boosted
/// top subjet ID analysis using the Johns Hopkins top tagger as
/// introduced in arXiv:0806.0848 (Kaplan, Rehermann, Schwartz
/// and Tweedie)
///
/// run it with    : ./13-boosted_top < data/boosted_top_event.dat
///
/// Source code: 13-boosted_top.cc
//----------------------------------------------------------------------


//STARTHEADER
// $Id$
//
// Copyright (c) 2005-2018, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER

#include <iostream> // needed for io
#include <sstream>  // needed for internal io 
#include <cmath>
#include <fstream>  // needed for file read write
#include <iomanip>  // for formatting output

#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/JHTopTagger.hh>

using namespace std;
using namespace fastjet;


//----------------------------------------------------------------------
// forward declaration for printing out info about a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream &, const PseudoJet &);
int jetform(vector<PseudoJet>,ofstream &);
int top=0,top1=0,total=0,total1=0;
int jetE[1000000],bin[12]={0,0,0,0,0,0,0,0,0,0,0,0},
              bint[12]={0,0,0,0,0,0,0,0,0,0,0,0};
//----------------------------------------------------------------------
// core of the program
//----------------------------------------------------------------------
int main(){

  vector<PseudoJet> particles;
  // read in data in format px py pz E b-tag [last of these is optional]
  // lines starting with "#" are considered as comments and discarded
  //----------------------------------------------------------
  string line;
  ifstream fin("Data/FASTJET-INPUT3_bck_ttbar_fastjet_new");
  ofstream fout("output/06-boosted_top_multev_histg.out");//, ios_base::app
  //ofstream hout("output/06-boosted_top_multev.out");
  while (getline(fin,line)) {
    if (line.substr(0,1) == "#") {continue;}
    istringstream linestream(line);
    double px,py,pz,E,pdg,n;
    linestream >> px >> py >> pz >> E >> pdg >> n;
    if(px==0 and py==0 and pz==0 and E==0 and pdg==0 and n==0){
      // total++;
      jetform(particles,fout);
      particles.clear();
    }
    else{
    // construct the particle
    particles.push_back(PseudoJet(px,py,pz,E));
    }
  }
   jetform(particles,fout);




  fout<<"#**********************FINAL RESULT*****************************"<<endl;
  fout<<"#NO of top quarks found in hardest jets:"<<top<<endl;
  fout<<"#Number of total events                :"<<total<<endl;
  float  percentage;
  if(total!=0){
    percentage = ((float)top/(float)total)*100;
  }
  else percentage = 0;
  
  fout<<"#percentage of top tagged:"<<percentage<<"%"<<endl<<"#";

  //bin data out
  fout<<setw(9)<<"#bin"<<setw(15)<<"total_event"<<setw(15)<<"tagged"<<setw(15)<<"percentage"<<endl;
  for(int i=0;i<=10;i++){
    float percentage;
    if(bin[i]!=0){
      percentage = ((float)bint[i]/(float)bin[i])*100;
    }
    else percentage = 0;
     
   fout<<setw(10)<<(i*100)<<setw(15)<<bin[i]<<setw(15)<<bint[i]<<setw(15)<<percentage<<endl; 
  }

  //cout<<endl<<endl<<"**********************FINAL RESULT*****************************"<<endl;
  //cout<<"NO of top quarks found in hardest jets:"<<top<<endl;
  //cout<<"Number of total valid events:"<<total<<endl;
  //cout<<"percentage of top tagged    :"<<percentage<<"%"<<endl;
  fout.close();
  fin.close();
}


//----------------------------------------------------------------------
// does the actual work for printing out a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream & ostr, const PseudoJet & jet) {
  ostr << "pt, y, phi =" << setprecision(6)
       << " " << setw(9) << jet.perp() 
       << " " << setw(9)  <<  jet.rap()  
       << " " << setw(9)  <<  jet.phi()  
       << ", mass = " << setw(9) << jet.m();
  return ostr;
}


int jetform(vector<PseudoJet> particles, ofstream & fout){
    //cout <<endl<<endl<<endl<< "New Event*************************************************************"<<endl;
  // total1++;
  // cout<< "this is event no" << total << endl;
  // compute the parameters to be used through the analysis
  // ----------------------------------------------------------
  double Et=0;
  for (unsigned int i=0; i<particles.size(); i++)
    Et += particles[i].perp();

  double R, delta_p, delta_r;
  if      (Et>2600){ R=0.4; delta_p=0.05; delta_r=0.19;}
  else if (Et>1600){ R=0.6; delta_p=0.05; delta_r=0.19;}
  else if (Et>1000){ R=0.8; delta_p=0.10; delta_r=0.19;}
  else{ /*cerr << "Et has to be at least 1 TeV"<< endl; */ return 1;}
  total++;
  cout<< "this is valid event no" << total << endl;
  double ptmin = min(500.0, 0.7*Et/2);

  // find the jets
  // ----------------------------------------------------------
  JetDefinition jet_def(cambridge_algorithm, R);
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
  jetE[total]=jets[0].pt();
  // cout << "Ran: " << jet_def.description() << endl << endl;
  // cout << "2 Hardest jets: " << jets[0] << endl
  //      << "                " << jets[1] << endl << endl;

  if (jets[0].perp()<ptmin){
    //cout << "No jet above the ptmin threshold" << endl;
    return 2;
  }

  // now do jet tagging using the Johns Hopkins top tagger
  // For simplicity, we just apply it to the hardest jet.
  //
  // In addition to delta_p and delta_r, note that there are two
  // further parameters to the JH top tagger that here are implicitly
  // set to their defaults:
  //
  // - cos_theta_W_max (defaults to 0.7) 
  // - mW (defaults to 80.4). 
  //
  // The value for mW implicitly assumes that momenta are passed in
  // GeV.
  // ----------------------------------------------------------
  JHTopTagger top_tagger(delta_p, delta_r);
  top_tagger.set_top_selector(SelectorMassRange(150,200));
  top_tagger.set_W_selector  (SelectorMassRange( 65, 95));

  PseudoJet tagged = top_tagger(jets[0]);
  PseudoJet tagged1 = top_tagger(jets[1]);
  //cout << "Ran the following top tagger: " << top_tagger.description() << endl;

  //bin creation
  for(int i=0;i<=1000;i=i+100){
    if(i==0){
     bin[i]=0;
     bint[i]=0;
    }
    else{
      if((jetE[total]>(i-100)) and (jetE[total]<i)){
        // fout<<bin[i/100];
        bin[i/100]++;
        // fout<<bin[i/100];
        if(tagged!=0){
        bint[i/100]++;
        }
      }
      else continue;
    }
  }
  // }
  // else if (jetE[total]>100 and jetE[total]<=200){
  //   bin[2]++;
  //   if(tagged!=0){
  //   bint[2]++;
  //   }
  // }
  // else if (jetE[total]>200 and jetE[total]<=300){
  //   bin[3]++;
  //   if(tagged!=0){
  //   bint[3]++;
  //   }
  // }
  // else if (jetE[total]>300 and jetE[total]<=400){
  //   bin[4]++;
  //   if(tagged!=0){
  //   bint[4]++;
  //   }
  // }
  // else if (jetE[total]>400 and jetE[total]<=500){
  //   bin[5]++;
  //   if(tagged!=0){
  //   bint[5]++;
  //   }
  // }
  // else if (jetE[total]>500 and jetE[total]<=600){
  //   bin[6]++;
  //   if(tagged!=0){
  //   bint[6]++;
  //   }
  // }
  // else if (jetE[total]>600 and jetE[total]<=700){
  //   bin[7]++;
  //   if(tagged!=0){
  //   bint[7]++;
  //   }
  // }
  // else if (jetE[total]>700 and jetE[total]<=800){
  //   bin[8]++;
  //   if(tagged!=0){
  //   bint[8]++;
  //   }
  // }
  // else if (jetE[total]>800 and jetE[total]<=900){
  //   bin[9]++;
  //   if(tagged!=0){
  //   bint[9]++;
  //   }
  // }
  // else if (jetE[total]>900 and jetE[total]<=1000){
  //   bin[10]++;
  //   if(tagged!=0){
  //   bint[10]++;
  //   }
  // }

  if (tagged == 0 ){
    // cout << "No top substructure found" << endl;
    // cout << "No top substructure found" << endl;
    return 0;
  }
  top++;
  if (tagged1 == 0){
    return 0;
  }

  //cout << "Found top substructure from the hardest jet:" << endl;
  //cout << "Found";
  top1++;
  //cout << "  top candidate:     " << tagged << endl;
  //cout << "  |_ W   candidate:  " << tagged.structure_of<JHTopTagger>().W() << endl;
  //cout << "  |  |_  W subjet 1: " << tagged.structure_of<JHTopTagger>().W1() << endl;
  //cout << "  |  |_  W subjet 2: " << tagged.structure_of<JHTopTagger>().W2() << endl;
  //cout << "  |  cos(theta_W) =  " << tagged.structure_of<JHTopTagger>().cos_theta_W() << endl;
  //cout << "  |_ non-W subjet:   " << tagged.structure_of<JHTopTagger>().non_W() << endl;

return 0;
}