//----------------------------------------------------------------------
/// \file
/// \page Example01 01 - basic usage example
///
/// fastjet basic example program:
///   simplest illustration of the usage of the basic classes:
///   fastjet::PseudoJet, fastjet::JetDefinition and
///   fastjet::ClusterSequence
/// compile it with : g++ 01-basic_multev.cc -o 01-basic_multev `../fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`
/// run it with     : ./01-basic_multev < ../fastjet-3.4.0/example/data/single-event.dat
///
/// Source code: 01-basic.cc
//----------------------------------------------------------------------

// STARTHEADER
//  $Id$
//
//  Copyright (c) 2005-2018, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
//  This file is part of FastJet.
//
//   FastJet is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//
//   The algorithms that underlie FastJet have required considerable
//   development and are described in hep-ph/0512210. If you use
//   FastJet as part of work towards a scientific publication, please
//   include a citation to the FastJet paper.
//
//   FastJet is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
// ENDHEADER

#include "fastjet/ClusterSequence.hh"
#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <fstream>  // needed for file read write
#include <iomanip>  // for formatting output

using namespace std;
int jetform(fastjet::JetAlgorithm jet_algorithm,
            vector<fastjet::PseudoJet> input_particles, double, ofstream &, ofstream &);
int algow(int f, vector<fastjet::PseudoJet> input_particles);

// int i;
double R[30];
 
 

/// an example program showing how to use fastjet
int main()
{

  // read in input particles
  //----------------------------------------------------------
  vector<fastjet::PseudoJet> input_particles;

  // input_particles is the (object of type)vector of initial particles of any type ( PseudoJet , HepLorentzVector ,
  // etc.) that can be used to initialise a PseudoJet and jet def contains the full specification of the
  // clustering (see Section 3.2).
  // https://www.geeksforgeeks.org/vector-in-cpp-stl/
  // So input_particles object can store series of pseudojets. Where if you call the address of first pseudojet you
  // you can do that by either i=input_particles.begin(stores the address of first jet) Or directly by using array
  // pointer where input_particle is treated like a array which stores array of jets inside it, So i=0 would point
  // to the address of first jet.
  ifstream fin("Data/FASTJET-INPUT3_bck_ttbar_fastjet_new");
  // ofstreame opening
  //ofstream fout("output/05-constituents_hardest_itaphi_multev.dat" );
  double px, py, pz, E, Etot = 0, ETtot = 0, pdg, n;
  valarray<double> fourvec(4);
  int index=0, f = 0;
  string base("05-constituents_hardest_itaphi_multev.dat");
  string loc("output/05-constituents_hardest_itaphi_multev/");
  // string base("/mafd");
  while (fin >> fourvec[0] >> fourvec[1] >> fourvec[2] >> fourvec[3] >> pdg >> n) {
    //cout<<fourvec[0]<<endl;
    if(fourvec[0]==0 and fourvec[1]==0 and fourvec[2]==0 and fourvec[3]==0 and pdg==0 and n==0){
        double Et=0; 
        for (unsigned int i=0; i<input_particles.size(); i++)
          Et+= input_particles[i].perp();
        if(Et<=1000) {
          input_particles.clear();
          index=0;
          }
        else{
          f++;
          cout<<f<<endl;
          algow(f,input_particles);
          input_particles.clear();
          index=0;
          }
        
        }
    else{
        fastjet::PseudoJet particle(fourvec);
        particle.set_user_index(index);
        //pushing new creeated pseudojets into array of addresses from end.
        input_particles.push_back(particle); 
        index++;
        }
     }

    // f++;
    algow(f,input_particles);
    input_particles.clear();
    index=0;









  fin.close();
  //closing output file
  //fout.close();
  
  return 0;
}

int jetform(fastjet::JetAlgorithm jet_algorithm,
            vector<fastjet::PseudoJet> input_particles, double R, ofstream &fout, ofstream &gout)
{
 


  // create a jet definition:
  // a jet algorithm with a given radius parameter
  //----------------------------------------------------------
  // double R = 0.6;

  fastjet::JetDefinition jet_def(jet_algorithm, R);

  // run the jet clustering with the above jet definition
  //----------------------------------------------------------
  fastjet::ClusterSequence clust_seq(input_particles, jet_def);

  // get the resulting jets ordered in pt
  //----------------------------------------------------------
  double ptmin = 5.0;
  vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));

  
  
  // fout<< "#Value of R is "<<R<<endl;
  // label the columns and print the details of hardest jet
    // fout<<"#"<<setw(5)<<"jet #"<<setw(15)<<"rapidity"<<setw(15)<<"phi"<<setw(15)
    // <<"pt"<<setw(15)<<"n constituents"<<endl;
 

  // get the constituents of the hardest jet
  vector<fastjet::PseudoJet> constituents = inclusive_jets[0].constituents();

  // fout<<"#"<<setw(5)<<0<<setw(15)<<inclusive_jets[0].rap()<<setw(15)<<inclusive_jets[0].phi()
  // <<setw(15)<<inclusive_jets[0].perp()<<setw(8)<<(unsigned int) constituents.size()<<endl;

  //lebel the column and print the ita and phi of constituents
  // fout<<setw(8)<<"#indices"<<setw(15)<< "rapidity(ita)"
  //       <<setw(15)<< "phi"<<setw(15)<<"pt"<<endl;
  for (unsigned int j=0; j<constituents.size(); j++){
      fout<<setw(8)<< constituents[j].user_index()<<setw(15)<< constituents[j].rap()
        <<setw(15)<< constituents[j].phi()<<setw(15)<<constituents[j].perp()<<endl;
    
  }
  constituents.clear();

   // get the constituents of the next hardest jet
  constituents = inclusive_jets[1].constituents();

  // fout<<"#"<<setw(5)<<0<<setw(15)<<inclusive_jets[0].rap()<<setw(15)<<inclusive_jets[0].phi()
  // <<setw(15)<<inclusive_jets[0].perp()<<setw(8)<<(unsigned int) constituents.size()<<endl;

  //lebel the column and print the ita and phi of constituents
  // fout<<setw(8)<<"#indices"<<setw(15)<< "rapidity(ita)"
  //       <<setw(15)<< "phi"<<setw(15)<<"pt"<<endl;
  for (unsigned int j=0; j<constituents.size(); j++){
      gout<<setw(8)<< constituents[j].user_index()<<setw(15)<< constituents[j].rap()
        <<setw(15)<< constituents[j].phi()<<setw(15)<<constituents[j].perp()<<endl;
    
  }
  
  
    // fout<<'$'<<endl;
  
  return 0;
}


int algow(int f, vector<fastjet::PseudoJet> input_particles ){
  cout<< "inside";
      cout<<f<<endl;
  // kt algorithm delta r and no of jets
    string base("05-constituents_hardest_itaphi_multev.dat");
    string loc("output/05-constituents_hardest_itaphi_multev/");
    string algo[] = {"kt","akt","cam"};
    fastjet::JetAlgorithm Algorithm[] ={
                                          fastjet::kt_algorithm,
                                          fastjet::cambridge_algorithm,
                                          fastjet::antikt_algorithm
                                        };
    for(int j=0;j<3;j++){
      R[0] = 0.6;
      ostringstream r;
      r<<R[0];
      string r_i=r.str();
      ofstream fout(loc+"H0/"+to_string(f)+"_"+r_i+"_"+algo[j]+"_"+base);
      ofstream gout(loc+"H1/"+to_string(f)+"_"+r_i+"_"+algo[j]+"_"+base);
      r_i.clear();
      r.clear();
      int i=0;
      jetform(Algorithm[j], input_particles, R[0], fout, gout);
      fout.close();
      gout.close();
      while (R[i] < 0.6){
        cout<<"inin";
        i++;
        R[i] = 0.2 * i+ 0.6;
        r<<R[i];
        r_i=r.str();
        ofstream fout(loc+"H0"+to_string(f)+"_"+r_i+"_"+algo[i]+"_"+base);
        ofstream gout(loc+"H1"+to_string(f)+"_"+r_i+"_"+algo[i]+"_"+base);
        jetform(Algorithm[j], input_particles, R[i], fout, gout);
        fout.close();
        gout.close();
        r_i.clear();
        r.clear();
      }
    }
    cout<<"out";
  return 0;
}