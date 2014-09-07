//Sean M. Law
//Aaron T. Frank
    
/*
This file is part of MoleTools.

MoleTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MoleTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MoleTools.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Analyze.hpp"
#include "Molecule.hpp"
#include "Chain.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Coor.hpp"
#include "Constants.hpp"
#include "DTree.hpp"
#include "PCASSO.hpp"
#include "LARMORCAP.hpp"
#include "LARMORCA.hpp"
#include "Misc.hpp"

#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <stdlib.h>

Analyze::Analyze (){
  sel.clear();
  mol.clear();
  tdata.clear();
  fdata.clear();
  ndata=0;
  resel=false;
  ifile.clear();
  ofile.clear();
  verbose=false;
}

Analyze::~Analyze (){
  //Do nothing
}

AnalyzeDihedral::AnalyzeDihedral(DihedralEnum thetain){
  this->setDihedralType(thetain);
}

AnalyzePcasso::AnalyzePcasso (std::string delim){
  std::vector<std::string> tokens;
  pout=PREDICT;
  t.clear();
  t.resize(PCASSO::getNTree());
  for (unsigned int i=0; i< PCASSO::getNTree(); i++){
    t.at(i)=new DTree;
    Misc::splitStr(Misc::trim(PCASSO::getTree(i)), " \t", tokens, false);
    t.at(i)->genDTree(tokens, delim);
  }
}

AnalyzeLarmorca::AnalyzeLarmorca (std::string delim){
  unsigned int ntree;
  /* Tree - Nucleus
  t1 - H
  t2 - HA
  t3 - CA
  t4 - C
  t5 - CB
  t6 - N 
  */
  std::vector<std::string> tokens;

  t1.clear();
  t2.clear();
  t3.clear();
  t4.clear();
  t5.clear();
  t6.clear();
  
  ntree=LARMORCAP::getNTree();
  t1.resize(ntree);
  t2.resize(ntree);
  t3.resize(ntree);
  t4.resize(ntree);
  t5.resize(ntree);
  t6.resize(ntree);
  
  for (unsigned int i=0; i< ntree; i++){
    t1.at(i)=new DTree;
    Misc::splitStr(LARMORCAP::getTree(i, "H"), " \t", tokens, false);
    t1.at(i)->genDTree(tokens, delim);
    //tokens.clear();

    t2.at(i)=new DTree;
    Misc::splitStr(LARMORCAP::getTree(i, "HA"), " \t", tokens, false);
    t2.at(i)->genDTree(tokens, delim);
    //tokens.clear();

    t3.at(i)=new DTree;
    Misc::splitStr(LARMORCAP::getTree(i, "CA"), " \t", tokens, false);
    t3.at(i)->genDTree(tokens, delim);
    //tokens.clear();

    t4.at(i)=new DTree;
    Misc::splitStr(LARMORCAP::getTree(i, "C"), " \t", tokens, false);
    t4.at(i)->genDTree(tokens, delim);
    //tokens.clear();

    t5.at(i)=new DTree;
    Misc::splitStr(LARMORCAP::getTree(i, "CB"), " \t", tokens, false);
    t5.at(i)->genDTree(tokens, delim);
    //tokens.clear();
    
    t6.at(i)=new DTree;
    Misc::splitStr(LARMORCAP::getTree(i, "N"), " \t", tokens, false);
    t6.at(i)->genDTree(tokens, delim);
    //tokens.clear();
  }
}

void Analyze::addSel(const std::string& selin){
  this->sel.push_back(selin);
}

std::string Analyze::getSel(const int& element){
  return this->sel.at(element);
}

unsigned int Analyze::getNSel(){
  return this->sel.size();
}

void Analyze::addMol(Molecule* molin){
  this->mol.push_back(molin);
}

void Analyze::setMol(const int& element, Molecule* molin){
  this->mol.at(element)=molin;
}

void Analyze::clearMol(){
  this->mol.clear();
}

void Analyze::resizeNMol(const int sizein){
  this->mol.resize(sizein);
}

void Analyze::readTopology(Molecule* molin, std::string topin){
  if (topin.length() > 0){
    molin->readTopology(topin);
    molin->setMass();
    molin->setCharge();
  }
}

void Analyze::setupMolSel(Molecule* molin){
  Molecule* tmpmol;

  molin->select(this->getSel(0));
  tmpmol=molin->copy();
  this->addMol(tmpmol);
}

void AnalyzeDistance::setupMolSel(Molecule* molin){
  unsigned int i;
  Molecule* tmpmol;

  for (i=0; i< this->getNSel(); i++){
    molin->select(this->getSel(i));
    tmpmol=molin->copy();
    this->addMol(tmpmol);
  }
}

void AnalyzeAngle::setupMolSel(Molecule* molin){
  unsigned int i;
  Molecule* tmpmol;

  for (i=0; i< this->getNSel(); i++){
    molin->select(this->getSel(i));
    tmpmol=molin->copy();
    this->addMol(tmpmol);
  }
}

void AnalyzeDihedral::setupMolSel(Molecule* molin){
  unsigned int i;
  Molecule *tmpmol, *mol1, *mol2, *mol3, *mol4;
  Residue *res;
  std::string tmpsel;
  std::stringstream resid, lastResId, nextResId;

  if (this->getDihedralType() == GENERAL){
    for (i=0; i< this->getNSel(); i++){
      molin->select(this->getSel(i));
      tmpmol=molin->copy();
      this->addMol(tmpmol);
    }
  }
  else if (this->getDihedralType() == PHI){
    molin->select(this->getSel(0));
    tmpmol=molin->copy();
    for (i=0; i<tmpmol->getResVecSize(); i++){
      res=tmpmol->getResidue(i);
      //Construct selection expression
      resid.str(""); //Clear stringstream
      lastResId.str(""); //Clear stringstream
      resid << res->getResId();
      lastResId << res->getResId()-1;
      tmpsel=res->getChainId()+":"+lastResId.str()+"."+"C_"+res->getChainId()+":"+resid.str()+"."+"CY";
      molin->select(tmpsel, false, false);
      mol1=molin->copy();
      this->addMol(mol1);
      tmpsel=res->getChainId()+":"+resid.str()+"."+"N";
      molin->select(tmpsel, false, false);
      mol2=molin->copy();
      this->addMol(mol2);
      tmpsel=res->getChainId()+":"+resid.str()+"."+"CA";
      molin->select(tmpsel, false, false);
      mol3=molin->copy();
      this->addMol(mol3);
      tmpsel=res->getChainId()+":"+resid.str()+"."+"C";
      molin->select(tmpsel, false, false);
      mol4=molin->copy(); 
      this->addMol(mol4);
    }
  }
  else if (this->getDihedralType() == PSI){
    molin->select(this->getSel(0));
    tmpmol=molin->copy();
    for (i=0; i<tmpmol->getResVecSize(); i++){
      res=tmpmol->getResidue(i);
      //Construct selection expression
      resid.str(""); //Clear stringstream
      nextResId.str(""); //Clear stringstream
      resid << res->getResId();
      nextResId << res->getResId()+1;
      tmpsel=res->getChainId()+":"+resid.str()+"."+"N";
      molin->select(tmpsel, false, false);
      mol1=molin->copy();
      this->addMol(mol1);
      tmpsel=res->getChainId()+":"+resid.str()+"."+"CA";
      molin->select(tmpsel, false, false);
      mol2=molin->copy();
      this->addMol(mol2);
      tmpsel=res->getChainId()+":"+resid.str()+"."+"C";
      molin->select(tmpsel, false, false);
      mol3=molin->copy();
      this->addMol(mol3);
      tmpsel=res->getChainId()+":"+nextResId.str()+"."+"N_"+res->getChainId()+":"+resid.str()+"."+"NT";
      molin->select(tmpsel, false, false);
      mol4=molin->copy();
      this->addMol(mol4);
    }
  }
  else if (this->getDihedralType() == OMEGA){
    molin->select(this->getSel(0));
    tmpmol=molin->copy();
    for (i=0; i<tmpmol->getResVecSize(); i++){
      res=tmpmol->getResidue(i);
      //Construct selection expression
      resid.str(""); //Clear stringstream
      nextResId.str(""); //Clear stringstream
      resid << res->getResId();
      nextResId << res->getResId()+1;
      tmpsel=res->getChainId()+":"+resid.str()+"."+"CA";
      molin->select(tmpsel, false, false);
      mol1=molin->copy();
      this->addMol(mol1);
      tmpsel=res->getChainId()+":"+resid.str()+"."+"C";
      molin->select(tmpsel, false, false);
      mol2=molin->copy();
      this->addMol(mol2);
      tmpsel=res->getChainId()+":"+nextResId.str()+"."+"N";
      molin->select(tmpsel, false, false);
      mol3=molin->copy();
      this->addMol(mol3);
      tmpsel=res->getChainId()+":"+nextResId.str()+"."+"CA";
      molin->select(tmpsel, false, false);
      mol4=molin->copy();
      this->addMol(mol4);
    }
  }
  else if (this->getDihedralType() == CHI1){
    molin->select(this->getSel(0));
    tmpmol=molin->copy();
    for (i=0; i<tmpmol->getResVecSize(); i++){
      res=tmpmol->getResidue(i);
      //Construct selection expression
      resid.str(""); //Clear stringstream
      resid << res->getResId();
      tmpsel=res->getChainId()+":"+resid.str()+"."+"N";
      molin->select(tmpsel, false, false);
      mol1=molin->copy();
      this->addMol(mol1);
      tmpsel=res->getChainId()+":"+resid.str()+"."+"CA";
      molin->select(tmpsel, false, false);
      mol2=molin->copy();
      this->addMol(mol2);
      tmpsel=res->getChainId()+":"+resid.str()+"."+"CB";
      molin->select(tmpsel, false, false);
      mol3=molin->copy();
      this->addMol(mol3);
      tmpsel=res->getChainId()+":"+resid.str()+".";
      if (res->getResName().compare("SER") == 0){
        tmpsel+="OG";
      }
      else if (res->getResName().compare("THR") == 0){
        tmpsel+="OG1";
      }
      else if (res->getResName().compare("VAL") == 0 || res->getResName().compare("ILE") == 0){
        tmpsel+="CG1";
      }
      else if (res->getResName().compare("CYS") == 0){
        tmpsel+="SG";
      }
      else{
        tmpsel+="CG";
      } 
      molin->select(tmpsel, false, false);
      mol4=molin->copy();
      this->addMol(mol4);
    }
  }
  else if (this->getDihedralType() == CHI2){
    molin->select(this->getSel(0));
    tmpmol=molin->copy();
    for (i=0; i<tmpmol->getResVecSize(); i++){
      res=tmpmol->getResidue(i);
      //Construct selection expression
      resid.str(""); //Clear stringstream
      resid << res->getResId();
      tmpsel=res->getChainId()+":"+resid.str()+"."+"CA";
      molin->select(tmpsel, false, false);
      mol1=molin->copy(); 
      this->addMol(mol1);
      tmpsel=res->getChainId()+":"+resid.str()+"."+"CB";
      molin->select(tmpsel, false, false); 
      mol2=molin->copy();
      this->addMol(mol2);
      tmpsel=res->getChainId()+":"+resid.str()+".";
      if (res->getResName().compare("ILE") == 0){
        tmpsel+="CG1";
      }
      else{
        tmpsel+="CG";
      }
      molin->select(tmpsel, false, false);
      mol3=molin->copy();
      this->addMol(mol3);
      tmpsel=res->getChainId()+":"+resid.str()+".";
      if (res->getResName().compare("LEU") == 0 || res->getResName().compare("PHE") == 0 || res->getResName().compare("TRP") == 0 || res->getResName().compare("TYR") == 0){
        tmpsel+="CD1";
      }
      else if (res->getResName().compare("ASN") == 0 || res->getResName().compare("ASP") == 0){
        tmpsel+="OD1";
      }
      else if (res->getResName().compare("HIS") == 0 || res->getResName().compare("HSD") == 0 || res->getResName().compare("HSE") == 0 || res->getResName().compare("HSP") == 0){
        tmpsel+="ND1";
      }
      else if (res->getResName().compare("MET") == 0){
        tmpsel+="SD";
      }
      else{
        tmpsel+="CD";
      }
      molin->select(tmpsel, false, false);
      mol4=molin->copy();
      this->addMol(mol4);
    }
  }
  else{
    std::cerr << "Error: Unrecognized dihedral type" << std::endl;
  }
}

Molecule* Analyze::getMol(const int& element){
  return this->mol.at(element);
}

unsigned int Analyze::getNMol(){
  return this->mol.size();
}

void Analyze::setNData(const int& ndatain){
  ndata=ndatain;
}

int& Analyze::getNData(){
  return ndata;
}

std::vector<double>& Analyze::getTDataVec(){
  return tdata;
}

std::vector<std::vector<double> >& Analyze::getFDataVec(){
  return fdata;
}


void Analyze::setVerbose(bool verbosein){
  verbose=verbosein;
}

bool Analyze::getVerbose(){
  return verbose;
}


//All preAnalysis Functions
void Analyze::preAnalysis(){
  //Do nothing
}

void Analyze::preAnalysis(Molecule* molin, std::string topin){
  this->readTopology(molin, topin);
  this->setupMolSel(molin);
}

void AnalyzeAverage::preAnalysis(Molecule* molin, std::string topin){
  this->readTopology(molin, topin);
  this->setupMolSel(molin);
  Molecule* refmol=molin->clone();
  this->readTopology(refmol, topin);
  this->setupMolSel(refmol);
  //Zero coordinates for future appending
  this->getMol(1)->zeroCoor();
}

void AnalyzePcasso::preAnalysis(Molecule* molin, std::string fin){
  this->setupMolSel(molin);
  //Resize FData in Analyze::pcasso
}

void AnalyzeLarmorca::preAnalysis(Molecule* molin, std::string fin){
  this->setupMolSel(molin);
}


//All runAnalysis Functions

void AnalyzeCOG::runAnalysis(){
  Coor xyz=Analyze::centerOfGeometry(this->getMol(0));
  std::cout << std::fixed;
  std::cout << std::setw(9) << std::right << std::setprecision(3) << xyz.x();
  std::cout << std::setw(9) << std::right << std::setprecision(3) << xyz.y();
  std::cout << std::setw(9) << std::right << std::setprecision(3) << xyz.z();
}

void AnalyzeAverage::runAnalysis(){
  Analyze::averageMol(this->getMol(0), this->getMol(1), getNData());
}

void AnalyzeDistance::runAnalysis(){
  std::cout << std::fixed;
  std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::distance(Analyze::centerOfGeometry(this->getMol(0)), Analyze::centerOfGeometry(this->getMol(1)));
}

void AnalyzeAngle::runAnalysis(){
  std::cout << std::fixed;
  std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::angle(Analyze::centerOfGeometry(this->getMol(0)), Analyze::centerOfGeometry(this->getMol(1)), Analyze::centerOfGeometry(this->getMol(2)));
}

void AnalyzeDihedral::runAnalysis(){
  unsigned int i;

  std::cout << std::fixed;

  for (i=0; i< this->getNMol()-3; i=i+4){
    std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::dihedral(Analyze::centerOfGeometry(this->getMol(i)), Analyze::centerOfGeometry(this->getMol(i+1)), Analyze::centerOfGeometry(this->getMol(i+2)), Analyze::centerOfGeometry(this->getMol(i+3)));
  }
}

void AnalyzePairwiseDistance::runAnalysis(){
//  std::map<std::pair<Atom*, Atom*>, double> pdist;
  std::vector<std::vector<double> > pdist;
  Atom* ai;
  Atom* aj;
  
  this->getMol(0)->assignAtmInx();
  Analyze::pairwiseDistance(this->getMol(0), pdist);

  for (unsigned int i=0; i< this->getMol(0)->getAtmVecSize(); i++){
    ai=this->getMol(0)->getAtom(i);
    for (unsigned int j=i+1; j< this->getMol(0)->getAtmVecSize(); j++){
      aj=this->getMol(0)->getAtom(j);
//      std::cout << "  " << ai->getSummary() << "-" << aj->getSummary();
      std::cout << std::fixed;
      std::cout << std::setw(9) << std::right << std::setprecision(3) << pdist.at(ai->getAtmInx()).at(aj->getAtmInx());
    }
  }
}

void AnalyzePcasso::runAnalysis(){
  std::string line;
  std::vector<std::string> s;
  unsigned int natom;
  bool last;

  Analyze::pcasso(this->getMol(0), this->getFDataVec()); //PCASSO features get stored in the second argument (a 2-D vector double)

  natom=0;
  last=false;

  std::vector<std::vector<double> > &feat=this->getFDataVec();
  std::map<std::string, unsigned int> vote; 
  std::map<std::string, unsigned int>::iterator iter;
  std::string tmpClass;
  std::string maxClass;
  unsigned int maxVote;
  bool majority;
  unsigned int ntree;

  ntree=PCASSO::getNTree();

  if (this->getOutType() == PREDICT || this->getOutType() == PREDICTION){
    for (unsigned int i=0; i< feat.size(); i++){
      vote.clear();
      maxVote=0;
      majority=false;
      for (unsigned int j=0; j< ntree && majority == false; j++){
        tmpClass=t.at(j)->getDTreeClass(feat.at(i));
        if (vote.find(tmpClass) != vote.end()){
          vote.at(tmpClass)++;
          if (vote.at(tmpClass) > maxVote){
            //Find majority vote, method adapted from openCV
            //Is pseudo-random since it depends on the order of the trees
            maxVote=vote.at(tmpClass);
            maxClass=tmpClass;
            if (static_cast<float>(maxVote)/ntree > 0.5){ //Unsigned integer division!
              majority=true;
            }
          }
        }
        else{
          vote.insert(std::pair<std::string, unsigned int>(tmpClass,1));
        }
      }
      if (this->getVerbose() == true && this->getMol(0)->getNAtom() == feat.size()){
        std::cout << " " << this->getMol(0)->getAtom(i)->getSummary();
        std::cout << " " << maxClass ; //Print majority vote
      }
      else{
        std::cout << " " << maxClass; //Print majority vote
      }
    }
  }
  else{
    std::cerr << "Warning: Unrecognized PCASSO output type" << std::endl;
  }

}

void AnalyzeLarmorca::runAnalysis(){
}

void AnalyzeLarmorca::runAnalysisTest(LARMORCA* larm, unsigned int frame, std::string fchemshift, std::string identification, bool analyzeError, bool printError,std::string errorType, bool accuracyAtom){
  Molecule *mol;
  Molecule *mol1;
  Molecule *mol2;
  Molecule *mol3;
  Molecule *mol4;
  Molecule *mol5;
  Molecule *mol6;
  
  std::stringstream fout;
  std::ofstream pdbout1;
  std::ofstream pdbout2;
  std::ofstream pdbout3;
  std::ofstream pdbout4;
  std::ofstream pdbout5;
  std::ofstream pdbout6;
  
  //LARMORCA *larm;
  std::vector<std::vector<double> > &feat=this->getFDataVec();
  double p1;
  double p2;
  double p3;
  double p4;
  double p5;
  double p6;
  std::vector<double> e1;
  std::vector<double> e2;
  std::vector<double> e3;
  std::vector<double> e4;
  std::vector<double> e5;
  std::vector<double> e6;
  
  double expcs;
  double predcs;
  double randcs;
  double errorW;
  double error;
  unsigned int natom;
  std::stringstream resid;
  unsigned int ntree;
  std::string resname;
  std::string key;
  

  mol = NULL;
  //larm = NULL;
  errorW=1.0;
  if (analyzeError == true){
    mol1 = NULL;
    mol2 = NULL;
    mol3 = NULL;
    mol4 = NULL;
    mol5 = NULL;
    mol6 = NULL;
    
    fout << identification << "_error_H.pdb";
    pdbout1.open(fout.str().c_str(), std::ios::out);
    fout.str(std::string());

    fout << identification << "_error_HA.pdb";
    pdbout2.open(fout.str().c_str(), std::ios::out);
    fout.str(std::string());

    fout << identification << "_error_CA.pdb";
    pdbout3.open(fout.str().c_str(), std::ios::out);
    fout.str(std::string());

    fout << identification << "_error_C.pdb";
    pdbout4.open(fout.str().c_str(), std::ios::out);
    fout.str(std::string());

    fout << identification << "_error_CB.pdb";
    pdbout5.open(fout.str().c_str(), std::ios::out);
    fout.str(std::string());

    fout << identification << "_error_N.pdb";
    pdbout6.open(fout.str().c_str(), std::ios::out);
    fout.str(std::string());
  }

  e1.clear();
  e2.clear();
  e3.clear();
  e4.clear();
  e5.clear();
  e6.clear();
  
  

  mol = this->getMol(0);

  /* rename HID, HIE => HIS */
  mol->renameRes("HID","HIS");
  mol->renameRes("HIE","HIS");
  mol->renameRes("HSD","HIS");

  if (analyzeError == true){
    mol1 = mol->clone();
    mol2 = mol->clone();
    mol3 = mol->clone();
    mol4 = mol->clone();
    mol5 = mol->clone();
    mol6 = mol->clone();
  }
  
  Analyze::pcasso(mol, this->getFDataVec()); //PCASSO features get stored in the second argument (a 2-D vector double)
  //larm = new LARMORCA(fchemshift);
  natom=0;
  
  ntree=LARMORCAP::getNTree();
  for (unsigned int i=0; i< mol->getNAtomSelected(); i++){
    p1=p2=p3=p4=p5=p6=0.0;
    resname = mol->getResidue(natom)->getResName();
    resid << mol->getResidue(natom)->getResId();
    feat[i].insert(feat[i].begin(),0.0);
    for (unsigned int j=0; j< ntree; j++){
      feat[i].at(0) = larm->getRandomShift("H:"+resname); // NOTE: adding RandomCoil reference chemical shifts to the front of the feature vectors
      p1 += atof(t1.at(j)->getDTreeClass(feat.at(i)).c_str());
      
      feat[i].at(0) = larm->getRandomShift("HA:"+resname);
      p2 += atof(t2.at(j)->getDTreeClass(feat.at(i)).c_str());
      
      feat[i].at(0) = larm->getRandomShift("CA:"+resname);
      p3 += atof(t3.at(j)->getDTreeClass(feat.at(i)).c_str());
      
      feat[i].at(0) = larm->getRandomShift("C:"+resname);
      p4 += atof(t4.at(j)->getDTreeClass(feat.at(i)).c_str());
      
      feat[i].at(0) = larm->getRandomShift("CB:"+resname);
      p5 += atof(t5.at(j)->getDTreeClass(feat.at(i)).c_str());
      
      feat[i].at(0) = larm->getRandomShift("N:"+resname);
      p6 += atof(t6.at(j)->getDTreeClass(feat.at(i)).c_str());
    }
    p1 /= ntree;
    p2 /= ntree;
    p3 /= ntree;
    p4 /= ntree;
    p5 /= ntree;
    p6 /= ntree;

    key = resid.str()+":H";
    expcs = larm->getExperimentalCS(key);
    randcs = larm->getRandomShift("H:"+resname);
    predcs = p1;
    if((randcs!=999.0) && (fchemshift.length() == 0 || expcs != 0.0)){
      if (printError == false){
        std::cout << frame << " "  << resid.str() << " H " << resname << " " << randcs << " " << predcs  << " " << expcs << " " << identification << std::endl;
      } else if (errorType=="W1MAE" || errorType=="W1RMSE" ){
        if(accuracyAtom==false){
            errorW = 1.0/larm->getMAE("H");
        } else {
            errorW = 1.0/larm->getMAE2("H:"+resname);
        }
        error = errorW*(expcs - predcs);
        e1.push_back(error);
      } else if (errorType=="W2MAE" || errorType=="W2RMSE" ){
        if(accuracyAtom==false){
            errorW = fabs(larm->getR("H"))/larm->getMAE("H");
        } else {
            errorW = fabs(larm->getR2("H:"+resname))/larm->getMAE2("H:"+resname);
        }
        error = errorW*(expcs - predcs);
        e1.push_back(error);
      } else if (errorType=="W3MAE" || errorType=="W3RMSE" ){
        if (fabs(expcs - predcs) <= 3*larm->getR("H")){
            if(accuracyAtom==false){
                errorW = 1.0/larm->getMAE("H");
            } else {
                errorW = 1.0/larm->getMAE2("H:"+resname);
            }
            error = errorW*(expcs - predcs);
            e1.push_back(error);
        }
      } else if (errorType=="W4MAE" || errorType=="W4RMSE" ){
        if (fabs(expcs - predcs) <= 3*larm->getR("H")){
            if(accuracyAtom==false){
                errorW = fabs(larm->getR("H"))/larm->getMAE("H");
            } else {
                errorW = fabs(larm->getR2("H:"+resname))/larm->getMAE2("H:"+resname);
            }
            error = errorW*(expcs - predcs);
            e1.push_back(error);
        }
      } else {
        error = errorW*(expcs - predcs);
        e1.push_back(error);
      }
      if (analyzeError == true){
        mol1->getAtom(i)->setBFac(fabs(expcs - predcs));
      }
    }


    key = resid.str()+":HA";
    expcs = larm->getExperimentalCS(key);
    randcs = larm->getRandomShift("HA:"+resname);
    predcs = p2;
    if((randcs!=999.0) && (fchemshift.length() == 0 || expcs != 0.0)){
      if (printError == false){
        std::cout << frame << " "  << resid.str() << " HA " << resname << " " << randcs << " " << predcs  << " " << expcs << " " << identification << std::endl;
      } else if (errorType=="W1MAE" || errorType=="W1RMSE" ){
        if(accuracyAtom==false){
            errorW = 1.0/larm->getMAE("HA");
        } else {
            errorW = 1.0/larm->getMAE2("HA:"+resname);
        }
        error = errorW*(expcs - predcs);
        e2.push_back(error);
      } else if (errorType=="W2MAE" || errorType=="W2RMSE" ){
        if(accuracyAtom==false){
            errorW = fabs(larm->getR("HA"))/larm->getMAE("HA");
        } else {
            errorW = fabs(larm->getR2("HA:"+resname))/larm->getMAE2("HA:"+resname);
        }
        error = errorW*(expcs - predcs);
        e2.push_back(error);
      } else if (errorType=="W3MAE" || errorType=="W3RMSE" ){
        if (fabs(expcs - predcs) <= 3*larm->getR("HA")){
            if(accuracyAtom==false){
                errorW = 1.0/larm->getMAE("HA");
            } else {
                errorW = 1.0/larm->getMAE2("HA:"+resname);
            }
            error = errorW*(expcs - predcs);
            e2.push_back(error);
        }
      } else if (errorType=="W4MAE" || errorType=="W4RMSE" ){
        if (fabs(expcs - predcs) <= 3*larm->getR("HA")){
            if(accuracyAtom==false){
                errorW = fabs(larm->getR("HA"))/larm->getMAE("HA");
            } else {
                errorW = fabs(larm->getR2("HA:"+resname))/larm->getMAE2("HA:"+resname);
            }
            error = errorW*(expcs - predcs);
            e2.push_back(error);
        }
      } else {
        error = errorW*(expcs - predcs);
        e2.push_back(error);
      }
      if (analyzeError == true){
        mol2->getAtom(i)->setBFac(fabs(expcs - predcs));
      }
    }


    key = resid.str()+":CA";
    expcs = larm->getExperimentalCS(key);
    randcs = larm->getRandomShift("CA:"+resname);
    predcs = p3;
    if((randcs!=999.0) && (fchemshift.length() == 0 || expcs != 0.0)){
      if (printError == false){
        std::cout << frame << " "  << resid.str() << " CA " << resname << " " << randcs << " " << predcs  << " " << expcs << " " << identification << std::endl;
      } else if (errorType=="W1MAE" || errorType=="W1RMSE" ){
        if(accuracyAtom==false){
            errorW = 1.0/larm->getMAE("CA");
        } else {
            errorW = 1.0/larm->getMAE2("CA:"+resname);
        }
        error = errorW*(expcs - predcs);
        e3.push_back(error);
      } else if (errorType=="W2MAE" || errorType=="W2RMSE" ){
        if(accuracyAtom==false){
            errorW = fabs(larm->getR("CA"))/larm->getMAE("CA");
        } else {
            errorW = fabs(larm->getR2("CA:"+resname))/larm->getMAE2("CA:"+resname);
        }
        error = errorW*(expcs - predcs);
        e3.push_back(error);
      } else if (errorType=="W3MAE" || errorType=="W3RMSE" ){
        if (fabs(expcs - predcs) <= 3*larm->getR("CA")){
            if(accuracyAtom==false){
                errorW = 1.0/larm->getMAE("CA");
            } else {
                errorW = 1.0/larm->getMAE2("CA:"+resname);
            }
            error = errorW*(expcs - predcs);
            e3.push_back(error);
        }
      } else if (errorType=="W4MAE" || errorType=="W4RMSE" ){
        if (fabs(expcs - predcs) <= 3*larm->getR("CA")){
            if(accuracyAtom==false){
                errorW = fabs(larm->getR("CA"))/larm->getMAE("CA");
            } else {
                errorW = fabs(larm->getR2("CA:"+resname))/larm->getMAE2("CA:"+resname);
            }
            error = errorW*(expcs - predcs);
            e3.push_back(error);
        }
      } else {
        error = errorW*(expcs - predcs);
        e3.push_back(error);
      }
      if (analyzeError == true){
        mol3->getAtom(i)->setBFac(fabs(expcs - predcs));
      }
    }

    key = resid.str()+":C";
    expcs = larm->getExperimentalCS(key);
    randcs = larm->getRandomShift("C:"+resname);
    predcs = p4;
    if((randcs!=999.0) && (fchemshift.length() == 0 || expcs != 0.0)){
      if (printError == false){
        std::cout << frame << " "  << resid.str() << " C " << resname << " " << randcs << " " << predcs  << " " << expcs << " " << identification << std::endl;
      } else if (errorType=="W1MAE" || errorType=="W1RMSE" ){
        if(accuracyAtom==false){
            errorW = 1.0/larm->getMAE("C");
        } else {
            errorW = 1.0/larm->getMAE2("C:"+resname);
        }
        error = errorW*(expcs - predcs);
        e4.push_back(error);
      } else if (errorType=="W2MAE" || errorType=="W2RMSE" ){
        if(accuracyAtom==false){
            errorW = fabs(larm->getR("C"))/larm->getMAE("C");
        } else {
            errorW = fabs(larm->getR2("C:"+resname))/larm->getMAE2("C:"+resname);
        }
        error = errorW*(expcs - predcs);
        e4.push_back(error);
      } else if (errorType=="W3MAE" || errorType=="W3RMSE" ){
        if (fabs(expcs - predcs) <= 3*larm->getR("C")){
            if(accuracyAtom==false){
                errorW = 1.0/larm->getMAE("C");
            } else {
                errorW = 1.0/larm->getMAE2("C:"+resname);
            }
            error = errorW*(expcs - predcs);
            e4.push_back(error);
        }
      } else if (errorType=="W4MAE" || errorType=="W4RMSE" ){
        if (fabs(expcs - predcs) <= 3*larm->getR("C")){
            if(accuracyAtom==false){
                errorW = fabs(larm->getR("C"))/larm->getMAE("C");
            } else {
                errorW = fabs(larm->getR2("C:"+resname))/larm->getMAE2("C:"+resname);
            }
            error = errorW*(expcs - predcs);
            e4.push_back(error);
        }
      } else {
        error = errorW*(expcs - predcs);
        e4.push_back(error);
      }
      if (analyzeError == true){
        mol4->getAtom(i)->setBFac(fabs(expcs - predcs));
      }
    }

    key = resid.str()+":CB";
    expcs = larm->getExperimentalCS(key);
    randcs = larm->getRandomShift("CB:"+resname);
    predcs = p5;
    if((randcs!=999.0) && (fchemshift.length() == 0 || expcs != 0.0)){
      if (printError == false){
        std::cout << frame << " "  << resid.str() << " CB " << resname << " " << randcs << " " << predcs  << " " << expcs << " " << identification << std::endl;
      } else if (errorType=="W1MAE" || errorType=="W1RMSE" ){
        if(accuracyAtom==false){
            errorW = 1.0/larm->getMAE("CB");
        } else {
            errorW = 1.0/larm->getMAE2("CB:"+resname);
        }
        error = errorW*(expcs - predcs);
        e5.push_back(error);
      } else if (errorType=="W2MAE" || errorType=="W2RMSE" ){
        if(accuracyAtom==false){
            errorW = fabs(larm->getR("CB"))/larm->getMAE("CB");
        } else {
            errorW = fabs(larm->getR2("CB:"+resname))/larm->getMAE2("CB:"+resname);
        }
        error = errorW*(expcs - predcs);
        e5.push_back(error);
      } else if (errorType=="W3MAE" || errorType=="W3RMSE" ){
        if (fabs(expcs - predcs) <= 3*larm->getR("CB")){
            if(accuracyAtom==false){
                errorW = 1.0/larm->getMAE("CB");
            } else {
                errorW = 1.0/larm->getMAE2("CB:"+resname);
            }
            error = errorW*(expcs - predcs);
            e5.push_back(error);
        }
      } else if (errorType=="W4MAE" || errorType=="W4RMSE" ){
        if (fabs(expcs - predcs) <= 3*larm->getR("CB")){
            if(accuracyAtom==false){
                errorW = fabs(larm->getR("CB"))/larm->getMAE("CB");
            } else {
                errorW = fabs(larm->getR2("CB:"+resname))/larm->getMAE2("CB:"+resname);
            }
            error = errorW*(expcs - predcs);
            e5.push_back(error);
        }
      } else {
        error = errorW*(expcs - predcs);
        e5.push_back(error);
      }
      if (analyzeError == true){
        mol5->getAtom(i)->setBFac(fabs(expcs - predcs));
      }
    }

    key = resid.str()+":N";
    expcs = larm->getExperimentalCS(key);
    randcs = larm->getRandomShift("N:"+resname);
    predcs = p6;
    if((randcs!=999.0) && (fchemshift.length() == 0 || expcs != 0.0)){
      if (printError == false){
        std::cout << frame << " "  << resid.str() << " N " << resname << " " << randcs << " " << predcs  << " " << expcs << " " << identification << std::endl;
      } else if (errorType=="W1MAE" || errorType=="W1RMSE" ){
        if(accuracyAtom==false){
            errorW = 1.0/larm->getMAE("N");
        } else {
            errorW = 1.0/larm->getMAE2("N:"+resname);
        }
        error = errorW*(expcs - predcs);
        e6.push_back(error);
      } else if (errorType=="W2MAE" || errorType=="W2RMSE" ){
        if(accuracyAtom==false){
            errorW = fabs(larm->getR("N"))/larm->getMAE("N");
        } else {
            errorW = fabs(larm->getR2("N:"+resname))/larm->getMAE2("N:"+resname);
        }
        error = errorW*(expcs - predcs);
        e6.push_back(error);
      } else if (errorType=="W3MAE" || errorType=="W3RMSE" ){
        if (fabs(expcs - predcs) <= 3*larm->getR("N")){
            if(accuracyAtom==false){
                errorW = 1.0/larm->getMAE("N");
            } else {
                errorW = 1.0/larm->getMAE2("N:"+resname);
            }
            error = errorW*(expcs - predcs);
            e6.push_back(error);
        }
      } else if (errorType=="W4MAE" || errorType=="W4RMSE" ){
        if (fabs(expcs - predcs) <= 3*larm->getR("N")){
            if(accuracyAtom==false){
                errorW = fabs(larm->getR("N"))/larm->getMAE("N");
            } else {
                errorW = fabs(larm->getR2("N:"+resname))/larm->getMAE2("N:"+resname);
            }
            error = errorW*(expcs - predcs);
            e6.push_back(error);
        }
      } else {
        error = errorW*(expcs - predcs);
        e6.push_back(error);
      }
      if (analyzeError == true){
        mol6->getAtom(i)->setBFac(fabs(expcs - predcs));
      }
    }
    natom ++;
    resid.str("");
  }
  if (analyzeError == true){
    pdbout1 << mol1->writePDB(true,false);
    pdbout2 << mol2->writePDB(true,false);
    pdbout3 << mol3->writePDB(true,false);
    pdbout4 << mol4->writePDB(true,false);
    pdbout5 << mol5->writePDB(true,false);
    pdbout6 << mol6->writePDB(true,false);

    pdbout1.close();
    pdbout2.close();
    pdbout3.close();
    pdbout4.close();
    pdbout5.close();
    pdbout6.close();
  }
  if (printError == true){
    if (errorType.find("RMSE")==std::string::npos){
      std::cout << frame << " " <<  Misc::mae(e1) << " " << Misc::mae(e2) << " " << Misc::mae(e3)  << " " << Misc::mae(e4)  << " " << Misc::mae(e5)  << " " << Misc::mae(e6) << " " << Misc::mae(e1)+Misc::mae(e2)+Misc::mae(e3)+Misc::mae(e4)+Misc::mae(e5)+Misc::mae(e6) << " " << identification << std::endl;
    } else {
      std::cout << frame << " " <<  Misc::rmse(e1) << " " << Misc::rmse(e2) << " " << Misc::rmse(e3)  << " " << Misc::rmse(e4)  << " " << Misc::rmse(e5)  << " " << Misc::rmse(e6) << " " << Misc::rmse(e1)+Misc::rmse(e2)+Misc::rmse(e3)+Misc::rmse(e4)+Misc::rmse(e5)+Misc::rmse(e6) << " " << identification << std::endl;
    }
  }
}

//All postAnalysis functions

void Analyze::postAnalysis(){
  //Do nothing
}

void AnalyzeAverage::postAnalysis(){
  unsigned int i;
  Atom *atm;

  std::cout << std::fixed;
  for(i=0; i< this->getMol(1)->getNAtomSelected(); i++){
    atm=this->getMol(1)->getAtom(i);
    if (atm->getSel() == false){
      continue;
    }
    atm->setCoor(atm->getCoor()/getNData());
  }
  this->getMol(1)->writePDB();
}

void AnalyzePcasso::postAnalysis(){
  while (!t.empty()){
    t.back()->delDTree();
    t.pop_back();
  }
}

//Basic analysis functions

Coor Analyze::centerOfGeometry(Molecule *mol, bool selFlag){
  Coor cog=Coor(0.0, 0.0, 0.0);

  for (unsigned int i=0; i< mol->getAtmVecSize(); i++){
    if (selFlag == true && mol->getAtom(i)->getSel() == false){
      continue;
    }
    cog+=mol->getAtom(i)->getCoor();
  }

  if (selFlag == true){
    cog/=mol->getNAtomSelected();
  }
  else{
    cog/=mol->getNAtom();
  }

  return cog;
}


void Analyze::averageMol (Molecule* cmpmol, Molecule* refmol, int &ndataIO){
  unsigned int i,j;
  Atom *atm;
  std::vector<Coor> coor;

  //Check selection sizes and resize matrices
  if (cmpmol->getNAtomSelected() != refmol->getNAtomSelected()){
    std::cerr << std::endl << "Error: Atom number mismatch in RMSD calculation" << std::endl;
    std::cerr << "CMP-NATOM: " << cmpmol->getNAtomSelected() << ", ";
    std::cerr << "REF-NATOM: " << refmol->getNAtomSelected() << std::endl;
  }
  else{
    for (i=0; i< cmpmol->getNAtom(); i++){
      atm=cmpmol->getAtom(i);
      if (atm->getSel() == false){
        continue;
      }
      coor.push_back(atm->getCoor());
    }

    j=0;
    for (i=0; i< refmol->getNAtom(); i++){
      atm=refmol->getAtom(i);
      if (atm->getSel() == false){
        continue;
      }
      atm->setCoor(coor.at(j)+atm->getCoor());
      j++;
    }
    ndataIO++;
  }
}

double Analyze::distance (const Coor& u, const Coor& v){
  Coor d=u-v;
  return d.norm();
}

double Analyze::angle (const Coor& u, const Coor& v, const Coor& w){
  double angle;
  Coor dx, dy;
  double dp, nx, ny;

  dx=u-v;
  dy=w-v;

  nx=dx.norm();
  ny=dy.norm();

  dp=dx.dot(dy);

  angle=acos(dp/(nx*ny));

  return angle/PI*180.0;
}

double Analyze::dihedral (const Coor& t, const Coor& u, const Coor& v, const Coor& w) {
  double dihedral;
  Coor dx, dy, dz, p1, p2, p3;
  double np1, np2, dp1, dp2, ts;

  dx=t-u;
  dy=u-v;
  dz=w-v; //This is correct!

  p1=dx.cross(dy);

  np1=p1.norm();
  p1=p1/np1;

  p2=dz.cross(dy);
  np2=p2.norm();
  p2=p2/np2;

  dp1=p1.dot(p2); //Dot product

  ts=1.0-dp1*dp1;
  ts=(ts<0.0)?0.0:sqrt(ts);
  dihedral=PI/2.0-atan2(dp1,ts);

  p3=p1.cross(p2);

  dp2=p3.dot(dy); //Dot product

  if (dp2 > 0.0){
    dihedral=-dihedral;
  }

  return dihedral/PI*180.0;
}

double Analyze::distance (Molecule* sel1, Molecule* sel2, bool selFlag){
  return Analyze::distance(Analyze::centerOfGeometry(sel1,selFlag), Analyze::centerOfGeometry(sel2,selFlag));
}

double Analyze::angle (Molecule* sel1, Molecule* sel2, Molecule* sel3, bool selFlag){
  return Analyze::angle(Analyze::centerOfGeometry(sel1,selFlag), Analyze::centerOfGeometry(sel2,selFlag), Analyze::centerOfGeometry(sel3,selFlag));
}

double Analyze::dihedral (Molecule* sel1, Molecule* sel2, Molecule* sel3, Molecule* sel4, bool selFlag){
  return Analyze::dihedral(Analyze::centerOfGeometry(sel1,selFlag), Analyze::centerOfGeometry(sel2,selFlag), Analyze::centerOfGeometry(sel3,selFlag), Analyze::centerOfGeometry(sel4,selFlag));
}

void Analyze::pairwiseDistance(Molecule *mol, std::vector<std::vector<double> >& pdin){
  std::vector<Atom*>::iterator ai;
  std::vector<Atom*>::iterator aj;
  unsigned int natom;
  unsigned int aiInx;
  bool flag;

  natom=mol->getAtmVecSize();

  pdin.clear();
  pdin.resize(natom);

  for (ai=mol->getAtmVec().begin(); ai != mol->getAtmVec().end(); ++ai){
    aiInx=(*ai)->getAtmInx();
    pdin.at(aiInx).resize(natom);
    pdin.at(aiInx).at(aiInx)=0.0; //Zero diagonal
  }

  for (ai=mol->getAtmVec().begin(); ai != mol->getAtmVec().end(); ++ai){
    aiInx=(*ai)->getAtmInx();
    if ((*ai)->getX() < 9999.9){
      flag=true;
    }
    else{
      flag=false;
    }
    //Lower Triangle
    for (aj=mol->getAtmVec().begin(); aj != ai; ++aj){
      if (flag && (*aj)->getX() < 9999.9){
        pdin.at(aiInx).at((*aj)->getAtmInx())=Analyze::distance((*ai)->getCoor(), (*aj)->getCoor());
      }
      else{
        pdin.at(aiInx).at((*aj)->getAtmInx())=9999.9;
      }
    }
  
    //Upper Triangle
    for (aj=ai+1; aj != mol->getAtmVec().end(); ++aj){
      if (flag && (*aj)->getX() < 9999.9){
        pdin.at(aiInx).at((*aj)->getAtmInx())=Analyze::distance((*ai)->getCoor(), (*aj)->getCoor());
      }
      else{
        pdin.at(aiInx).at((*aj)->getAtmInx())=9999.9;
      }
    }
  }

}

void Analyze::allAnglesDihedrals(Molecule *mol, std::vector<std::vector<double> >& anglesin){
  unsigned int i, j;
  Chain *c;
  Atom *atmI, *iMinusTwo, *iMinusOne, *iPlusOne, *iPlusTwo, *iPlusThree;
  unsigned int size;
  std::vector<double> angles;

  anglesin.resize(mol->getAtmVecSize());

  for (i=0; i< mol->getChnVecSize(); i++){
    c=mol->getChain(i);
    size=c->getAtmVecSize();
    for (j=0; j< c->getAtmVecSize(); j++){
      atmI=c->getAtom(j);
      iMinusTwo=NULL;
      iMinusOne=NULL;
      iPlusOne=NULL;
      iPlusTwo=NULL;
      iPlusThree=NULL;

      if (j > 1 && atmI->getResId()-2 == c->getAtom(j-2)->getResId()){
        iMinusTwo=c->getAtom(j-2);
      }
      else{
        if (j > 1 && atmI->getResId() == c->getAtom(j-2)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j-1)->getICode(),0,1) != 0)){
          iMinusTwo=c->getAtom(j-2);
        }
      }

      if (j > 0 && atmI->getResId()-1 == c->getAtom(j-1)->getResId()){
        iMinusOne=c->getAtom(j-1);
      }
      else{
        if (j > 0 && atmI->getResId() == c->getAtom(j-1)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j-1)->getICode(),0,1) != 0)){
          iMinusOne=c->getAtom(j-1);
        }
      }

      if (j+1 < size && atmI->getResId()+1 == c->getAtom(j+1)->getResId()){
        iPlusOne=c->getAtom(j+1);
      }
      else{
        if (j+1 < size && atmI->getResId() == c->getAtom(j+1)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j+1)->getICode(),0,1) != 0)){
          iPlusOne=c->getAtom(j+1);
        }
      }

      if (j+2 < size && atmI->getResId()+2 == c->getAtom(j+2)->getResId()){
        iPlusTwo=c->getAtom(j+2);
      }
      else{
        if (j+2 < size && atmI->getResId() == c->getAtom(j+2)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j+2)->getICode(),0,1) != 0)){
          iPlusTwo=c->getAtom(j+2);
        }
      }

      if (j+3 < size && atmI->getResId()+3 == c->getAtom(j+3)->getResId()){
        iPlusThree=c->getAtom(j+3);
      }
      else{
        if (j+3 < size && atmI->getResId() == c->getAtom(j+3)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j+3)->getICode(),0,1) !=0)){
          iPlusThree=c->getAtom(j+3);
        }
      }
      
      angles.clear();
      angles.resize(3, 9999.9);
      if (iMinusOne != NULL && iPlusOne != NULL){
        //Get Angle
        angles.at(0)=Analyze::angle(iMinusOne->getCoor(), atmI->getCoor(), iPlusOne->getCoor());
      }
      if (iPlusOne != NULL && iPlusTwo != NULL && iPlusThree != NULL){
        //Get Dihedral
        angles.at(1)=Analyze::dihedral(atmI->getCoor(), iPlusOne->getCoor(), iPlusTwo->getCoor(), iPlusThree->getCoor());
      }
      if (iMinusTwo != NULL && iPlusTwo != NULL){
        //Get Wide Angle (i-2, i, i+2)
        angles.at(2)=Analyze::angle(iMinusTwo->getCoor(), atmI->getCoor(), iPlusTwo->getCoor());
      }
      anglesin.at(atmI->getAtmInx()).resize(3);
      anglesin.at(atmI->getAtmInx())=angles;
    }
  }
}

void Analyze::pcasso(Molecule* mol, std::vector<std::vector<double> > &fdataIO){
  Chain *c;
  Atom *ai, *aj, *ak;
  unsigned int i, j, start;
  double defVal;
  double iMinus6;//For non-local contacts
  double iPlus6; //For non-local contacts
  unsigned minx;
  unsigned pinx;
  int diffResId;
  std::vector<std::vector<double> > caPairDist; //Ca-Ca Distances
  std::vector<std::vector<double> > pcPairDist; //Pc-Pc Distances
  std::vector<std::vector<double> > caAngles; //Ca-Ca Angle/Diehdral/Wide
  std::vector<std::vector<double> > pcAngles; //Pc-Pc Angle/Dihedral/Wide
  Molecule *camol, *pcmol;
  unsigned int natom;
  double dist;

  defVal=9999.9;
  camol=NULL;
  pcmol=NULL;

  mol->storeSel();
  mol->select(":.CA");
  camol=mol->clone(true,true); //Copy selection, keep original
  pcmol=mol->clone(true,true);
  mol->recallSel(); //Restore original selection
  mol->eraseSel();

  camol->assignAtmInx();
  pcmol->assignAtmInx();

  //Analyze all C-alpha first
  Analyze::pairwiseDistance(camol, caPairDist);
  Analyze::allAnglesDihedrals(camol, caAngles);

  natom=camol->getAtmVecSize();

  for (unsigned int ichain=0; ichain < camol->getChnVecSize(); ichain++){
    c=camol->getChain(ichain);
    for (unsigned int iatom=0; iatom < c->getAtmVecSize(); iatom++){
      ai=c->getAtom(iatom);
      ai->clearData();
      i=ai->getAtmInx();
      //i-5, i-4, i-3, i-2, i-1, i+1, i+2, i+3, i+4, i+5 Distances
      //Deal with unsigned int subtraction from zero
      if (iatom == 0){
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-0;
      }
      else if (iatom == 1){
        //std::cout << defVal << " ";
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-1;
      }
      else if (iatom == 2){
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-2;
      }
      else if (iatom == 3){
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-3;
      }
      else if (iatom == 4){
        ai->addData(defVal);  
        start=iatom-4;
      }
      else{
        start=iatom-5;
      }
      for (j=start; j<= iatom+5; j++){
        aj=c->getAtom(j);
        if (aj == NULL){
          ai->addData(defVal);
        }
        else if (j == iatom){
          //Distance == 0
          continue;
        }
        else{
          ai->addData(caPairDist.at(i).at(aj->getAtmInx()));
        }
      }

      //Angles and Dihedrals
      for (j=0; j< caAngles.at(i).size(); j++){
        ai->addData(caAngles.at(i).at(j));
      }
    
      //Shortest non-local contact distance, >= i+6 and <= i-6
      iPlus6=1E10;
      iMinus6=1E10;
      pinx=natom;
      minx=natom;

      for (j=0; j< natom; j++){
        aj=camol->getAtom(j);

        if (ai == aj){
          continue;
        }
        
        dist=caPairDist.at(i).at(j);

        //i+6
        //Assess distance first to avoid unnecessary string comparison
        if (dist < iPlus6){
          if (ai->getChainId().compare(aj->getChainId()) != 0){
            //atom i and atom j are on different chains
            pinx=j;
            iPlus6=dist;
          }
          else{
            //atom i and atom j are on the same chain
            diffResId=aj->getResId() - ai->getResId();
            if (diffResId >= 6){
              pinx=j;
              iPlus6=dist;
            }
            else{
              if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ai->getICode()) >= 6)){
                pinx=j;
                iPlus6=dist;
              }
            }
          }
        }

        //i-6
        //Assess distance first to avoid unnecessary string comparison
        if (dist < iMinus6){
          if (ai->getChainId().compare(aj->getChainId()) != 0){
          //atom i and atom j are on different chains
            minx=j;
            iMinus6=dist;
          }
          else{
            //atom i and atom j are on the same chain
            diffResId=aj->getResId() - ai->getResId();
            if (diffResId <= -6){
              minx=j;
              iMinus6=dist;
            }
            else{
              if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ai->getICode()) <= -6)){
                minx=j;
                iMinus6=dist;
              }
            }
          }
        }
      }

      int k;
      unsigned int q;
      int max;
      max=10;
      if (iatom == 0){
        for (k=0; k< max; k++){
          ai->addData(defVal);
        }
        start=0;
      }
      else{
        start=iatom-1;
      }
      for (q=start; q<= iatom+1; q++){
        if (q < c->getAtmVecSize()){
          ak=c->getAtom(q);
          for (k=static_cast<int>(pinx)-2; k<=static_cast<int>(pinx)+2; k++){
            if (k >= 0 && k < static_cast<int>(natom)){
              aj=camol->getAtom(k);
              ai->addData(caPairDist.at(ak->getAtmInx()).at(aj->getAtmInx()));
            }
            else{
              ai->addData(defVal);
            }
          }
          for (k=static_cast<int>(minx)-2; k<=static_cast<int>(minx)+2; k++){
            if (k >= 0 && k < static_cast<int>(natom)){
              aj=camol->getAtom(k);
              ai->addData(caPairDist.at(ak->getAtmInx()).at(aj->getAtmInx()));
            }
            else{
              ai->addData(defVal);
            } 
          }
        }
        else{
          for (k=0; k< max; k++){
            ai->addData(defVal);
          }
        }
      }

    } //Loop through atoms
  }//Loop through chains

  //Analyze all pseudocenter
  pcmol->modPseudoCenter();
  Analyze::pairwiseDistance(pcmol, pcPairDist);
  Analyze::allAnglesDihedrals(pcmol, pcAngles);

  natom=pcmol->getAtmVecSize();

  for (unsigned int ichain=0; ichain < pcmol->getChnVecSize(); ichain++){
    c=pcmol->getChain(ichain);
    for (unsigned int iatom=0; iatom < c->getAtmVecSize(); iatom++){
      ai=camol->getChain(ichain)->getAtom(iatom); //From C-alpha
      ak=c->getAtom(iatom);
      i=ak->getAtmInx();
      //i-5, i-4, i-3, i-2, i-1, i+1, i+2, i+3, i+4, i+5 Distances
      //Deal with unsigned int subtraction from zero
      if (iatom == 0){
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-0;
      }
      else if (iatom == 1){
        //std::cout << defVal << " ";
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-1;
      }
      else if (iatom == 2){
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-2;
      }
      else if (iatom == 3){
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-3;
      }
      else if (iatom == 4){
        ai->addData(defVal);  
        start=iatom-4;
      }
      else{
        start=iatom-5;
      }
      for (j=start; j<= iatom+5; j++){
        aj=c->getAtom(j);
        if (aj == NULL){
          ai->addData(defVal);
        }
        else if (j == iatom){
          //Distance == 0
          continue;
        }
        else{
          ai->addData(pcPairDist.at(i).at(aj->getAtmInx()));
        }
      }

      //Angles and Dihedrals
      for (j=0; j< pcAngles.at(ak->getAtmInx()).size(); j++){
        ai->addData(pcAngles.at(ak->getAtmInx()).at(j));
      }
    
      //Shortest non-local contact distance, >= i+6 and <= i-6
      iPlus6=1E10;
      iMinus6=1E10;
      pinx=natom;
      minx=natom;
      for (j=0; j< natom; j++){
        aj=pcmol->getAtom(j);
        if (ak == aj){
          continue;
        }

        dist=pcPairDist.at(i).at(j);

        //i+6
        //Assess distance first to avoid unnecessary string comparison
        if (dist < iPlus6){
          if (ak->getChainId().compare(aj->getChainId()) != 0){
            //atom i and atom j are on different chains
            pinx=j;
            iPlus6=dist;
          }
          else{
            //atom i and atom j are on the same chain
            diffResId=aj->getResId() - ak->getResId();
            if (diffResId >= 6){
              pinx=j;
              iPlus6=dist;
            }
            else{
              if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ak->getICode()) >= 6)){
                pinx=j;
                iPlus6=dist;
              }
            }
          }
        }

        //i-6
        //Assess distance first to avoid unnecessary string comparison
        if (dist < iMinus6){
          if (ak->getChainId().compare(aj->getChainId()) != 0){
            //atom i and atom j are on different chains
            minx=j;
            iMinus6=dist;
          }
          else{
            //atom i and atom j are on the same chain
            diffResId=aj->getResId() - ak->getResId();
            if (diffResId <= -6){
              minx=j;
              iMinus6=dist;
            }
            else{
              if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ak->getICode()) <= -6)){
                minx=j;
                iMinus6=dist;
              }
            }
          }
        }
      }

      int k;
      unsigned int q;
      int max;
      max=10;
      if (iatom == 0){
        for (k=0; k< max; k++){
          ai->addData(defVal);
        }
        start=0;
      }
      else{
        start=iatom-1;
      }
      for (q=start; q<= iatom+1; q++){
        if (q < c->getAtmVecSize()){
          ak=c->getAtom(q);
          for (k=static_cast<int>(pinx)-2; k<=static_cast<int>(pinx)+2; k++){
            if (k >= 0 && k < static_cast<int>(natom)){
              aj=pcmol->getAtom(k);
              ai->addData(pcPairDist.at(ak->getAtmInx()).at(aj->getAtmInx()));
            }
            else{
              ai->addData(defVal);
            }
          }
          for (k=static_cast<int>(minx)-2; k<=static_cast<int>(minx)+2; k++){
            if (k >= 0 && k < static_cast<int>(natom)){
              aj=pcmol->getAtom(k);
              ai->addData(pcPairDist.at(ak->getAtmInx()).at(aj->getAtmInx()));
            }
            else{
              ai->addData(defVal);
            } 
          }
        }
        else{
          for (k=0; k< max; k++){
            ai->addData(defVal);
          }
        }
      }

    } //Loop through atoms
  }//Loop through chains

  //Always clear and resize!
  fdataIO.clear();
  fdataIO.resize(camol->getNAtom());

  //Store features
  natom=0; //This is needed since natom is used for other things above
  for (unsigned int ichain=0; ichain < camol->getChnVecSize(); ichain++){
    c=camol->getChain(ichain);
    for (unsigned int iatom=0; iatom < c->getAtmVecSize(); iatom++){
      ai=c->getAtom(iatom);
      fdataIO.at(natom).reserve(3*ai->getDataSize());
  
      //Store S(i)
      for (j=0; j< ai->getDataSize(); j++){
        fdataIO.at(natom).push_back(ai->getDataPoint(j));
      }

      //Store S(i-1)
      if (iatom > 0){
        //Not first atom of chain
        aj=c->getAtom(iatom-1);
        for (j=0; j< ai->getDataSize(); j++){
          if (aj != NULL){
            fdataIO.at(natom).push_back(aj->getDataPoint(j));
          }
          else{
            fdataIO.at(natom).push_back(defVal);
          }
        }
      }
      else{
        //Store S(i-1) which is all default values
        for (j=0; j< ai->getDataSize(); j++){
          fdataIO.at(natom).push_back(defVal);
        }
      }

      //Store S(i+1)
      aj=c->getAtom(iatom+1);
      for (j=0; j< ai->getDataSize(); j++){
        if (aj != NULL){
          fdataIO.at(natom).push_back(aj->getDataPoint(j));
        }
        else{
          fdataIO.at(natom).push_back(defVal);
        }
      }
      natom++;
    } //Loop through atoms
  } //Loop through chains


  if (camol != NULL){
    delete camol;
  }
  if (pcmol != NULL){
    delete pcmol;
  }
}

void AnalyzePcasso::setOutType(PcassoOutEnum pin){
  pout=pin;
}

PcassoOutEnum AnalyzePcasso::getOutType(){
  return pout;
}

void AnalyzeDihedral::setDihedralType(DihedralEnum thetain){
  theta=thetain;
}

DihedralEnum AnalyzeDihedral::getDihedralType(){
  return theta;
}
