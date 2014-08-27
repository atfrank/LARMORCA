//Aaron T. Frank
//Sean M. Law

    
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

//Code generated using: awk '{print "{\""$1":"$2":"$3"\","$4"},"}' larmorD_both.dat | tr '\n' ' '

#include "CASA.hpp"
#include "Molecule.hpp"
#include "Misc.hpp"

#include <fstream>
#include <stdlib.h>


CASA::CASA (Molecule *mol){
    this->initializeRefSASs();
    this->initializeSASA();
    this->initializePolarity();
    this->initializeRadius();
}

double CASA::getRefSASA(const std::string &key){
    if (this->refSASA.find (key) == this->refSASA.end()){
        return 0.0;
    } else {
        return (this->refSASA.at(key));
    }
}

double CASA::getSASA(const std::string &key){
    if (this->SASA.find (key) == this->SASA.end()){
        return 0.0;
    } else {
        return (this->SASA.at(key));
    }
}

double CASA::getPolarity(const std::string &key){
    if (this->polarity.find (key) == this->polarity.end()){
        return 0.0;
    } else {
        return (this->polarity.at(key));
    }
}

double CASA::getRadius(const std::string &key){
    if (this->radius.find (key) == this->radius.end()){
        return 0.0;
    } else {
        return (this->radius.at(key));
    }
}


void CASA::initializeRefSASs(){
  this->refSASA.insert(std::pair<std::string,double>("GLY",47.155));
  this->refSASA.insert(std::pair<std::string,double>("ALA",61.063));
  this->refSASA.insert(std::pair<std::string,double>("SER",69.073));
  this->refSASA.insert(std::pair<std::string,double>("CYS",78.62));
  this->refSASA.insert(std::pair<std::string,double>("PRO",81.03));
  this->refSASA.insert(std::pair<std::string,double>("THR",82.92));
  this->refSASA.insert(std::pair<std::string,double>("ASP",86.142));
  this->refSASA.insert(std::pair<std::string,double>("VAL",90.264));
  this->refSASA.insert(std::pair<std::string,double>("ASN",90.541));
  this->refSASA.insert(std::pair<std::string,double>("GLU",102.57));
  this->refSASA.insert(std::pair<std::string,double>("ILE",105.307));
  this->refSASA.insert(std::pair<std::string,double>("LEU",105.842));
  this->refSASA.insert(std::pair<std::string,double>("GLN",106.534));
  this->refSASA.insert(std::pair<std::string,double>("HIS",110.864));
  this->refSASA.insert(std::pair<std::string,double>("MET",112.01));
  this->refSASA.insert(std::pair<std::string,double>("LYS",120.502));
  this->refSASA.insert(std::pair<std::string,double>("PHE",124.714));
  this->refSASA.insert(std::pair<std::string,double>("TYR",130.885));
  this->refSASA.insert(std::pair<std::string,double>("ARG",138.595));
  this->refSASA.insert(std::pair<std::string,double>("TRP",148.915));
}

void CASA::initializeSASA(){
  this->SASA.insert(std::pair<std::string,double>("ALA",0.99042));
  this->SASA.insert(std::pair<std::string,double>("ARG",0.71552));
  this->SASA.insert(std::pair<std::string,double>("ASN",0.87896));
  this->SASA.insert(std::pair<std::string,double>("ASP",0.88912));
  this->SASA.insert(std::pair<std::string,double>("CYS",0.95123));
  this->SASA.insert(std::pair<std::string,double>("GLN",0.79978));
  this->SASA.insert(std::pair<std::string,double>("GLU",0.77849));
  this->SASA.insert(std::pair<std::string,double>("GLY",1.05203));
  this->SASA.insert(std::pair<std::string,double>("HIS",0.79982));
  this->SASA.insert(std::pair<std::string,double>("ILE",1.04589));
  this->SASA.insert(std::pair<std::string,double>("LEU",1.02160));
  this->SASA.insert(std::pair<std::string,double>("LYS",0.65241));
  this->SASA.insert(std::pair<std::string,double>("MET",0.97725));
  this->SASA.insert(std::pair<std::string,double>("PHE",1.02315));
  this->SASA.insert(std::pair<std::string,double>("PRO",0.98890));
  this->SASA.insert(std::pair<std::string,double>("SER",0.98230));
  this->SASA.insert(std::pair<std::string,double>("THR",0.89542));
  this->SASA.insert(std::pair<std::string,double>("TRP",1.01999));
  this->SASA.insert(std::pair<std::string,double>("TYR",0.86185));
  this->SASA.insert(std::pair<std::string,double>("VAL",0.98227));
}

void CASA::initializePolarity(){
  this->polarity.insert(std::pair<std::string,double>("ALA",0));
  this->polarity.insert(std::pair<std::string,double>("ARG",1));
  this->polarity.insert(std::pair<std::string,double>("ASN",1));
  this->polarity.insert(std::pair<std::string,double>("ASP",1));
  this->polarity.insert(std::pair<std::string,double>("CYS",1));
  this->polarity.insert(std::pair<std::string,double>("GLN",1));
  this->polarity.insert(std::pair<std::string,double>("GLU",1));
  this->polarity.insert(std::pair<std::string,double>("GLY",0));
  this->polarity.insert(std::pair<std::string,double>("HIS",1));
  this->polarity.insert(std::pair<std::string,double>("ILE",0));
  this->polarity.insert(std::pair<std::string,double>("LEU",0));
  this->polarity.insert(std::pair<std::string,double>("LYS",1));
  this->polarity.insert(std::pair<std::string,double>("MET",0));
  this->polarity.insert(std::pair<std::string,double>("PHE",0));
  this->polarity.insert(std::pair<std::string,double>("PRO",0));
  this->polarity.insert(std::pair<std::string,double>("SER",1));
  this->polarity.insert(std::pair<std::string,double>("THR",1));
  this->polarity.insert(std::pair<std::string,double>("TRP",1));
  this->polarity.insert(std::pair<std::string,double>("TYR",1));
  this->polarity.insert(std::pair<std::string,double>("VAL",0));
}

void CASA::initializeRadius(){
  this->radius.insert(std::pair<std::string,double>("ALA",4.01256));
  this->radius.insert(std::pair<std::string,double>("ARG",4.62059));
  this->radius.insert(std::pair<std::string,double>("ASN",4.28348));
  this->radius.insert(std::pair<std::string,double>("ASP",4.27733));
  this->radius.insert(std::pair<std::string,double>("CYS",3.74826));
  this->radius.insert(std::pair<std::string,double>("GLN",4.35527));
  this->radius.insert(std::pair<std::string,double>("GLU",4.29564));
  this->radius.insert(std::pair<std::string,double>("GLY",3.86518));
  this->radius.insert(std::pair<std::string,double>("HIS",4.32923));
  this->radius.insert(std::pair<std::string,double>("ILE",4.57110));
  this->radius.insert(std::pair<std::string,double>("LEU",4.55731));
  this->radius.insert(std::pair<std::string,double>("LYS",4.16101));
  this->radius.insert(std::pair<std::string,double>("MET",4.49251));
  this->radius.insert(std::pair<std::string,double>("PHE",4.69671));
  this->radius.insert(std::pair<std::string,double>("PRO",4.51078));
  this->radius.insert(std::pair<std::string,double>("SER",4.12788));
  this->radius.insert(std::pair<std::string,double>("THR",4.15639));
  this->radius.insert(std::pair<std::string,double>("TRP",5.05033));
  this->radius.insert(std::pair<std::string,double>("TYR",4.51586));
  this->radius.insert(std::pair<std::string,double>("VAL",4.28142));
}
