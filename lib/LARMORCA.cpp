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

#include "LARMORD.hpp"
#include "LARMORCA.hpp"
#include "Molecule.hpp"
#include "Misc.hpp"

#include <fstream>
#include <stdlib.h>


LARMORCA::LARMORCA (Molecule *mol, const std::string fchemshift){
    /* set random coil chemical shifts */
    this->initializeRandomShifts();

    /* load chemical shifts from file */
    if (fchemshift.length() > 0){
        this->loadCSFile(fchemshift);
    }
    
    /* initialize accuracy measure */
    this->initializeExpectedAccuracy();
}

void LARMORCA::loadCSFile(const std::string fchemshift){
    std::ifstream csFile;
    std::istream* csinp;
    std::string line;
    std::vector<std::string> s;
    std::string n;
    if (fchemshift.length() > 0){
        csFile.open(fchemshift.c_str(), std::ios::in);
        csinp=&csFile;
        while (csinp->good() && !(csinp->eof())){
            getline(*csinp, line);
            Misc::splitStr(line, " ", s, true);
            if (s.size() >= 4){
                n = Misc::trim(s.at(2));
                if( n=="H" || n=="HA" || n=="C" ||n=="CA" || n=="CB" || n=="N" ){
                    this->experimentalCS.insert(std::pair<std::string,double>(Misc::trim(s.at(1))+":"+Misc::trim(s.at(2)),atof(Misc::trim(s.at(3)).c_str())));
                }
            }
        }
    }
}


int LARMORCA::getNShiftAtoms(){
    return (this->shiftAtoms.size());
}

double LARMORCA::getRandomShift(const std::string &key){
    if (this->randomShifts.find (key) == this->randomShifts.end()){
        return 0.0;
    } else {
        return (this->randomShifts.at(key));
    }
}

double LARMORCA::getR(const std::string &key){
    if (this->accuracyR.find (key) == this->accuracyR.end()){
        return 0.0;
    } else {
        return (this->accuracyR.at(key));
    }
}

double LARMORCA::getMAE(const std::string &key){
    if (this->accuracyMAE.find (key) == this->accuracyMAE.end()){
        return 0.0;
    } else {
        return (this->accuracyMAE.at(key));
    }
}


double LARMORCA::getExperimentalCS(const std::string &key){
    if (this->experimentalCS.find (key) != this->experimentalCS.end()){
        return this->experimentalCS.at(key);
    } else {
        return 0.0;
    }
}

void LARMORCA::initializeExpectedAccuracy(){
    this->accuracyMAE.insert(std::pair<std::string,double>("C:ALA",0.938));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:ARG",0.902));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:ASN",0.932));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:ASP",0.883));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:CYS",0.957));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:GLN",0.844));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:GLU",0.866));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:GLY",0.999));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:HIS",1.458));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:ILE",0.970));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:LEU",1.037));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:LYS",0.835));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:MET",2.665));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:PHE",1.092));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:PRO",1.016));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:SER",0.952));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:THR",0.944));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:TRP",0.940));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:TYR",0.960));
    this->accuracyMAE.insert(std::pair<std::string,double>("C:VAL",0.911));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:ALA",1.010));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:ARG",0.996));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:ASN",0.967));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:ASP",0.999));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:CYS",2.195));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:GLN",0.913));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:GLU",1.069));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:GLY",1.145));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:HIS",1.369));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:ILE",1.298));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:LEU",0.940));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:LYS",0.913));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:MET",1.265));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:PHE",1.342));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:PRO",3.610));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:SER",1.097));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:THR",1.378));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:TRP",1.295));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:TYR",1.205));
    this->accuracyMAE.insert(std::pair<std::string,double>("CA:VAL",1.413));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:ALA",1.800));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:ARG",1.164));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:ASN",1.448));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:ASP",1.179));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:CYS",6.932));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:GLN",1.132));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:GLU",1.244));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:HIS",2.305));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:ILE",1.851));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:LEU",1.491));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:LYS",1.161));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:MET",1.832));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:PHE",1.213));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:PRO",1.991));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:SER",2.343));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:THR",2.093));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:TRP",0.947));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:TYR",1.575));
    this->accuracyMAE.insert(std::pair<std::string,double>("CB:VAL",1.318));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:ALA",0.390));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:ARG",0.402));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:ASN",0.396));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:ASP",0.412));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:CYS",0.460));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:GLN",0.370));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:GLU",0.408));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:GLY",0.418));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:HIS",0.461));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:ILE",0.434));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:LEU",0.375));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:LYS",0.373));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:MET",0.439));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:PHE",0.428));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:SER",0.408));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:THR",0.455));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:TRP",0.480));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:TYR",0.455));
    this->accuracyMAE.insert(std::pair<std::string,double>("H:VAL",0.436));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:ALA",0.205));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:ARG",0.244));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:ASN",0.224));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:ASP",0.197));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:CYS",0.344));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:GLN",0.238));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:GLU",0.213));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:GLY",0.841));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:HIS",0.329));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:ILE",0.284));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:LEU",0.241));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:LYS",0.202));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:MET",0.291));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:PHE",0.339));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:PRO",0.270));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:SER",0.240));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:THR",0.254));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:TRP",0.326));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:TYR",0.317));
    this->accuracyMAE.insert(std::pair<std::string,double>("HA:VAL",0.284));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:ALA",2.315));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:ARG",2.410));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:ASN",2.568));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:ASP",2.341));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:CYS",2.690));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:GLN",2.252));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:GLU",2.265));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:GLY",2.642));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:HIS",2.753));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:ILE",2.991));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:LEU",2.448));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:LYS",2.209));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:MET",2.238));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:PHE",2.528));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:SER",2.618));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:THR",3.749));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:TRP",2.240));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:TYR",3.183));
    this->accuracyMAE.insert(std::pair<std::string,double>("N:VAL",3.141));



    this->accuracyR.insert(std::pair<std::string,double>("C:ALA",0.813));
    this->accuracyR.insert(std::pair<std::string,double>("C:ARG",0.814));
    this->accuracyR.insert(std::pair<std::string,double>("C:ASN",0.727));
    this->accuracyR.insert(std::pair<std::string,double>("C:ASP",0.715));
    this->accuracyR.insert(std::pair<std::string,double>("C:CYS",0.698));
    this->accuracyR.insert(std::pair<std::string,double>("C:GLN",0.819));
    this->accuracyR.insert(std::pair<std::string,double>("C:GLU",0.823));
    this->accuracyR.insert(std::pair<std::string,double>("C:GLY",0.640));
    this->accuracyR.insert(std::pair<std::string,double>("C:HIS",0.793));
    this->accuracyR.insert(std::pair<std::string,double>("C:ILE",0.766));
    this->accuracyR.insert(std::pair<std::string,double>("C:LEU",0.773));
    this->accuracyR.insert(std::pair<std::string,double>("C:LYS",0.814));
    this->accuracyR.insert(std::pair<std::string,double>("C:MET",0.771));
    this->accuracyR.insert(std::pair<std::string,double>("C:PHE",0.725));
    this->accuracyR.insert(std::pair<std::string,double>("C:PRO",0.583));
    this->accuracyR.insert(std::pair<std::string,double>("C:SER",0.701));
    this->accuracyR.insert(std::pair<std::string,double>("C:THR",0.685));
    this->accuracyR.insert(std::pair<std::string,double>("C:TRP",0.749));
    this->accuracyR.insert(std::pair<std::string,double>("C:TYR",0.749));
    this->accuracyR.insert(std::pair<std::string,double>("C:VAL",0.754));
    this->accuracyR.insert(std::pair<std::string,double>("CA:ALA",0.749));
    this->accuracyR.insert(std::pair<std::string,double>("CA:ARG",0.795));
    this->accuracyR.insert(std::pair<std::string,double>("CA:ASN",0.740));
    this->accuracyR.insert(std::pair<std::string,double>("CA:ASP",0.707));
    this->accuracyR.insert(std::pair<std::string,double>("CA:CYS",0.634));
    this->accuracyR.insert(std::pair<std::string,double>("CA:GLN",0.783));
    this->accuracyR.insert(std::pair<std::string,double>("CA:GLU",0.734));
    this->accuracyR.insert(std::pair<std::string,double>("CA:GLY",0.274));
    this->accuracyR.insert(std::pair<std::string,double>("CA:HIS",0.681));
    this->accuracyR.insert(std::pair<std::string,double>("CA:ILE",0.704));
    this->accuracyR.insert(std::pair<std::string,double>("CA:LEU",0.772));
    this->accuracyR.insert(std::pair<std::string,double>("CA:LYS",0.772));
    this->accuracyR.insert(std::pair<std::string,double>("CA:MET",0.655));
    this->accuracyR.insert(std::pair<std::string,double>("CA:PHE",0.708));
    this->accuracyR.insert(std::pair<std::string,double>("CA:PRO",0.466));
    this->accuracyR.insert(std::pair<std::string,double>("CA:SER",0.685));
    this->accuracyR.insert(std::pair<std::string,double>("CA:THR",0.795));
    this->accuracyR.insert(std::pair<std::string,double>("CA:TRP",0.760));
    this->accuracyR.insert(std::pair<std::string,double>("CA:TYR",0.626));
    this->accuracyR.insert(std::pair<std::string,double>("CA:VAL",0.756));
    this->accuracyR.insert(std::pair<std::string,double>("CB:ALA",0.317));
    this->accuracyR.insert(std::pair<std::string,double>("CB:ARG",0.446));
    this->accuracyR.insert(std::pair<std::string,double>("CB:ASN",0.324));
    this->accuracyR.insert(std::pair<std::string,double>("CB:ASP",0.265));
    this->accuracyR.insert(std::pair<std::string,double>("CB:CYS",0.227));
    this->accuracyR.insert(std::pair<std::string,double>("CB:GLN",0.547));
    this->accuracyR.insert(std::pair<std::string,double>("CB:GLU",0.512));
    this->accuracyR.insert(std::pair<std::string,double>("CB:HIS",0.126));
    this->accuracyR.insert(std::pair<std::string,double>("CB:ILE",0.307));
    this->accuracyR.insert(std::pair<std::string,double>("CB:LEU",0.346));
    this->accuracyR.insert(std::pair<std::string,double>("CB:LYS",0.357));
    this->accuracyR.insert(std::pair<std::string,double>("CB:MET",0.454));
    this->accuracyR.insert(std::pair<std::string,double>("CB:PHE",0.653));
    this->accuracyR.insert(std::pair<std::string,double>("CB:PRO",0.083));
    this->accuracyR.insert(std::pair<std::string,double>("CB:SER",0.101));
    this->accuracyR.insert(std::pair<std::string,double>("CB:THR",0.075));
    this->accuracyR.insert(std::pair<std::string,double>("CB:TRP",0.744));
    this->accuracyR.insert(std::pair<std::string,double>("CB:TYR",0.716));
    this->accuracyR.insert(std::pair<std::string,double>("CB:VAL",0.430));
    this->accuracyR.insert(std::pair<std::string,double>("H:ALA",0.694));
    this->accuracyR.insert(std::pair<std::string,double>("H:ARG",0.640));
    this->accuracyR.insert(std::pair<std::string,double>("H:ASN",0.572));
    this->accuracyR.insert(std::pair<std::string,double>("H:ASP",0.577));
    this->accuracyR.insert(std::pair<std::string,double>("H:CYS",0.591));
    this->accuracyR.insert(std::pair<std::string,double>("H:GLN",0.633));
    this->accuracyR.insert(std::pair<std::string,double>("H:GLU",0.631));
    this->accuracyR.insert(std::pair<std::string,double>("H:GLY",0.589));
    this->accuracyR.insert(std::pair<std::string,double>("H:HIS",0.652));
    this->accuracyR.insert(std::pair<std::string,double>("H:ILE",0.692));
    this->accuracyR.insert(std::pair<std::string,double>("H:LEU",0.714));
    this->accuracyR.insert(std::pair<std::string,double>("H:LYS",0.692));
    this->accuracyR.insert(std::pair<std::string,double>("H:MET",0.654));
    this->accuracyR.insert(std::pair<std::string,double>("H:PHE",0.691));
    this->accuracyR.insert(std::pair<std::string,double>("H:SER",0.567));
    this->accuracyR.insert(std::pair<std::string,double>("H:THR",0.526));
    this->accuracyR.insert(std::pair<std::string,double>("H:TRP",0.722));
    this->accuracyR.insert(std::pair<std::string,double>("H:TYR",0.701));
    this->accuracyR.insert(std::pair<std::string,double>("H:VAL",0.673));
    this->accuracyR.insert(std::pair<std::string,double>("HA:ALA",0.826));
    this->accuracyR.insert(std::pair<std::string,double>("HA:ARG",0.734));
    this->accuracyR.insert(std::pair<std::string,double>("HA:ASN",0.668));
    this->accuracyR.insert(std::pair<std::string,double>("HA:ASP",0.657));
    this->accuracyR.insert(std::pair<std::string,double>("HA:CYS",0.720));
    this->accuracyR.insert(std::pair<std::string,double>("HA:GLN",0.681));
    this->accuracyR.insert(std::pair<std::string,double>("HA:GLU",0.785));
    this->accuracyR.insert(std::pair<std::string,double>("HA:GLY",-0.601));
    this->accuracyR.insert(std::pair<std::string,double>("HA:HIS",0.688));
    this->accuracyR.insert(std::pair<std::string,double>("HA:ILE",0.763));
    this->accuracyR.insert(std::pair<std::string,double>("HA:LEU",0.787));
    this->accuracyR.insert(std::pair<std::string,double>("HA:LYS",0.781));
    this->accuracyR.insert(std::pair<std::string,double>("HA:MET",0.746));
    this->accuracyR.insert(std::pair<std::string,double>("HA:PHE",0.711));
    this->accuracyR.insert(std::pair<std::string,double>("HA:PRO",0.457));
    this->accuracyR.insert(std::pair<std::string,double>("HA:SER",0.677));
    this->accuracyR.insert(std::pair<std::string,double>("HA:THR",0.762));
    this->accuracyR.insert(std::pair<std::string,double>("HA:TRP",0.651));
    this->accuracyR.insert(std::pair<std::string,double>("HA:TYR",0.702));
    this->accuracyR.insert(std::pair<std::string,double>("HA:VAL",0.796));
    this->accuracyR.insert(std::pair<std::string,double>("N:ALA",0.675));
    this->accuracyR.insert(std::pair<std::string,double>("N:ARG",0.582));
    this->accuracyR.insert(std::pair<std::string,double>("N:ASN",0.655));
    this->accuracyR.insert(std::pair<std::string,double>("N:ASP",0.678));
    this->accuracyR.insert(std::pair<std::string,double>("N:CYS",0.640));
    this->accuracyR.insert(std::pair<std::string,double>("N:GLN",0.682));
    this->accuracyR.insert(std::pair<std::string,double>("N:GLU",0.618));
    this->accuracyR.insert(std::pair<std::string,double>("N:GLY",0.551));
    this->accuracyR.insert(std::pair<std::string,double>("N:HIS",0.682));
    this->accuracyR.insert(std::pair<std::string,double>("N:ILE",0.611));
    this->accuracyR.insert(std::pair<std::string,double>("N:LEU",0.724));
    this->accuracyR.insert(std::pair<std::string,double>("N:LYS",0.699));
    this->accuracyR.insert(std::pair<std::string,double>("N:MET",0.665));
    this->accuracyR.insert(std::pair<std::string,double>("N:PHE",0.662));
    this->accuracyR.insert(std::pair<std::string,double>("N:SER",0.596));
    this->accuracyR.insert(std::pair<std::string,double>("N:THR",0.514));
    this->accuracyR.insert(std::pair<std::string,double>("N:TRP",0.796));
    this->accuracyR.insert(std::pair<std::string,double>("N:TYR",0.595));
    this->accuracyR.insert(std::pair<std::string,double>("N:VAL",0.578));

}
void LARMORCA::initializeRandomShifts(){
    /* generated using: awk '{print "this->randomShifts.insert(std::pair<std::string,double>(\""$1":"$2"\","$3"));"}' randomcoil.dat */
    this->randomShifts.insert(std::pair<std::string,double>("C:ILE",176.4));
    this->randomShifts.insert(std::pair<std::string,double>("C:GLN",176.0));
    this->randomShifts.insert(std::pair<std::string,double>("C:GLY",174.9));
    this->randomShifts.insert(std::pair<std::string,double>("C:GLU",176.6));
    this->randomShifts.insert(std::pair<std::string,double>("C:CYS",174.6));
    this->randomShifts.insert(std::pair<std::string,double>("C:ASP",176.3));
    this->randomShifts.insert(std::pair<std::string,double>("C:SER",174.6));
    this->randomShifts.insert(std::pair<std::string,double>("C:LYS",176.6));
    this->randomShifts.insert(std::pair<std::string,double>("C:PRO",177.3));
    this->randomShifts.insert(std::pair<std::string,double>("C:HID",174.1));
    this->randomShifts.insert(std::pair<std::string,double>("C:HIE",174.1));
    this->randomShifts.insert(std::pair<std::string,double>("C:ASN",175.2));
    this->randomShifts.insert(std::pair<std::string,double>("C:VAL",176.3));
    this->randomShifts.insert(std::pair<std::string,double>("C:THR",174.7));
    this->randomShifts.insert(std::pair<std::string,double>("C:HIS",174.1));
    this->randomShifts.insert(std::pair<std::string,double>("C:TRP",176.1));
    this->randomShifts.insert(std::pair<std::string,double>("C:PHE",175.8));
    this->randomShifts.insert(std::pair<std::string,double>("C:ALA",177.8));
    this->randomShifts.insert(std::pair<std::string,double>("C:MET",173.6));
    this->randomShifts.insert(std::pair<std::string,double>("C:LEU",177.6));
    this->randomShifts.insert(std::pair<std::string,double>("C:ARG",176.3));
    this->randomShifts.insert(std::pair<std::string,double>("C:TYR",175.9));
    this->randomShifts.insert(std::pair<std::string,double>("CB:ILE",38.8));
    this->randomShifts.insert(std::pair<std::string,double>("CB:GLN",29.4));
    this->randomShifts.insert(std::pair<std::string,double>("CB:GLY",999));
    this->randomShifts.insert(std::pair<std::string,double>("CB:GLU",29.9));
    this->randomShifts.insert(std::pair<std::string,double>("CB:CYS",34.5));
    this->randomShifts.insert(std::pair<std::string,double>("CB:ASP",41.1));
    this->randomShifts.insert(std::pair<std::string,double>("CB:SER",63.8));
    this->randomShifts.insert(std::pair<std::string,double>("CB:LYS",33.1));
    this->randomShifts.insert(std::pair<std::string,double>("CB:PRO",33.3));
    this->randomShifts.insert(std::pair<std::string,double>("CB:HID",29.0));
    this->randomShifts.insert(std::pair<std::string,double>("CB:HIE",29.0));
    this->randomShifts.insert(std::pair<std::string,double>("CB:ASN",38.9));
    this->randomShifts.insert(std::pair<std::string,double>("CB:VAL",32.9));
    this->randomShifts.insert(std::pair<std::string,double>("CB:THR",69.8));
    this->randomShifts.insert(std::pair<std::string,double>("CB:HIS",29.0));
    this->randomShifts.insert(std::pair<std::string,double>("CB:TRP",29.6));
    this->randomShifts.insert(std::pair<std::string,double>("CB:PHE",39.6));
    this->randomShifts.insert(std::pair<std::string,double>("CB:ALA",19.1));
    this->randomShifts.insert(std::pair<std::string,double>("CB:MET",32.9));
    this->randomShifts.insert(std::pair<std::string,double>("CB:LEU",42.4));
    this->randomShifts.insert(std::pair<std::string,double>("CB:ARG",30.9));
    this->randomShifts.insert(std::pair<std::string,double>("CB:TYR",37.8));
    this->randomShifts.insert(std::pair<std::string,double>("CA:ILE",61.1));
    this->randomShifts.insert(std::pair<std::string,double>("CA:GLN",55.7));
    this->randomShifts.insert(std::pair<std::string,double>("CA:GLY",45.1));
    this->randomShifts.insert(std::pair<std::string,double>("CA:GLU",56.6));
    this->randomShifts.insert(std::pair<std::string,double>("CA:CYS",56.8));
    this->randomShifts.insert(std::pair<std::string,double>("CA:ASP",54.2));
    this->randomShifts.insert(std::pair<std::string,double>("CA:SER",58.3));
    this->randomShifts.insert(std::pair<std::string,double>("CA:LYS",56.2));
    this->randomShifts.insert(std::pair<std::string,double>("CA:PRO",58.05));
    this->randomShifts.insert(std::pair<std::string,double>("CA:HID",55.0));
    this->randomShifts.insert(std::pair<std::string,double>("CA:HIE",55.0));
    this->randomShifts.insert(std::pair<std::string,double>("CA:ASN",53.1));
    this->randomShifts.insert(std::pair<std::string,double>("CA:VAL",62.2));
    this->randomShifts.insert(std::pair<std::string,double>("CA:THR",61.8));
    this->randomShifts.insert(std::pair<std::string,double>("CA:HIS",55.0));
    this->randomShifts.insert(std::pair<std::string,double>("CA:TRP",57.5));
    this->randomShifts.insert(std::pair<std::string,double>("CA:PHE",57.7));
    this->randomShifts.insert(std::pair<std::string,double>("CA:ALA",52.5));
    this->randomShifts.insert(std::pair<std::string,double>("CA:MET",55.4));
    this->randomShifts.insert(std::pair<std::string,double>("CA:LEU",55.1));
    this->randomShifts.insert(std::pair<std::string,double>("CA:ARG",56.0));
    this->randomShifts.insert(std::pair<std::string,double>("CA:TYR",57.9));
    this->randomShifts.insert(std::pair<std::string,double>("N:ILE",119.9));
    this->randomShifts.insert(std::pair<std::string,double>("N:GLN",119.8));
    this->randomShifts.insert(std::pair<std::string,double>("N:GLY",108.8));
    this->randomShifts.insert(std::pair<std::string,double>("N:GLU",120.2));
    this->randomShifts.insert(std::pair<std::string,double>("N:CYS",118.7));
    this->randomShifts.insert(std::pair<std::string,double>("N:ASP",120.4));
    this->randomShifts.insert(std::pair<std::string,double>("N:SER",115.7));
    this->randomShifts.insert(std::pair<std::string,double>("N:LYS",120.4));
    this->randomShifts.insert(std::pair<std::string,double>("N:PRO",999));
    this->randomShifts.insert(std::pair<std::string,double>("N:HID",118.2));
    this->randomShifts.insert(std::pair<std::string,double>("N:HIE",118.2));
    this->randomShifts.insert(std::pair<std::string,double>("N:ASN",118.7));
    this->randomShifts.insert(std::pair<std::string,double>("N:VAL",119.2));
    this->randomShifts.insert(std::pair<std::string,double>("N:THR",113.6));
    this->randomShifts.insert(std::pair<std::string,double>("N:HIS",118.2));
    this->randomShifts.insert(std::pair<std::string,double>("N:TRP",121.3));
    this->randomShifts.insert(std::pair<std::string,double>("N:PHE",120.3));
    this->randomShifts.insert(std::pair<std::string,double>("N:ALA",123.8));
    this->randomShifts.insert(std::pair<std::string,double>("N:MET",119.6));
    this->randomShifts.insert(std::pair<std::string,double>("N:LEU",121.8));
    this->randomShifts.insert(std::pair<std::string,double>("N:ARG",120.5));
    this->randomShifts.insert(std::pair<std::string,double>("N:TYR",120.3));
    this->randomShifts.insert(std::pair<std::string,double>("H:ILE",8.0));
    this->randomShifts.insert(std::pair<std::string,double>("H:GLN",8.32));
    this->randomShifts.insert(std::pair<std::string,double>("H:GLY",8.3));
    this->randomShifts.insert(std::pair<std::string,double>("H:GLU",8.42));
    this->randomShifts.insert(std::pair<std::string,double>("H:CYS",8.375));
    this->randomShifts.insert(std::pair<std::string,double>("H:ASP",8.34));
    this->randomShifts.insert(std::pair<std::string,double>("H:SER",8.31));
    this->randomShifts.insert(std::pair<std::string,double>("H:LYS",8.29));
    this->randomShifts.insert(std::pair<std::string,double>("H:PRO",999));
    this->randomShifts.insert(std::pair<std::string,double>("H:HID",8.42));
    this->randomShifts.insert(std::pair<std::string,double>("H:HIE",8.42));
    this->randomShifts.insert(std::pair<std::string,double>("H:ASN",8.4));
    this->randomShifts.insert(std::pair<std::string,double>("H:VAL",8.03));
    this->randomShifts.insert(std::pair<std::string,double>("H:THR",8.15));
    this->randomShifts.insert(std::pair<std::string,double>("H:HIS",8.42));
    this->randomShifts.insert(std::pair<std::string,double>("H:TRP",8.25));
    this->randomShifts.insert(std::pair<std::string,double>("H:PHE",8.3));
    this->randomShifts.insert(std::pair<std::string,double>("H:ALA",8.24));
    this->randomShifts.insert(std::pair<std::string,double>("H:MET",8.28));
    this->randomShifts.insert(std::pair<std::string,double>("H:LEU",8.16));
    this->randomShifts.insert(std::pair<std::string,double>("H:ARG",8.23));
    this->randomShifts.insert(std::pair<std::string,double>("H:TYR",8.12));
    this->randomShifts.insert(std::pair<std::string,double>("HA:ILE",4.17));
    this->randomShifts.insert(std::pair<std::string,double>("HA:GLN",4.34));
    this->randomShifts.insert(std::pair<std::string,double>("HA:GLY",3.96));
    this->randomShifts.insert(std::pair<std::string,double>("HA:GLU",4.35));
    this->randomShifts.insert(std::pair<std::string,double>("HA:CYS",4.55));
    this->randomShifts.insert(std::pair<std::string,double>("HA:ASP",4.64));
    this->randomShifts.insert(std::pair<std::string,double>("HA:SER",4.47));
    this->randomShifts.insert(std::pair<std::string,double>("HA:LYS",4.32));
    this->randomShifts.insert(std::pair<std::string,double>("HA:PRO",4.42));
    this->randomShifts.insert(std::pair<std::string,double>("HA:HID",4.73));
    this->randomShifts.insert(std::pair<std::string,double>("HA:HIE",4.73));
    this->randomShifts.insert(std::pair<std::string,double>("HA:ASN",4.74));
    this->randomShifts.insert(std::pair<std::string,double>("HA:VAL",4.12));
    this->randomShifts.insert(std::pair<std::string,double>("HA:THR",4.35));
    this->randomShifts.insert(std::pair<std::string,double>("HA:HIS",4.73));
    this->randomShifts.insert(std::pair<std::string,double>("HA:TRP",4.66));
    this->randomShifts.insert(std::pair<std::string,double>("HA:PHE",4.62));
    this->randomShifts.insert(std::pair<std::string,double>("HA:ALA",4.32));
    this->randomShifts.insert(std::pair<std::string,double>("HA:MET",4.48));
    this->randomShifts.insert(std::pair<std::string,double>("HA:LEU",4.34));
    this->randomShifts.insert(std::pair<std::string,double>("HA:ARG",4.34));
    this->randomShifts.insert(std::pair<std::string,double>("HA:TYR",4.55));
}
