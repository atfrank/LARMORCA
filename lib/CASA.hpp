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

#ifndef CASA_H
#define CASA_H

#include <string>
#include <vector>
#include <map>


//#include "Molecule.hpp"

/* Forward Declaration  (only valid for pointers and references) */
class Molecule;

class CASA {
    private:
         std::map<std::string,double> refSASA;
         std::map<std::string,double> SASA;
         std::map<std::string,double> polarity;
         std::map<std::string,double> radius;
    public:
        CASA (Molecule *mol=NULL);
        void initializeRefSASs();
        void initializeSASA();
        void initializePolarity();
        void initializeRadius();
        
        double getRefSASA(const std::string &key);
        double getSASA(const std::string &key);
        double getPolarity(const std::string &key);
        double getRadius(const std::string &key);
        
};
#endif
