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

#ifndef LARMORCA_H
#define LARMORCA_H

#include <string>
#include <vector>
#include <map>


//#include "Molecule.hpp"

/* Forward Declaration  (only valid for pointers and references) */
 class Molecule;

class LARMORCA {
    private:
         std::map<std::string,bool> shiftAtoms;
         std::map<std::string,double> randomShifts;
         std::map<std::string,double> alphas;
         std::map<std::string,double> experimentalCS; 
         std::map<std::string,double> accuracyMAE;
         std::map<std::string,double> accuracyR;
    public:
        LARMORCA (Molecule *mol=NULL, const std::string fchemshift="");
        void initializeShiftAtoms();
        void initializeExpectedAccuracy();
        void initializeRandomShifts();
        void initializeAlpha();
        double getRandomShift(const std::string &key);
        double getExperimentalCS(const std::string &key);
        double getR(const std::string &key);
        double getMAE(const std::string &key);
        int getNShiftAtoms();
        void loadCSFile(const std::string fchemshift);
};
#endif
