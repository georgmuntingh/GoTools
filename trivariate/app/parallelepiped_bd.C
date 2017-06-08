/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/Parallelepiped.h"
#include "GoTools/geometry/SplineSurface.h"
#include <memory>

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;


int main(int argc, char* argv[] )
{
  if (argc != 3)
      cout << "Usage: " << "infile outfile " << endl;

  // Open input surface file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  // Open outfile
  ofstream os(argv[2]);
  ALWAYS_ERROR_IF(os.bad(), "Bad or no output filename");

  // Read parallelepiped spacification
  double x, y, z;
  is >> x >> y >> z;
  Point corner(x,y,z);
  is >> x >> y >> z;
  Point dir_u(x,y,z);
  is >> x >> y >> z;
  Point dir_v(x,y,z);
  is >> x >> y >> z;
  Point dir_w(x,y,z);;
  double len_u, len_v, len_w;
  is >> len_u >> len_v >> len_w;

  Parallelepiped vol(corner, dir_u, dir_v, dir_w, len_u, len_v, len_w);
  shared_ptr<SplineVolume> vol2(vol.geometryVolume());

  vector<shared_ptr<SplineSurface> > bd_sfs = vol2->getBoundarySurfaces();
  for (size_t ki=0; ki<bd_sfs.size(); ++ki)
    {
      bd_sfs[ki]->writeStandardHeader(os);
      bd_sfs[ki]->write(os);
    }
}


  
