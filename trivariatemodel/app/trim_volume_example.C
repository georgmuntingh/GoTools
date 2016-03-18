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
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/ftVolumeTools.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/Utils.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::vector;

int main(int argc, char* argv[] )
{
  if (argc != 2)
      cout << "Usage: " << "g2 Brep file" << endl;

  ifstream infile(argv[1]);
  ALWAYS_ERROR_IF(infile.bad(), "Bad or no input filename");

  // The tolerances must be set according to the properties of the model.
  // The neighbour tolerance must be smaller than the smallest entity in the
  // model, but larger than the largest gap.
  // The gap tolerance must be smaller than the neighbour tolerance
  double gap = 0.001; //0.001;
  double neighbour = 0.01; //0.01;
  double kink = 0.01;
  double approxtol = 0.001;
  int degree = 3;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromG2(infile);

  shared_ptr<SurfaceModel> sfmodel = 
    shared_ptr<SurfaceModel>(dynamic_cast<SurfaceModel*>(model));
  if (!sfmodel.get())
    {
      std::cout << "No input model read" << std::endl;
      exit(-1);
    }
 
  if (sfmodel->nmbBoundaries() > 0)
    {
      std::cout << "Not a brep solid. Consider increasing the neighbour tolerance" << std::endl;
      exit(-1);
    }
      
  bool isOK = sfmodel->checkShellTopology();
  std::cout << "Shell topology: " << isOK << std::endl;

  // A spline volume will be created, but no attempt to orient it in space or adapt to
  // certain surfaces is performed
  shared_ptr<ftVolume> vol(new ftVolume(sfmodel));

  // Check if the volume can be represented without any trimming
  if (vol->isRegularized())
    {
      int degree = 3;  // cubic spline volume
      bool untrimmed = vol->untrimRegular(degree);  // An approximation is performed
      std::cout << "Untrimming performed: " << untrimmed << std::endl;
    }

  // Check surface category as volume boundary surfaces
  // First fetch outer boundary
  shared_ptr<SurfaceModel> shell = vol->getOuterShell();

  // Number of surfaces in shell
  int nmb = shell->nmbEntities();
  for (int ki=0; ki<nmb; ++ki)
    {
      // Fetch face
      shared_ptr<ftSurface> face = shell->getFace(ki);

      // Check category
      // = -1 : not a volume boundary (trimming face)
      // = 0 : umin (the surface may be trimmed, it may not cover the entire volume boundar,
      // it may not have the same representation as the volume boundary surface, 
      // but it is coincident
      // = 1 : umax
      // = 2 : vmin
      // = 3 : vmax
      // = 4 : wmin
      // = 5 : wmax
      int bd_status = ftVolumeTools::boundaryStatus(vol, face, gap);

      // Fetch associated surface
      shared_ptr<ParamSurface> surf = face->surface();

      // Parameter domain containing the surface. In case of trimming, it might be larger
      RectDomain dom = surf->containingDomain();

      // Define grid in the internal of the surface (non-trimmed)
      int nmb_u = 5;
      int nmb_v = 5;
      vector<double> param_u(nmb_u);
      vector<double> param_v(nmb_v);
      double del_u = (dom.umax() - dom.umin())/(double)(nmb_u+1);
      double del_v = (dom.vmax() - dom.vmin())/(double)(nmb_v+1);
      int kj, kr;
      double par;
      for (kj=0, par=dom.umin()+del_u; kj<nmb_u; ++kj, par+=del_u)
	param_u[kj] = par;
      for (kj=0, par=dom.vmin()+del_v; kj<nmb_v; ++kj, par+=del_v)
	param_v[kj] = par;

      // Check if the grid points are inside the trimmed surface. In that case
      // evaluate basis functions, point and surface normal. The surface is oriented
      // such that the normal points out of the volume.
      // Fetch associated spline surface. If the surface is trimmed, the underlying
      // non-trimmed surface is returned. If the surface is not a spline surface (this
      // can happen in theory), a null pointer is returned.
      SplineSurface* spline_surf = surf->asSplineSurface();
      for (kj=0; kj<nmb_v; ++kj)
	for (kr=0; kr<nmb_u; ++kr)
	  {
	  }

      int stop_debug = 1;
    }
  

}

