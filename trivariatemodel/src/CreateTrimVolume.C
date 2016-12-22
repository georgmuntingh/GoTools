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

#include "GoTools/trivariatemodel/CreateTrimVolume.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/compositemodel/Loop.h"
#include "GoTools/compositemodel/ftEdge.h"
#include <fstream>
#include <cstdlib>

#define DEBUG

using std::vector;
using std::set;
using std::make_pair;

using namespace Go;

//==========================================================================
CreateTrimVolume::CreateTrimVolume(shared_ptr<SurfaceModel> model)
//==========================================================================
{
  model_ = model;
}

//==========================================================================
CreateTrimVolume::~CreateTrimVolume()
//==========================================================================
{

}

//==========================================================================
shared_ptr<ftVolume> CreateTrimVolume::fetchOneTrimVol()
//==========================================================================
{
  // Distinguish between boundary surfaces and trimming surfaces
  vector<shared_ptr<ftSurface> > bd_faces;
  vector<shared_ptr<ftSurface> > trim_faces;
  identifyBoundaryFaces(bd_faces, trim_faces);

  // Create boundary fitted volume
  shared_ptr<ParamVolume> vol = createVolume(bd_faces);

  // Trim volume with the remaining faces
  shared_ptr<ftVolume> trimvol = createTrimVolume(vol, trim_faces);
  return trimvol;
}

//==========================================================================
void 
CreateTrimVolume::identifyBoundaryFaces(vector<shared_ptr<ftSurface> >& bd_faces,
					vector<shared_ptr<ftSurface> >& trim_faces)
//==========================================================================
{
  // Initially all faces can be boundary faces
  bd_faces = model_->allFaces();

  // Remove faces connected to inner trim curves from the pool of boundary
  // faces
  identifyInnerTrim(bd_faces, trim_faces);

#ifdef DEBUG
  std::ofstream of1("bd_faces.g2");
for (size_t ki=0; ki<bd_faces.size(); ++ki)
  {
    bd_faces[ki]->surface()->writeStandardHeader(of1);
    bd_faces[ki]->surface()->write(of1);
  }
#endif

  // Identify faces belonging to the same underlying surface

  // Probably a too simple solution
  if (bd_faces.size() <= 6)
    return;
}

//==========================================================================
shared_ptr<ParamVolume> 
CreateTrimVolume::createVolume(vector<shared_ptr<ftSurface> >& bd_faces)
//==========================================================================
{
}

//==========================================================================
shared_ptr<ftVolume> 
CreateTrimVolume::createTrimVolume(shared_ptr<ParamVolume> vol, 
				   vector<shared_ptr<ftSurface> >& trim_faces)
//==========================================================================
{
}

//==========================================================================
void 
CreateTrimVolume::identifyInnerTrim(vector<shared_ptr<ftSurface> >& bd_faces,
				    vector<shared_ptr<ftSurface> >& trim_faces)
//==========================================================================
{
  // Collect all edges belonging to inner trimming loops
  vector<shared_ptr<ftEdgeBase> > trim_edg;
  for (size_t ki=0; ki<bd_faces.size(); ++ki)
    {
      int nmb_loops = bd_faces[ki]->nmbBoundaryLoops();
      if (nmb_loops <= 1)
	continue;

      for (int ka=1; ka<nmb_loops; ++ka)
	{
	  shared_ptr<Loop> loop = bd_faces[ki]->getBoundaryLoop(ka);
	  vector<shared_ptr<ftEdgeBase> > tmp_edg = loop->getEdges();
	  trim_edg.insert(trim_edg.end(), tmp_edg.begin(), tmp_edg.end());
	}
    }

  // Classify surfaces with outer boundary loop edges being connected
  // to the inner trim edges as trim faces
  for (size_t ki=0; ki<bd_faces.size(); )
    {
      // Fetch outer boundary loop
      shared_ptr<Loop> loop = bd_faces[ki]->getBoundaryLoop(0);
      int nmb = loop->size();
      int ka = 0;
      for (; ka<nmb; ++ka)
	{
	  shared_ptr<ftEdgeBase> edg = loop->getEdge(ka);
	  size_t kj = 0;
	  for (; kj<trim_edg.size(); ++kj)
	    {
	      if (trim_edg[kj].get() == edg.get() ||
		  trim_edg[kj].get() == edg->twin())
		break;
	    }
	  if (kj < trim_edg.size())
	    {
	      trim_faces.push_back(bd_faces[ki]);
	      bd_faces.erase(bd_faces.begin()+ki);
	      break;
	    }
	}
      if (ka >= nmb)
	++ki;
    }
}
