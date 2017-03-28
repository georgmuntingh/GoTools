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

#include "GoTools/compositemodel/RegularizeFace.h"
#include "GoTools/compositemodel/ModifyFaceSet.h"
#include "GoTools/compositemodel/RegularizeUtils.h"
#include "GoTools/compositemodel/ftEdge.h"
#include <fstream>
#include <cstdlib>

#define DEBUG

using std::vector;
using std::set;
using std::make_pair;

using namespace Go;

//==========================================================================
    ModifyFaceSet::ModifyFaceSet(shared_ptr<SurfaceModel> model)
      : model_(model)
//==========================================================================
{
}

//==========================================================================
ModifyFaceSet::~ModifyFaceSet()
//==========================================================================
{

}

//==========================================================================
shared_ptr<SurfaceModel> ModifyFaceSet::getModifiedModel(int& nmb)
//==========================================================================
{
  nmb = divide();
  return model_;
}

  // Collect all faces
//==========================================================================
int ModifyFaceSet::divide()
//==========================================================================
{
  // Collect all faces
  vector<shared_ptr<ftSurface> > faces = model_->allFaces();
  int nmb_faces = (int)faces.size();

  tpTolerances tptol = model_->getTolerances();

  // Fetch all sharp edges
  vector<ftEdge*> sharp_edges = fetchSharpEdges();

  int nmb = 0;
  for (size_t ki=0; ki<sharp_edges.size(); ++ki)
    {
      nmb++;
      shared_ptr<Vertex> vxs[2];
      sharp_edges[ki]->getVertices(vxs[0], vxs[1]);
      ftEdge* curr_edges[2];
      curr_edges[0] = curr_edges[1] = sharp_edges[ki];
      ftSurface* adj_face[2];
      adj_face[0] = adj_face[1] = NULL;
      for (size_t ka=0; ka<faces.size(); ++ka)
	{
	  // Split faces until a full edge loop including the sharp edge
	  // is obtained. The computatition is not expected to involve all
	  // faces, but the loop is limited to avoid to go infinite
	  bool finish = false;
	  ftEdge* next_edge[2];
	  next_edge[0] = next_edge[1] = NULL;
	  double angle[2];
	  for (int kj=0; kj<2; ++kj)
	    {
	      if (curr_edges[kj])
		adj_face[kj] = fetchNextFace(curr_edges[kj], vxs[kj].get(),
					     tptol.bend, next_edge[kj], angle[kj]);
	    }
	  if (adj_face[0] == adj_face[1])
	    {
	      adj_face[1] = NULL;   // Same face, split only once
	      finish = true;
	    }

	  for (int kj=0; kj<2; ++kj)
	    {
	      if (adj_face[kj] && M_PI-angle[kj] > tptol.bend)
		{
		  // Regularize face
		  shared_ptr<ftSurface> curr = 
		    model_->fetchAsSharedPtr(adj_face[kj]);
		  if (!curr.get())
		    {
		      finish = true;
		      continue;
		    }
		  RegularizeFace regularize(curr, model_);
		  vector<shared_ptr<Vertex> > vx_pri;
		  vx_pri.push_back(vxs[kj]);
		  regularize.setPriVx(vx_pri);

		  vector<shared_ptr<Vertex> > corner = 
		    curr->getCornerVertices(tptol.bend);
		  if (corner.size() > 4)
		    regularize.setDivideInT(false);
		  vector<shared_ptr<ftSurface> > faces2 = 
		    regularize.getRegularFaces();
#ifdef DEBUG
		  std::ofstream of2("post_regface.g2");
		  for (size_t kr=0; kr<faces2.size(); ++kr)
		    {
		      faces2[kr]->surface()->writeStandardHeader(of2);
		      faces2[kr]->surface()->write(of2);
		    }
#endif
		  // Identify edge continuation. First fetch candidates
		  vector<shared_ptr<ftEdge> > cand_edg;
		  for (size_t kr=0; kr<faces2.size(); ++kr)
		    for (size_t kh=kr+1; kh<faces2.size(); ++kh)
		      {
			vector<shared_ptr<ftEdge> > tmp_edg = 
			  faces2[kr]->getCommonEdges(faces2[kh].get());
			for (size_t kb=0; kb<tmp_edg.size(); ++kb)
			  {
			    if (curr_edges[kj]->commonVertex(tmp_edg[kb].get()))
			      cand_edg.push_back(tmp_edg[kb]);
			  }
		      }
		  if (cand_edg.size() == 1)
		    next_edge[kj] = cand_edg[0].get();
 		}

	      // Fetch the next vertex in the chain
	      if (next_edge[kj])
		{
		  shared_ptr<Vertex> common_vx = 
		    next_edge[kj]->getCommonVertex(curr_edges[kj]);
		  vxs[kj] = next_edge[kj]->getOtherVertex(common_vx.get());
		}

	      curr_edges[kj] = next_edge[kj];
	      nmb++;
	    }

	  if (adj_face[0] == adj_face[1])
	    nmb--;
	  if (finish || 
	      vxs[0]->getVertexPoint().dist(vxs[1]->getVertexPoint()) < tptol.neighbour)
	    break;
	}
    }
  return nmb;
}

//==========================================================================
vector<ftEdge*>  ModifyFaceSet::fetchSharpEdges()
//==========================================================================
{
  vector<ftEdge*> sharp_edges;
  model_->getCorners(sharp_edges);

  // Remove concex edges
  for (size_t ki=0; ki<sharp_edges.size(); )
    {
      if (!sharp_edges[ki]->twin())
	{
	  sharp_edges.erase(sharp_edges.begin()+ki);
	  continue;
	}
      int nmb_samples = 10;
      int nmb_concave = 0;
      int nmb_convex = 0;
      double t1 = sharp_edges[ki]->tMin();
      double t2 = sharp_edges[ki]->tMax();
      double del = (t2 - t1)/(double)(nmb_samples-1);
      double par;
      int kr;
      for (kr=0, par=t1; kr<nmb_samples; ++kr, par+=del)
	{
	  Point pos1 = sharp_edges[ki]->point(par);
	  Point norm1 = sharp_edges[ki]->normal(par);
	  Point tan1 = sharp_edges[ki]->tangent(par);
	  Point vec1 = tan1.cross(norm1);

	  double clo_t, clo_dist;
	  Point clo_pt;
	  sharp_edges[ki]->twin()->closestPoint(pos1, clo_t, clo_pt,
					       clo_dist);
	  Point norm2 = sharp_edges[ki]->twin()->normal(clo_t);
	  Point tan2 = sharp_edges[ki]->twin()->tangent(clo_t);
	  Point vec2 = tan2.cross(norm2);
	  double ang2 = vec1.angle(norm2);
	  double ang3 = vec2.angle(norm1);
	  if (ang2 + ang3 > M_PI)
	    nmb_concave++;
	  else
	    nmb_convex++;
	  int stop_break;
	}
      if (nmb_convex >= nmb_concave)
	  sharp_edges.erase(sharp_edges.begin()+ki);
      else
	++ki;
    }

#ifdef DEBUG
  std::ofstream of("convex_crvs.g2");
  for (size_t ki=0; ki<sharp_edges.size(); ++ki)
    {
      shared_ptr<ParamCurve> cv = 
	shared_ptr<ParamCurve>(sharp_edges[ki]->geomCurve()->subCurve(sharp_edges[ki]->tMin(),
								      sharp_edges[ki]->tMax()));
      cv->geometryCurve()->writeStandardHeader(of);
      cv->geometryCurve()->write(of);
    }
#endif
  return sharp_edges;
}

//==========================================================================
ftSurface*  ModifyFaceSet::fetchNextFace(ftEdge* edge, Vertex* vx, double angtol,
					 ftEdge*& next_edge, double& angle)
//==========================================================================
{
  next_edge = NULL;  // Initially
  ftFaceBase *edg_faces[2];
  edg_faces[0] = edge->face();
  edg_faces[1] = edge->twin()->face();
  vector<ftSurface*> vx_faces = vx->faces();
  Point vxpoint = vx->getVertexPoint();
  shared_ptr<Vertex> vx2 = edge->getOtherVertex(vx);
  Point vxdir = vxpoint - vx2->getVertexPoint();
  angle = M_PI;
  for (size_t kr=0; kr<vx_faces.size();)
    {
      int kh;
      for (kh=0; kh<2; ++kh)
	if (vx_faces[kr] == edg_faces[kh])
	  break;
      if (kh < 2)
	vx_faces.erase(vx_faces.begin()+kr);
      else
	++kr;
    }
  if (vx_faces.size() == 1)
    {
#ifdef DEBUG
      std::ofstream of1("adj_faces.g2");
      vx_faces[0]->surface()->writeStandardHeader(of1);
      vx_faces[0]->surface()->write(of1);
#endif

      // Check angles between associated edges
      Point tan1 = edge->tangent(edge->parAtVertex(vx));
      vector<ftEdge*> edgs = vx->getEdges(vx_faces[0]);
      vector<double> scpr(edgs.size(), 0.0);
      for (size_t kr=0; kr<edgs.size(); ++kr)
	{
	  Point tan2 = edgs[kr]->tangent(edgs[kr]->parAtVertex(vx));
	  double ang = tan1.angle(tan2);
	  angle = std::min(angle, ang);
	  if (M_PI-ang < angtol)
	    {
	      angle = ang;   // Continue along edge chain without face split

	      // Find direction of edge continuation
	      shared_ptr<Vertex> vx3 = edgs[kr]->getOtherVertex(vx);
	      Point vxdir3 = vx3->getVertexPoint() - vxpoint;
	      scpr[kr] = vxdir*vxdir3;
	    }
	}
      double maxpr = 0.0;
      int ix = -1;
      for (size_t kr=0; kr<scpr.size(); ++kr)
	{
	  if (scpr[kr] > maxpr)
	    {
	      maxpr = scpr[kr];
	      ix = (int)kr;
	    }
	}
      if (ix >= 0)
	next_edge = edgs[ix];
      return vx_faces[0];
    }
  return NULL;
}


