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
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/intersections/Identity.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Torus.h"

#include <fstream>

#define DEBUG

using std::vector;
using namespace Go;

//===========================================================================
vector<shared_ptr<ParamSurface> > 
SurfaceModelUtils::checkClosedFaces(shared_ptr<ParamSurface> surface, double tol)
//===========================================================================
{
  vector<shared_ptr<ParamSurface> > sfs;

#ifdef DEBUG
  std::ofstream of("close_sf.g2");
  surface->writeStandardHeader(of);
  surface->write(of);
#endif

  // Fetch non-trimmed surface to test
  shared_ptr<ParamSurface> sf;
  RectDomain dom = surface->containingDomain();
  shared_ptr<BoundedSurface> bd_sf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surface);
  if (bd_sf.get())
    {
      shared_ptr<ParamSurface> tmp = bd_sf->underlyingSurface();
      vector<shared_ptr<ParamSurface> > sub_sfs;
      try {
	sub_sfs = tmp->subSurfaces(dom.umin(), dom.vmin(), 
				   dom.umax(), dom.vmax());
      }
      catch (...)
	{
	  sf = tmp;
	}

      if (sub_sfs.size() > 1)
	sf = tmp;
      else if (sub_sfs.size() == 1)
	sf = sub_sfs[0];
    }
  else
    sf = surface;
  
  // Fetch opposite boundary curves and check coincidence
  vector<shared_ptr<ParamCurve> > cvs_u1 = sf->constParamCurves(dom.umin(), false);
  vector<shared_ptr<ParamCurve> > cvs_u2 = sf->constParamCurves(dom.umax(), false);
  vector<shared_ptr<ParamCurve> > cvs_v1 = sf->constParamCurves(dom.vmin(), true);
  vector<shared_ptr<ParamCurve> > cvs_v2 = sf->constParamCurves(dom.vmax(), true);

  Identity ident;
  int coinc1 = ident.identicalCvs(cvs_u1[0], cvs_u1[0]->startparam(), cvs_u1[0]->endparam(),
				  cvs_u2[0], cvs_u2[0]->startparam(), cvs_u2[0]->endparam(),
				  tol);
  int coinc2 = ident.identicalCvs(cvs_v1[0], cvs_v1[0]->startparam(), cvs_v1[0]->endparam(),
				  cvs_v2[0], cvs_v2[0]->startparam(), cvs_v2[0]->endparam(),
				  tol);
  vector<shared_ptr<ParamSurface> > sub_sfs1;
  if (coinc1)
    {
      // Split in the first parameter direction
      vector<shared_ptr<ParamSurface> > sub_sfs2;
      double mid = 0.5*(dom.umin()+dom.umax());
      try {
	sub_sfs1 = surface->subSurfaces(dom.umin(), dom.vmin(), 
					mid, dom.vmax());
      }
      catch (...)
	{
	  sfs.push_back(surface);
	  return sfs;
	}
	
      try {
	sub_sfs2 = surface->subSurfaces(mid, dom.vmin(), 
					dom.umax(), dom.vmax());
      }
      catch(...)
	{
	  sfs.push_back(surface);
	  return sfs;
	}
      sub_sfs1.insert(sub_sfs1.end(), sub_sfs2.begin(), sub_sfs2.end());
    }
  else
    sub_sfs1.push_back(surface);
	
  if (coinc2)
    {
      for (size_t ki=0; ki<sub_sfs1.size(); ++ki)
	{
	  // Split in the second parameter direction
	  RectDomain dom2 = sub_sfs1[ki]->containingDomain();
	  vector<shared_ptr<ParamSurface> > sub_sfs2;
	  double mid = 0.5*(dom2.vmin()+dom2.vmax());
	  try {
	    sub_sfs2 = sub_sfs1[ki]->subSurfaces(dom2.umin(), dom2.vmin(), 
						 dom2.umax(), mid);
	  }
	  catch (...)
	    {
	      return sub_sfs1;
	    }
	    
	  sfs.insert(sfs.end(), sub_sfs2.begin(), sub_sfs2.end());
	  sub_sfs2.clear();
	  try {
	    sub_sfs2 = sub_sfs1[ki]->subSurfaces(dom2.umin(), mid,
						 dom2.umax(), dom2.vmax());
	  }
	  catch (...)
	    {
	      return sub_sfs1;
	    }
	  sfs.insert(sfs.end(), sub_sfs2.begin(), sub_sfs2.end());
	}
#ifdef DEBUG
      for (size_t kj=0; kj<sfs.size(); ++kj)
	{
	  sfs[kj]->writeStandardHeader(of);
	  sfs[kj]->write(of);
	}
#endif
      return sfs;
    }
  else
    {
#ifdef DEBUG
      for (size_t kj=0; kj<sub_sfs1.size(); ++kj)
	{
	  sub_sfs1[kj]->writeStandardHeader(of);
	  sub_sfs1[kj]->write(of);
	}
#endif
      return sub_sfs1;
    }
}

//===========================================================================
void
SurfaceModelUtils::sameUnderlyingSurf(vector<shared_ptr<ftSurface> >& sf_set,
				      double tol, double angtol,
				      vector<vector<shared_ptr<ftSurface> > >& faces,
				      vector<shared_ptr<ParamSurface> >& under_sfs)
//===========================================================================
{
  // Extract bounded surfaces
  vector<shared_ptr<ParamSurface> > cand_sfs;
  for (size_t ki=0; ki<sf_set.size(); ++ki)
    {
      shared_ptr<ParamSurface> sf = sf_set[ki]->surface();
      if (sf->instanceType() == Class_BoundedSurface)
	cand_sfs.push_back(sf);
    }

  // Make pairwise check of candidate surfaces for identical underlying surfaces
  for (size_t ki=0; ki<cand_sfs.size();)
    {
      size_t incr = 1;
      shared_ptr<BoundedSurface> bd_sf1 =
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(cand_sfs[ki]);
      shared_ptr<ParamSurface> under1 = bd_sf1->underlyingSurface();
      for (size_t kj=ki+1; kj<cand_sfs.size(); ++kj)
	{
	  shared_ptr<BoundedSurface> bd_sf2 =
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(cand_sfs[kj]);
	  shared_ptr<ParamSurface> under2 = bd_sf2->underlyingSurface();
	  bool same = false;
	  if (under1.get() == under2.get())
	    // Same underlying surface
	    same = true;
	  else
	    {
	      // Check for equality of elementary surfaces
	      shared_ptr<ElementarySurface> elem1 = 
		dynamic_pointer_cast<ElementarySurface, ParamSurface>(under1);
	      if (!elem1.get() && under1->isSpline())
		{
		  shared_ptr<SplineSurface> spline = 
		    dynamic_pointer_cast<SplineSurface, ParamSurface>(under1);
		  if (spline.get())
		    elem1 = spline->getElementarySurface();
		}
	      shared_ptr<ElementarySurface> elem2 = 
		dynamic_pointer_cast<ElementarySurface, ParamSurface>(under2);
	      if (!elem2.get() && under2->isSpline())
		{
		  shared_ptr<SplineSurface> spline = 
		    dynamic_pointer_cast<SplineSurface, ParamSurface>(under2);
		  if (spline.get())
		    elem2 = spline->getElementarySurface();
		}
	      
	      if (elem1.get() && elem2.get())
		{
		  // Both surfaces are elementary. Check for equality
		  if (elem1->instanceType() == Class_Plane &&
		      elem2->instanceType() == Class_Plane)
		    {
		      shared_ptr<Plane> plane1 = 
			dynamic_pointer_cast<Plane, ElementarySurface>(elem1);
		      shared_ptr<Plane> plane2 = 
			dynamic_pointer_cast<Plane, ElementarySurface>(elem2);
		      Point pt1 = plane1->getPoint();
		      Point pt2 = plane2->getPoint();
		      Point norm1 = plane1->getNormal();
		      Point norm2 = plane2->getNormal();
		      double ang = norm1.angle(norm2);
		      if (ang < angtol || M_PI-ang < angtol)
			{
			  double len = fabs((pt2 - pt1)*norm1);
			  if (len < tol)
			    same = true;
			}
		    }
		  else if (elem1->instanceType() == Class_Cylinder &&
			   elem2->instanceType() == Class_Cylinder)
		    {
		      shared_ptr<Cylinder> cyl1 = 
			dynamic_pointer_cast<Cylinder, ElementarySurface>(elem1);
		      shared_ptr<Cylinder> cyl2 = 
			dynamic_pointer_cast<Cylinder, ElementarySurface>(elem2);
		      Point pt1 = cyl1->getLocation();
		      Point pt2 = cyl2->getLocation();
		      Point axis1 = cyl1->getAxis();
		      Point axis2 = cyl2->getAxis();
		      double rad1 = cyl1->getRadius();
		      double rad2 = cyl2->getRadius();
		      double ang = axis1.angle(axis2);
		      if (fabs(rad1-rad2) < tol && 
			  (ang < angtol || M_PI-ang < angtol))
			{
			  double len = fabs((pt2 - pt1)*axis1);
			  if (len < tol)
			    same = true;
			}
		    }
		  else if (elem1->instanceType() == Class_Cone &&
			   elem2->instanceType() == Class_Cone)
		    {
		      shared_ptr<Cone> cone1 = 
			dynamic_pointer_cast<Cone, ElementarySurface>(elem1);
		      shared_ptr<Cone> cone2 = 
			dynamic_pointer_cast<Cone, ElementarySurface>(elem2);
		      Point pt1 = cone1->getLocation();
		      Point pt2 = cone2->getLocation();
		      Point axis1 = cone1->getAxis();
		      Point axis2 = cone2->getAxis();
		      double rad1 = cone1->getRadius();
		      double rad2 = cone2->getRadius();
		      double angle1 = cone1->getConeAngle();
		      double angle2 = cone2->getConeAngle();
		      double ang = axis1.angle(axis2);
		      if (fabs(angle1-angle2) < angtol && 
			  (ang < angtol || M_PI-ang < angtol))
			{
			  double len = fabs((pt2 - pt1)*axis1);
			  double d = pt1.dist(pt2);
			  double tanalpha = tan(angle1);
			  if (fabs(tanalpha) > tol)
			    {
			      double d1 = rad1/tanalpha;
			      double d2 = rad2/tanalpha - d;
			      if (len < tol && fabs(d1-d2) < tol)
				same = true;
			    }
			  else if (fabs(rad1-rad2) < tol)
			    same = true;
			}
		    }
		  else if (elem1->instanceType() == Class_Sphere &&
			   elem2->instanceType() == Class_Sphere)
		    {
		      shared_ptr<Sphere> sphere1 = 
			dynamic_pointer_cast<Sphere, ElementarySurface>(elem1);
		      shared_ptr<Sphere> sphere2 = 
			dynamic_pointer_cast<Sphere, ElementarySurface>(elem2);
		      Point pt1 = sphere1->getLocation();
		      Point pt2 = sphere2->getLocation();
		      double rad1 = sphere1->getRadius();
		      double rad2 = sphere2->getRadius();
		      double d = pt1.dist(pt2);
		      if (d < tol && fabs(rad1-rad2) < tol)
			same = true;
		    }
		  else if (elem1->instanceType() == Class_Torus &&
			   elem2->instanceType() == Class_Torus)
		    {
		      shared_ptr<Torus> tor1 = 
			dynamic_pointer_cast<Torus, ElementarySurface>(elem1);
		      shared_ptr<Torus> tor2 = 
			dynamic_pointer_cast<Torus, ElementarySurface>(elem2);
		      Point pt1 = tor1->getLocation();
		      Point pt2 = tor2->getLocation();
		      double radmin1 = tor1->getMinorRadius();
		      double radmin2 = tor2->getMinorRadius();
		      double radmax1 = tor1->getMajorRadius();
		      double radmax2 = tor2->getMajorRadius();
		      double d = pt1.dist(pt2);
		      Point x1, y1, z1, x2, y2, z2;
		      tor1->getCoordinateAxes(x1, y1, z1);
		      tor2->getCoordinateAxes(x2, y2, z2);
		      double ang1 = x1.angle(x2);
		      double ang2 = x1.angle(x2);
		      double ang3 = x1.angle(x2);
		      if (d < tol && fabs(radmin1-radmin2) < tol &&
			  fabs(radmax1-radmax2) < tol &&
			  (ang1 < angtol || M_PI-ang1 < angtol) &&
			  (ang2 < angtol || M_PI-ang2 < angtol) &&
			  (ang3 < angtol || M_PI-ang3 < angtol))
			same = true;
		    }
		}
	    }
	  if (!same)
	    {
	      // Check coincidence
	      Identity ident;
	      int coinc = ident.identicalSfs(under1, under2, tol);
	      if (coinc >= 1)
		same = true;
	    }
	  if (same)
	    {
	      std::swap(cand_sfs[ki+1], cand_sfs[kj]);
	      incr++;
	    }
	}
      if (incr > 1)
	{
	  // Identical underlying surfaces are found
	  vector<shared_ptr<ftSurface> > curr_faces;
	  for (size_t kj=ki; kj<ki+incr; ++kj)
	    {
	      for (size_t kr=0; kr<sf_set.size(); ++kr)
		{
		  if (sf_set[kr]->surface().get() == cand_sfs[kj].get())
		    {
		      curr_faces.push_back(sf_set[kr]);
		      break;
		    }
		}
	    }
	  shared_ptr<BoundedSurface> bd_sf =
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(cand_sfs[ki]);
	  if (bd_sf.get())
	    {
	      faces.push_back(curr_faces);
	      shared_ptr<ParamSurface> curr_under = bd_sf->underlyingSurface();
	      ElementarySurface *curr_elem = curr_under->elementarySurface();
	      if (curr_elem)
		{
		  // Construct underlying surface with a possibly extend 
		  // parameter domain covering all associated faces
		  // Case dependent
		  curr_under = extendedUnderlyingSurface(curr_faces, tol, angtol);
		}
	      under_sfs.push_back(curr_under);
	    }
	}
      ki += incr;
    }
}

//===========================================================================
shared_ptr<ParamSurface>
SurfaceModelUtils::extendedUnderlyingSurface(vector<shared_ptr<ftSurface> >& sf_set,
					     double tol, double angtol)
//===========================================================================
{
  shared_ptr<ParamSurface> surf;  // Resulting surface

  ElementarySurface *elem1 = sf_set[0]->surface()->elementarySurface();
  if (!elem1)
    return surf; // Function only applicable for elementary surfaces. No result

  // Check that all surfaces is of the same type
  // Check also equality of surface descriptions
  RectDomain dom = elem1->containingDomain();
  RectDomain dom1 = dom;
  Point loc1 = elem1->location();
  Point dir1 = elem1->direction();
  for (size_t ki=1; ki<sf_set.size(); ++ki)
    {
      ElementarySurface *elem2 = sf_set[ki]->surface()->elementarySurface();
      if (!elem2 || elem1->instanceType() != elem2->instanceType())
	return surf; 

      RectDomain dom2 = elem2->containingDomain();
      Point loc2 = elem2->location();
      Point dir2 = elem2->direction();
      double ang = (dir1.dimension() == 0) ? 0.0 : dir1.angle(dir2);
      if (ang > angtol && M_PI-ang > angtol)
	continue;   // Not the same surface

      // Case distinction
      if (elem1->instanceType() == Class_Plane)
	{
	  double len1 = fabs((loc1-loc2)*dir1);
	  if (len1 > tol)
	    continue;   // Not the same surface

	  Plane *plane1 = dynamic_cast<Plane*>(elem1);
	  Plane *plane2 = dynamic_cast<Plane*>(elem2);
	  if (plane1 == NULL || plane2 == NULL)
	    return surf;  // Something is wrong

	  double len2 = loc1.dist(loc2);
	  Point axis1_1, axis1_2, axis2_1, axis2_2;
	  plane1->getSpanningVectors(axis1_1, axis1_2);
	  plane2->getSpanningVectors(axis2_2, axis2_2);
	  if (plane1->isSwapped())
	      std::swap(axis1_1, axis1_2);
	  if (plane2->isSwapped())
	      std::swap(axis2_1, axis2_2);
	  double ang2 = axis1_1.angle(axis1_2);
	  if (len2 > tol || ang2 > angtol)
	    {
	      // Parameterization of the two planes differ.
	      // Modify parameter domain of surface 2 to match that of surface 1
	      Vector2D low = dom2.lowerLeft();
	      Vector2D high = dom2.upperRight();
	      if (len2 > tol)
		{
		  // Move domain
		  Point diff = loc2 - loc1;
		  low[0] += diff[0];
		  high[0] += diff[0];
		  low[1] += diff[1];
		  high[1] += diff[1];
		}
	      if (ang2 > angtol)
		{
		  // Extend domain to make sure to cover the rotated domain
		  double fac = 2.0/sqrt(2.0);
		  Vector2D mid = 0.5*(high + low);
		  low = mid - fac*(mid - low);
		  high += high + fac*(high - mid);
		}
	      RectDomain dom3(low, high);
	      dom1.addUnionWith(dom3);
	    }
	  else
	    {
	      // Simply extend the initial domain
	      dom1.addUnionWith(dom2);
	    }
	}
      else if (elem1->instanceType() == Class_Cylinder)
	{
	  double len1 = fabs((loc1-loc2)*dir1);
	  if (len1 > tol)
	    continue;   // Not the same surfac

	  Cylinder *cyl1 = dynamic_cast<Cylinder*>(elem1);
	  Cylinder *cyl2 = dynamic_cast<Cylinder*>(elem2);
	  if (cyl1 == NULL || cyl2 == NULL)
	    return surf;  // Something is wrong

	  double rad1 = cyl1->getRadius();
	  double rad2 = cyl2->getRadius();
	  if (fabs(rad2 - rad1) > tol)
	    continue;   // Not the same surface

	  Point axis1_1, axis1_2, axis1_3, axis2_1, axis2_2, axis2_3;;
	  cyl1->getCoordinateAxes(axis1_1, axis1_2, axis1_3);
	  cyl2->getCoordinateAxes(axis2_1, axis2_2, axis2_3);
	  bool swapped = ((cyl1->isSwapped() && !cyl2->isSwapped()) || 
			  (!cyl1->isSwapped() && cyl2->isSwapped()));
	  double len2 = loc1.dist(loc2);
	  double ang2 = axis1_1.angle(axis2_1);
	  if (len2 > tol || ang2 > angtol)
	    {
	      // Parameterization of the two cylinders differ.
	      // Modify parameter domain of surface 2 to match that of surface 1
	      Vector2D low = dom2.lowerLeft();
	      Vector2D high = dom2.upperRight();
	      if (len2 > tol)
		{
		  // Move domain
		  Point diff = loc2 - loc1;
		  int sgn = (diff * dir1 > 0.0);
		  low[1] += sgn*len2;
		  high[1] += sgn*len2;
		}
	      if (ang2 > angtol)
		{
		  Point v = axis1_1.cross(axis2_1);
		  std::cout << "Reparameterization of cylinder. To be continued" << std::endl;
		}
	      if (swapped)
		{
		  std::swap(low[0], low[1]);
		  std::swap(high[0], high[1]);
		}
	      RectDomain dom3(low, high);
	      dom1.addUnionWith(dom3);
	    }
	  else
	    {
	      if (swapped)
		{
		  Vector2D low = dom2.lowerLeft();
		  Vector2D high = dom2.upperRight();
		  std::swap(low[0], low[1]);
		  std::swap(high[0], high[1]);
		  RectDomain dom3(low, high);
		  dom1.addUnionWith(dom3);
		}
	      else
		dom1.addUnionWith(dom2);
	    }
	}
      else if (elem1->instanceType() == Class_Cone)
	{
	  Cone *cone1 = dynamic_cast<Cone*>(elem1);
	  Cone *cone2 = dynamic_cast<Cone*>(elem2);
	  if (cone1 == NULL || cone2 == NULL)
	    return surf;  // Something is wrong
	  double angle1 = cone1->getConeAngle();
	  double angle2 = cone2->getConeAngle();
	  double ang = dir1.angle(dir2);
	  if (fabs(angle1-angle2) > angtol || 
	      (ang > angtol && M_PI-ang > angtol))
	    continue; // Not the same surface

	  Point axis1_1, axis1_2, axis1_3, axis2_1, axis2_2, axis2_3;;
	  cone1->getCoordinateAxes(axis1_1, axis1_2, axis1_3);
	  cone2->getCoordinateAxes(axis2_1, axis2_2, axis2_3);
	  bool swapped = ((cone1->isSwapped() && !cone2->isSwapped()) || 
			  (!cone1->isSwapped() && cone2->isSwapped()));
	  double len2 = loc1.dist(loc2);
	  double ang2 = axis1_1.angle(axis2_1);
	  if (len2 > tol || ang2 > angtol)
	    {
	      // The parameterization differ
	      // Modify parameter domain of surface 2 to match that of surface 1
	      Vector2D low = dom2.lowerLeft();
	      Vector2D high = dom2.upperRight();
	      if (len2 > tol)
		{
		  // Move domain
		  Point diff = loc2 - loc1;
		  int sgn = (diff * dir1 > 0.0);
		  low[1] += sgn*len2;
		  high[1] += sgn*len2;
		}
	      if (ang2 > angtol)
		{
		  std::cout << "Reparameterization of cone. To be continued" << std::endl;

		}
	      if (swapped)
		{
		  std::swap(low[0], low[1]);
		  std::swap(high[0], high[1]);
		}
	      RectDomain dom3(low, high);
	      dom1.addUnionWith(dom3);
	    }
	  else
	    {
	      if (swapped)
		{
		  Vector2D low = dom2.lowerLeft();
		  Vector2D high = dom2.upperRight();
		  std::swap(low[0], low[1]);
		  std::swap(high[0], high[1]);
		  RectDomain dom3(low, high);
		  dom1.addUnionWith(dom3);
		}
	      else
		dom1.addUnionWith(dom2);
	    }
	}
      else if (elem1->instanceType() == Class_Sphere)
	{
	  std::cout << "Sphere" << std::endl;
	}
      else if (elem1->instanceType() == Class_Torus)
	{
	  std::cout << "Torus" << std::endl;
	}
      else
	return surf;  // Surface type not supported
    }

  // Create surface. Case distinction
  if (elem1->instanceType() == Class_Plane)
    {
      Plane *plane1 = dynamic_cast<Plane*>(elem1);
      Point axis1, axis2;
      plane1->getSpanningVectors(axis1, axis2);
      shared_ptr<Plane> plane2(new Plane(loc1, dir1, axis1, plane1->isSwapped()));
      plane2->setParameterBounds(dom1.umin(), dom1.vmin(), dom1.umax(), dom1.vmax());
      surf = plane2;
    }
  else if (elem1->instanceType() == Class_Cylinder)
    {
      Cylinder *cyl1 = dynamic_cast<Cylinder*>(elem1);
      double rad1 = cyl1->getRadius();
      Point axis1, axis2, axis3;
      cyl1->getCoordinateAxes(axis1, axis2, axis3);
      shared_ptr<Cylinder> cyl2(new Cylinder(rad1, loc1, dir1, axis1, 
					     cyl1->isSwapped()));
      cyl2->setParameterBounds(dom1.umin(), dom1.vmin(), dom1.umax(), dom1.vmax());
      surf = cyl2;
    }
  else if (elem1->instanceType() == Class_Cone)
    {
      Cone *cone1 = dynamic_cast<Cone*>(elem1);
      double rad = cone1->getRadius();
      double cone_angle = cone1->getConeAngle();
      Point axis1, axis2, axis3;
      cone1->getCoordinateAxes(axis1, axis2, axis3);
      shared_ptr<Cone> cone2(new Cone(rad, loc1, dir1, axis1, 
				      cone_angle, cone1->isSwapped()));
      cone2->setParameterBounds(dom1.umin(), dom1.vmin(), dom1.umax(), dom1.vmax());

      surf = cone2;
    }

#ifdef DEBUG
  std::ofstream of("elem_faces.g2");
  for (size_t ka=0; ka<sf_set.size(); ++ka)
    {
      sf_set[ka]->surface()->writeStandardHeader(of);
      sf_set[ka]->surface()->write(of);
    }
  surf->writeStandardHeader(of);
  surf->write(of);
#endif
  return surf;
}
