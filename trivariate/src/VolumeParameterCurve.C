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

#include "GoTools/trivariate/VolumeParameterCurve.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/Volumes.h"

using namespace Go;
using std::vector;

//===========================================================================
VolumeParameterCurve::VolumeParameterCurve(shared_ptr<ParamVolume> vol, 
					   shared_ptr<ParamCurve> crv)
  : vol_(vol), crv_(crv)
//===========================================================================
{
}

//===========================================================================
VolumeParameterCurve::VolumeParameterCurve(shared_ptr<ParamVolume> vol, 
					   shared_ptr<ParamCurve> crv,
					   shared_ptr<Point> par1,
					   shared_ptr<Point> par2)
  : vol_(vol), crv_(crv), par1_(par1), par2_(par2)
//===========================================================================
{
}

///===========================================================================
VolumeParameterCurve::~VolumeParameterCurve()
//===========================================================================
{
}
 
//===========================================================================
Point VolumeParameterCurve::eval(double t) const
//===========================================================================
{
  vector<Point> result(1);
  eval(t, 0, &result[0]);
  return result[0];
}

//===========================================================================
void VolumeParameterCurve::eval(double t, int n, Point der[]) const
//===========================================================================
{
  if (n > 1)
    n = 1;
  if (n < 0)
    n = 0;

  // Evaluate curve
  vector<Point> cv_der(n+1);
  crv_->point(cv_der, t, n);

  // Find closest point in volume
  double eps = 1.0e-6;
  double u1, v1, w1, dist;
  Point vol_pt;
  vol_->closestPoint(cv_der[0], u1, v1, w1, vol_pt, dist, eps);
   if (dist > 100.0*eps)
     std::cout << "Volume parameter curve, closest point: " << dist << std::endl;

  if (t == start() && par1_.get())
    {
      u1 = (*par1_)[0];
      v1 = (*par1_)[1];
      w1 = (*par1_)[2];
    }
  else if (t == end() && par2_.get())
    {
      u1 = (*par2_)[0];
      v1 = (*par2_)[1];
      w1 = (*par2_)[2];
    }

  der[0] = Point(u1, v1, w1);
  if (n == 1)
    {
      // Differentiate volume
      vector<Point> vol_der(4);
      vol_->point(vol_der, u1, v1, w1, 1);

      // Find the factors (r1, s1, t1) such that
      // r1*vol_der[1] + s1*vol_der[2] + t1*vol_der[3] = cv_der[1]
      // Solve by Cramers rule
      double det = vol_der[1][0]*(vol_der[2][1]*vol_der[3][2] -
				  vol_der[2][2]*vol_der[3][1]) -
	vol_der[1][1]*(vol_der[2][0]*vol_der[3][2] -
		       vol_der[2][2]*vol_der[3][0]) +
	vol_der[1][2]*(vol_der[2][0]*vol_der[3][1] -
		       vol_der[2][1]*vol_der[3][0]);
      double r1 = (cv_der[1][0]*(vol_der[2][1]*vol_der[3][2] -
				  vol_der[2][2]*vol_der[3][1]) -
		   cv_der[1][1]*(vol_der[2][0]*vol_der[3][2] -
				  vol_der[2][2]*vol_der[3][0]) +
		   cv_der[1][2]*(vol_der[2][0]*vol_der[3][1] -
				  vol_der[2][1]*vol_der[3][0]))/det;
      double s1 = (vol_der[1][0]*(cv_der[1][1]*vol_der[3][2] -
				  cv_der[1][2]*vol_der[3][1]) -
		   vol_der[1][1]*(cv_der[1][0]*vol_der[3][2] -
				  cv_der[1][2]*vol_der[3][0]) +
		   vol_der[1][2]*(cv_der[1][0]*vol_der[3][1] -
				  cv_der[1][1]*vol_der[3][0]))/det;
      double t1 = (vol_der[1][0]*(vol_der[2][1]*cv_der[1][2] -
				  vol_der[2][2]*cv_der[1][1]) -
		   vol_der[1][1]*(vol_der[2][0]*cv_der[1][2] -
				  vol_der[2][2]*cv_der[1][0]) +
		   vol_der[1][2]*(vol_der[2][0]*cv_der[1][1] -
				  vol_der[2][1]*cv_der[1][0]))/det;
//       Array<Point,3> a0(vol_der[1], vol_der[2], vol_der[3]);
//       Array<Point,3> a1(cv_der[1], vol_der[2], vol_der[3]);
//       Array<Point,3> a2(vol_der[1], cv_der[1], vol_der[3]);
//       Array<Point,3> a3(vol_der[1], vol_der[2], cv_der[1]);
//       double det = determinantOf(a0);
//       double r1 = determinantOf(a1)/det;
//       double s1 = determinantOf(a2)/det;
//       double t1 = determinantOf(a2)/det;

      // Dest
      Point tmp = r1*vol_der[1] + s1*vol_der[2] + t1*vol_der[3];
      dist = tmp.dist(cv_der[1]);
      if (dist > eps)
	std::cout << "Volume parameter curve, dist: " << dist << std::endl;

      der[1] = Point(r1, s1, t1);
    }
}

//===========================================================================
double VolumeParameterCurve::start() const
//===========================================================================
{
  return crv_->startparam();
}

//===========================================================================
double VolumeParameterCurve::end() const
//===========================================================================
{
  return crv_->endparam();
}

//===========================================================================
int VolumeParameterCurve::dim() const
//===========================================================================
{
// Trivariate parameter area
  return 3;
}

//===========================================================================
bool VolumeParameterCurve::approximationOK(double par, Point approxpos,
					   double tol1, double tol2) const
//===========================================================================
{
  Point pos = eval(par);
  double dist = pos.dist(approxpos);
  return (dist <= tol1);
}

