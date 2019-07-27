#ifndef INCLUDE_CURVEDRESLICE_HPP_
#define INCLUDE_CURVEDRESLICE_HPP_

/*
    Once we have a single lung we can use a new curved linear coordinate system
    and re-slice the lung. The goal is to provide a convenient planar view of both
    lungs in a single slice.

    We are using thin-plate splines with 4x4 control points that are adjusted
    to fit the lungs curved appearance in a quasy-sagittal view. Once we find the
    cuved best approximation we get the image stack by translating that surface.
    The direction of translation is not yet fully clear.

    It is not clear if the support by data of the spline in the top part is sufficient
    as the lung is basically a triangular shape. Should we remove the control points
    at the top? We could do a 2, 2, 3, 3 design initially and adjust that?
 */

#include "../mytypes.h" // get the types required to access the image data
#include "itkImageRegionIterator.h"
#include "linalg3d.h"
#include "ludecomposition.h"
#include <boost/numeric/ublas/matrix.hpp>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

using namespace boost::numeric::ublas;

#define GRID_W 4
#define GRID_H 4
static float grid[GRID_W][GRID_H];

std::vector<Vec> control_points;
int selected_cp = -1;

double regularization = 0.0;
double bending_energy = 0.0;

static double tps_base_func(double r) {
  if (r == 0.0)
    return 0.0;
  else
    return r * r * log(r);
}

/*
 *  Calculate Thin Plate Spline (TPS) weights from
 *  control points and build a new height grid by
 *  interpolating with them.
 */
static void calc_tps() {
  // You We need at least 3 points to define a plane
  if (control_points.size() < 3)
    return;

  unsigned p = control_points.size();

  // Allocate the matrix and vector
  matrix<double> mtx_l(p + 3, p + 3);
  matrix<double> mtx_v(p + 3, 1);
  matrix<double> mtx_orig_k(p, p);

  // Fill K (p x p, upper left of L) and calculate
  // mean edge length from control points
  //
  // K is symmetrical so we really have to
  // calculate only about half of the coefficients.
  double a = 0.0;
  for (unsigned i = 0; i < p; ++i) {
    for (unsigned j = i + 1; j < p; ++j) {
      Vec pt_i = control_points[i];
      Vec pt_j = control_points[j];
      pt_i.y = pt_j.y = 0;
      double elen = (pt_i - pt_j).len();
      mtx_l(i, j) = mtx_l(j, i) = mtx_orig_k(i, j) = mtx_orig_k(j, i) = tps_base_func(elen);
      a += elen * 2; // same for upper & lower tri
    }
  }
  a /= (double)(p * p);

  // Fill the rest of L
  for (unsigned i = 0; i < p; ++i) {
    // diagonal: reqularization parameters (lambda * a^2)
    mtx_l(i, i) = mtx_orig_k(i, i) = regularization * (a * a);

    // P (p x 3, upper right)
    mtx_l(i, p + 0) = 1.0;
    mtx_l(i, p + 1) = control_points[i].x;
    mtx_l(i, p + 2) = control_points[i].z;

    // P transposed (3 x p, bottom left)
    mtx_l(p + 0, i) = 1.0;
    mtx_l(p + 1, i) = control_points[i].x;
    mtx_l(p + 2, i) = control_points[i].z;
  }
  // O (3 x 3, lower right)
  for (unsigned i = p; i < p + 3; ++i)
    for (unsigned j = p; j < p + 3; ++j)
      mtx_l(i, j) = 0.0;

  // Fill the right hand vector V
  for (unsigned i = 0; i < p; ++i)
    mtx_v(i, 0) = control_points[i].y;
  mtx_v(p + 0, 0) = mtx_v(p + 1, 0) = mtx_v(p + 2, 0) = 0.0;

  // Solve the linear system "inplace"
  if (0 != LU_Solve(mtx_l, mtx_v)) {
    puts("Singular matrix! Aborting.");
    exit(1);
  }

  // Interpolate grid heights
  for (int x = -GRID_W / 2; x < GRID_W / 2; ++x) {
    for (int z = -GRID_H / 2; z < GRID_H / 2; ++z) {
      double h = mtx_v(p + 0, 0) + mtx_v(p + 1, 0) * x + mtx_v(p + 2, 0) * z;
      Vec pt_i, pt_cur(x, 0, z);
      for (unsigned i = 0; i < p; ++i) {
        pt_i = control_points[i];
        pt_i.y = 0;
        h += mtx_v(i, 0) * tps_base_func((pt_i - pt_cur).len());
      }
      grid[x + GRID_W / 2][z + GRID_H / 2] = h;
    }
  }

  // Calc bending energy
  matrix<double> w(p, 1);
  for (int i = 0; i < p; ++i)
    w(i, 0) = mtx_v(i, 0);
  matrix<double> be = prod(prod<matrix<double>>(trans(w), mtx_orig_k), w);
  bending_energy = be(0, 0);
}

// should return a fitted thin-plate spline to the lung data for each lung area
ImageType::Pointer computeReslice(ImageType::Pointer final, ImageType::Pointer labelField,
                                  int labelFromLabelField) {
  // we need to place the initial spline with its control points in roughly the right orientation
  // that would be a sagittal orientation ontop of the bounding box of this lung object
  /////////////////////////////////////////////////
  // compute bounding box in physical coordinates
  /////////////////////////////////////////////////
  int bb[6]; // xmin, ymin, zmin, xmax, ymax, zmax in physical coordinates - to make the sampling on
             // the final spline surface easier - geometrically correct
  bb[0] = bb[3] = final->GetOrigin()[0];
  bb[1] = bb[4] = final->GetOrigin()[1];
  bb[2] = bb[5] = final->GetOrigin()[2];
  ImageType::IndexType pixelIndex;
  using PointType = itk::Point<double, ImageType::ImageDimension>;
  PointType point;
  // iterate through all pixel of this label
  ImageType::RegionType lungRegion = labelField->GetLargestPossibleRegion();
  itk::ImageRegionIterator<ImageType> lungIterator(labelField, lungRegion);
  while (!lungIterator.IsAtEnd()) {
    if (lungIterator.Get() == labelFromLabelField) {
      pixelIndex = lungIterator.GetIndex();
      final->TransformIndexToPhysicalPoint(pixelIndex, point);
      if (point[0] < bb[0])
        bb[0] = point[0];
      if (point[1] < bb[1])
        bb[1] = point[1];
      if (point[2] < bb[2])
        bb[2] = point[2];
      if (point[0] > bb[3])
        bb[3] = point[0];
      if (point[1] > bb[4])
        bb[4] = point[1];
      if (point[2] > bb[5])
        bb[5] = point[2];
    }
  }

  //////////////////////////////////////////////////////////////////////////////////
  // we also need a sample of the voxels for this lung, like the once at the border
  //////////////////////////////////////////////////////////////////////////////////
  std::vector<Vec> samples;

  return ImageType::Pointer();
}

#endif // INCLUDE_CURVEDRESLICE_HPP_