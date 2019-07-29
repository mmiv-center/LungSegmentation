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

// use numerical recipes for the minimization
#include "numRec/ludcmp.h"
#include "numRec/nr3.h"
#include "numRec/qrdcmp.h"
#include "numRec/quasinewton.h"
#include "numRec/roots_multidim.h"

using namespace boost::numeric::ublas;

// This helper class preforms a thin-plate spline interpolation as a height of a 2-d grid.
class tpspline {
public:
  std::vector<Vec> control_points;

  double regularization = 0.0;
  double bending_energy = 0.0;

private:
  matrix<double> mtx_l;
  matrix<double> mtx_v;
  matrix<double> mtx_orig_k;

public:
  // do not allow to use the default constructor (C++11 feature)
  tpspline() = delete;

  tpspline(std::vector<Vec> cp) { // constructor, allocate memory
    // make a copy of the input vector
    for (int i = 0; i < cp.size(); i++)
      this->control_points.push_back(Vec(cp[i].x, cp[i].y, cp[i].z));
    int p = this->control_points.size();
    // create sufficient storage for the matrices
    this->mtx_l = matrix<double>(p + 3, p + 3);
    this->mtx_v = matrix<double>(p + 3, 1);
    this->mtx_orig_k = matrix<double>(p, p);
    // this->calc_tps(); // lets do this manually instead
  }

  static double tps_base_func(double r) {
    if (r == 0.0)
      return 0.0;
    else
      return r * r * log(r);
  }

  double interpolate_height(double x, double z) {
    int p = control_points.size();
    double h = mtx_v(p + 0, 0) + mtx_v(p + 1, 0) * x + mtx_v(p + 2, 0) * z;
    Vec pt_i, pt_cur(x, 0, z);
    for (unsigned i = 0; i < p; ++i) {
      pt_i = control_points[i];
      pt_i.y = 0;
      h += mtx_v(i, 0) * tps_base_func((pt_i - pt_cur).len());
    }
    return h;
  }

  /*
   *  Calculate Thin Plate Spline (TPS) weights from
   *  control points.
   */
  void calc_tps() {
    // You We need at least 3 points to define a plane
    if (control_points.size() < 3)
      return;

    unsigned p = this->control_points.size();

    // Allocate the matrix and vector
    // matrix<double> mtx_l(p + 3, p + 3);
    // matrix<double> mtx_v(p + 3, 1);
    // matrix<double> mtx_orig_k(p, p);

    // Fill K (p x p, upper left of L) and calculate
    // mean edge length from control points
    //
    // K is symmetrical so we really have to
    // calculate only about half of the coefficients.
    double a = 0.0;
    for (unsigned i = 0; i < p; ++i) {
      for (unsigned j = i + 1; j < p; ++j) {
        Vec pt_i = this->control_points[i];
        Vec pt_j = this->control_points[j];
        pt_i.y = pt_j.y = 0;
        double elen = (pt_i - pt_j).len();
        this->mtx_l(i, j) = this->mtx_l(j, i) = this->mtx_orig_k(i, j) = this->mtx_orig_k(j, i) = tps_base_func(elen);
        a += elen * 2; // same for upper & lower tri
      }
    }
    a /= (double)(p * p);

    // Fill the rest of L
    for (unsigned i = 0; i < p; ++i) {
      // diagonal: reqularization parameters (lambda * a^2)
      mtx_l(i, i) = mtx_orig_k(i, i) = this->regularization * (a * a);

      // P (p x 3, upper right)
      mtx_l(i, p + 0) = 1.0;
      mtx_l(i, p + 1) = this->control_points[i].x;
      mtx_l(i, p + 2) = this->control_points[i].z;

      // P transposed (3 x p, bottom left)
      mtx_l(p + 0, i) = 1.0;
      mtx_l(p + 1, i) = this->control_points[i].x;
      mtx_l(p + 2, i) = this->control_points[i].z;
    }
    // O (3 x 3, lower right)
    for (unsigned i = p; i < p + 3; ++i)
      for (unsigned j = p; j < p + 3; ++j)
        this->mtx_l(i, j) = 0.0;

    // Fill the right hand vector V
    for (unsigned i = 0; i < p; ++i)
      this->mtx_v(i, 0) = control_points[i].y;
    this->mtx_v(p + 0, 0) = mtx_v(p + 1, 0) = mtx_v(p + 2, 0) = 0.0;

    // Solve the linear system "inplace"
    if (0 != LU_Solve(this->mtx_l, this->mtx_v)) {
      puts("Singular matrix! Aborting.");
      exit(1); // todo
    }

    // Interpolate grid heights
    /*     for (int x = -GRID_W / 2; x < GRID_W / 2; ++x) {
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
        } */

    // Calc bending energy
    matrix<double> w(p, 1);
    for (int i = 0; i < p; ++i)
      w(i, 0) = this->mtx_v(i, 0);
    matrix<double> be = prod(prod<matrix<double>>(trans(w), this->mtx_orig_k), w);
    this->bending_energy = be(0, 0);
  }
};

// ok these are the resolutions of our thin-plate spline - not the number of control points used to define it
/*
#define GRID_W 100
#define GRID_H 100
#define CGRID_W 4
#define CGRID_H 4
static float grid[GRID_W][GRID_H];

std::vector<Vec> control_points;

double regularization = 0.0;
double bending_energy = 0.0;

static double tps_base_func(double r) {
  if (r == 0.0)
    return 0.0;
  else
    return r * r * log(r);
}

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
*/

#define V_MINUS(A, B)                                                                                                                                          \
  { A[0] - B[0], A[1] - B[1], A[2] - B[2] }
#define V_CROSS(A, B)                                                                                                                                          \
  { A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0] }

// should return a fitted thin-plate spline to the lung data for each lung area
ImageType::Pointer computeReslice(ImageType::Pointer final, ImageType::Pointer labelField, int labelFromLabelField, bool verbose) {
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
    ++lungIterator;
  }

  //////////////////////////////////////////////////////////////////////////////////
  // How do we regularize the spline sheet? It could collapse, or it could fold over.
  // Fold over is easy - because we can use the bending energy - if we limit that
  // we should get a stiff plate (needs some testing). The second collapse issue
  // can be prevented if we keep a minimum distance between the control points.
  // So our energy term is composed of three parts: minimize the distance towards the border,
  // keep a low bending energy and keep the distance between each pair of control
  // points the same.
  //////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////
  // Ok, we have everything. Now we can arrange the spline control points initially
  ///////////////////////////////////////////////////////////////////////////////
  int CGRID_W = 4;
  int CGRID_H = 4;
  std::vector<Vec> control_points;
  for (int w = 0; w < CGRID_W; w++) {   // this is a z-direction?
    for (int h = 0; h < CGRID_H; h++) { // this is the y-direction?
      // what is the orientation of the sheet? If z is up/down we should use y for front-back?
      control_points.push_back(
          Vec(bb[0] + (bb[3] - bb[0]) / 2.0, bb[1] + h * ((bb[4] - bb[1]) / (CGRID_H - 1)), bb[2] + w * ((bb[5] - bb[2]) / (CGRID_W - 1))));
    }
  }
  // this makes a copy of the control_points
  tpspline *tps = new tpspline(control_points);

  double expected_spacing[2] = {(bb[4] - bb[1]) / (CGRID_H - 1.0), (bb[5] - bb[2]) / (CGRID_W - 1.0)};
  tps->regularization = 0.025;
  tps->calc_tps();

  if (verbose)
    fprintf(stdout, "control points: %lu, reqularization: %2.3f, bending energy: %4.3f\n", control_points.size(), tps->regularization, tps->bending_energy);

  // http://numerical.recipes/forumarchive/index.php/t-1360.html
  const int NDIM = control_points.size() * 3;
  const Doub GTOL = 1.0e-4;
  int iter;
  Doub fret;
  VecDoub_IO p(NDIM);

  struct Ftor {
    double expected_spacing[2];
    float **grid;
    tpspline *tps; // needs to be set before calling the operator()
    int CGRID_W;
    int CGRID_H;

    Doub operator()(VecDoub_I &x) const {
      // lets copy the values back to get a control_points array
      std::vector<Vec> control_points(x.size() / 3);
      for (int i = 0; i < control_points.size(); i++) {
        control_points[i].x = x[i * 3 + 0];
        control_points[i].y = x[i * 3 + 1];
        control_points[i].z = x[i * 3 + 2];
      }
      tps->calc_tps();

      // FIRST PART: calculate the current spacing - should be expected_spacing!
      int count_spacings = 0;
      double sum_spacing = 0.0;
      for (int w = 1; w < CGRID_W; w++) {   // this is a z-direction?
        for (int h = 1; h < CGRID_H; h++) { // this is the y-direction?
          // x is always regularized by bending energy and can be ignored here
          int i = w * CGRID_H + h;
          double ds = sqrt((control_points[i].x - control_points[i - 1].x) * (control_points[i].x - control_points[i - 1].x) +
                           (control_points[i].y - control_points[i - 1].y) * (control_points[i].y - control_points[i - 1].y) +
                           (control_points[i].z - control_points[i - 1].z) * (control_points[i].z - control_points[i - 1].z));
          sum_spacing = (ds - expected_spacing[0]) * (ds - expected_spacing[0]);
          count_spacings++;
          int j = (w - 1) * CGRID_H + h;
          ds = sqrt((control_points[i].x - control_points[j].x) * (control_points[i].x - control_points[j].x) +
                    (control_points[i].y - control_points[j].y) * (control_points[i].y - control_points[j].y) +
                    (control_points[i].z - control_points[j].z) * (control_points[i].z - control_points[j].z));
          sum_spacing = (ds - expected_spacing[1]) * (ds - expected_spacing[1]);
          count_spacings++;
        }
      }

      // THIRD PART: calculate for a given spacing the normal of our thin-plate spline
      // we can look at all the grid points now and calculate the normals
      int GRID_W = 100; // the resolution at which we compute the normals
      int GRID_H = 100;
      for (int w = 0; w < GRID_W - 1; w++) {   // this is a z-direction?
        for (int h = 0; h < GRID_H - 1; h++) { // this is the y-direction?
          // compute the normal at this point
          // float a[] = {grid[h][w].x, grid[h][w].y, grid[h][w].z};
          float b[] = {};
          float c[] = {};
          // float ab[] = V_MINUS(a, b);
          // float cb[] = V_MINUS(c, b);
          // float n[] = V_CROSS(cb, ab); // this is the normal
        }
      }
      Doub metric = 0.0;
      // SECOND PART
      metric += tps->bending_energy;
      metric += sum_spacing / count_spacings;

      return metric;
    }
  };

  Ftor f;
  f.tps = tps;
  f.CGRID_W = CGRID_W;
  f.CGRID_H = CGRID_H;
  f.expected_spacing[0] = expected_spacing[0];
  f.expected_spacing[1] = expected_spacing[1];

  Funcd<Ftor> aa(f);
  for (int i = 0; i < control_points.size(); i++) {
    p[i * 3 + 0] = control_points[i].x;
    p[i * 3 + 1] = control_points[i].y;
    p[i * 3 + 2] = control_points[i].z;
  }
  dfpmin(p, GTOL, iter, fret, aa);
  // copy the updated control points back to the spline - we can use them later to sample the volume
  for (int i = 0; i < control_points.size(); i++) {
    control_points[i].x = p[i * 3 + 0];
    control_points[i].y = p[i * 3 + 1];
    control_points[i].z = p[i * 3 + 2];
  }
  tps->calc_tps();
  // now the final spline can be used to resample the data

  return ImageType::Pointer();
}

#endif // INCLUDE_CURVEDRESLICE_HPP_