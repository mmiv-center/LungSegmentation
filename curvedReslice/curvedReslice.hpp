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
#include "itkNearestNeighborInterpolateImageFunction.h"
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
    reset(cp);
  }

  void reset(std::vector<Vec> cp) {
    this->control_points.resize(0);
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
    int p = this->control_points.size();
    double h = this->mtx_v(p + 0, 0) + this->mtx_v(p + 1, 0) * x + this->mtx_v(p + 2, 0) * z;
    Vec pt_i, pt_cur(x, 0, z);
    for (unsigned i = 0; i < p; ++i) {
      pt_i = this->control_points[i];
      pt_i.y = 0;
      h += this->mtx_v(i, 0) * tps_base_func((pt_i - pt_cur).len());
    }
    return h;
  }

  /*
   *  Calculate Thin Plate Spline (TPS) weights from
   *  control points.
   */
  void calc_tps() {
    // You We need at least 3 points to define a plane
    if (this->control_points.size() < 3)
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
      this->mtx_v(i, 0) = this->control_points[i].y;
    this->mtx_v(p + 0, 0) = this->mtx_v(p + 1, 0) = this->mtx_v(p + 2, 0) = 0.0;

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

#define V_MINUS(A, B)                                                                                                                                          \
  { A[0] - B[0], A[1] - B[1], A[2] - B[2] }
#define V_CROSS(A, B)                                                                                                                                          \
  { A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0] }
#define V_NORMALIZE(A)                                                                                                                                         \
  {                                                                                                                                                            \
    A[0] / sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]), A[1] / sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]),                                                \
        A[2] / sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2])                                                                                                   \
  }

// should return a fitted thin-plate spline to the lung data for each lung area
ImageType::Pointer computeReslice(ImageType::Pointer final, ImageType::Pointer labelField, int labelFromLabelField, bool verbose) {
  // we need to place the initial spline with its control points in roughly the right orientation
  // that would be a sagittal orientation ontop of the bounding box of this lung object
  /////////////////////////////////////////////////
  // compute bounding box in physical coordinates
  /////////////////////////////////////////////////
  if (verbose)
    fprintf(stdout, "compute bounding box in physical coordinates...");
  float bb[6]; // xmin, ymin, zmin, xmax, ymax, zmax in physical coordinates - to make the sampling on
               // the final spline surface easier - geometrically correct
  ImageType::IndexType pixelIndex;
  using PointType = itk::Point<double, ImageType::ImageDimension>;
  PointType point;
  // iterate through all pixel of this label
  ImageType::RegionType lungRegion = labelField->GetLargestPossibleRegion();
  ImageType::RegionType finalRegion = final->GetLargestPossibleRegion();
  itk::ImageRegionIterator<ImageType> lungIterator(labelField, lungRegion);
  itk::ImageRegionIterator<ImageType> finalIterator(final, finalRegion);
  size_t counter = 0;
  while (!lungIterator.IsAtEnd() && !finalIterator.IsAtEnd()) {
    if (counter == 0 && lungIterator.Get() == labelFromLabelField) { // init the bounding box with the first point
      pixelIndex = finalIterator.GetIndex();
      final->TransformIndexToPhysicalPoint(pixelIndex, point);
      bb[0] = bb[3] = point[0];
      bb[1] = bb[4] = point[1];
      bb[2] = bb[5] = point[2];
      counter++; // only do this once
    }
    if (lungIterator.Get() == labelFromLabelField) {
      // finalIterator.Set(2000); // debug set the current label in final to single value
      pixelIndex = finalIterator.GetIndex();
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
    ++finalIterator;
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
  // In the tps the dimension we change is y with x and z defining the plane.
  // Here our control points are actually: h - x, tps - y, w - z
  ///////////////////////////////////////////////////////////////////////////////
  int CGRID_W = 5;
  int CGRID_H = 5;
  std::vector<Vec> control_points;
  for (int w = 0; w < CGRID_W; w++) {   // this is a z-direction?
    for (int h = 0; h < CGRID_H; h++) { // this is the y-direction?
      // what is the orientation of the sheet? If z is up/down we should use y for front-back?
      control_points.push_back(
          Vec(bb[1] + h * ((bb[4] - bb[1]) / (CGRID_H - 1)), bb[0] + (bb[3] - bb[0]) / 2.0, bb[2] + w * ((bb[5] - bb[2]) / (CGRID_W - 1))));
      fprintf(stdout, "control point: %f %f %f\n", control_points.back().x, control_points.back().y, control_points.back().z);
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  // We also need a grid (finite sampling) for the computation of the metric.
  // Strange thing here is that we will fix x and z. Only y will be provided by the
  // thin-plate spline.
  //////////////////////////////////////////////////////////////////////////////
  if (verbose)
    fprintf(stdout, "compute grid position...");
  int GRID_W = 60;
  int GRID_H = 60;
  // store the coordinates (x,z) but calculate the height value y later with interpolate_height from tps
  matrix<Vec> grid_pos(GRID_W, GRID_H);
  for (int w = 0; w < GRID_W; w++) {
    for (int h = 0; h < GRID_H; h++) {
      grid_pos(w, h).x = bb[1] + h * ((bb[4] - bb[1]) / (GRID_H - 1)); // initially we have a flat sheet here - hope this will curve
      grid_pos(w, h).y = bb[0] + ((bb[3] - bb[0]) / 2.0);              // flat value will be calculated!
      grid_pos(w, h).z = bb[2] + w * ((bb[5] - bb[2]) / (GRID_W - 1));
      // fprintf(stdout, "grid point: %f %f %f\n", grid_pos(w, h).x, grid_pos(w, h).y, grid_pos(w, h).z);
    }
  }

  // we simulate the tpspline in y direction (x,h,z) but we draw the tps in x direction (left-right, (h, x, z))
  // so sampling the data should also use the drawing direction - to solve our problem

  // this makes a copy of the control_points
  if (verbose)
    fprintf(stdout, "create tps...\n");

  tpspline *tps = new tpspline(control_points);
  // we don't use the spacing anymore, only height adjustments are allowed
  double expected_spacing[2] = {(bb[4] - bb[1]) / (CGRID_H - 1.0), (bb[5] - bb[2]) / (CGRID_W - 1.0)};
  tps->regularization = 0.025;
  tps->calc_tps();

  { // save a version of the data (interpolated we expect that all values are inside the volume - in physical coordinates)
    fprintf(stdout, "bounding box is here: %f %f %f %f %f %f\n", bb[0], bb[1], bb[2], bb[3], bb[4], bb[5]);
    fprintf(stdout, "are we at the right location?\n");
    PointType point;
    ImageType::IndexType pixelIndex;
    for (int w = 0; w < GRID_W; w++) {   // this is a z-direction?
      for (int h = 0; h < GRID_H; h++) { // this is the y-direction?
        point[1] = grid_pos(w, h).x;     // tps->interpolate_height(w, h); // grid_pos(w, h).y;
        point[2] = grid_pos(w, h).z;
        point[0] = tps->interpolate_height(point[1], point[2]); // grid_pos(w, h).y;
        bool ok = final->TransformPhysicalPointToIndex(point, pixelIndex);
        // now set this point in the data to some value
        if (ok) // is inside
          final->SetPixel(pixelIndex, 3000);
      }
    }

    if (1) { // debug output
      fprintf(stdout, "and save to temp file...\n");
      typedef itk::ImageFileWriter<ImageType> WriterType;
      WriterType::Pointer writer = WriterType::New();

      // std::string a(labelfieldfilename + "trachea.nii");
      writer->SetFileName("/tmp/debug_begin_tps.nii");
      writer->SetInput(final);

      std::cout << "Writing the final mask as " << std::endl;
      std::cout << "/tmp/debug_begin_tps.nii" << std::endl;

      try {
        writer->Update();
      } catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
        return ImageType::Pointer();
      }
    }
  }
  // exit(-1);

  if (verbose)
    fprintf(stdout, "control points: %lu, regularization: %2.3f, bending energy: %4.3f\n", control_points.size(), tps->regularization, tps->bending_energy);

  // http://numerical.recipes/forumarchive/index.php/t-1360.html
  const int NDIM = control_points.size(); // only height is adjusted
  const Doub GTOL = 1.0e-4;
  int iter;
  Doub fret;
  VecDoub_IO p(NDIM);

  struct Ftor {
    double expected_spacing[2];
    matrix<Vec> *grid_pos;
    std::vector<Vec> *orig_control_points;
    tpspline *tps; // needs to be set before calling the operator()
    int CGRID_W;
    int CGRID_H;
    // get access to the x and z position (they are not adjusted by the algorithm)
    float *bb;
    int labelFromLabelField;
    // we  need access to the image volume as well - for sampling of intensities
    ImageType::Pointer labelField;

    Doub operator()(VecDoub_I &x) const {
      // lets copy the values back to get a control_points array
      std::vector<Vec> control_points(x.size()); // only the heights in here (y) are adjusted
      for (int i = 0; i < control_points.size(); i++) {
        control_points[i].x = (*orig_control_points)[i].x; // this should be h
        control_points[i].y = x[i];                        // this should be the interpolated height (as y value for tps
        control_points[i].z = (*orig_control_points)[i].z; // this should be w
      }
      // fprintf(stdout, "y[0] %f\n", x[0]); // should change
      tps->reset(control_points);
      tps->calc_tps();

      // now we need a new grid_pos with values!
      // but out grid_pos values only change in y, and we can calc those using by hand (so we don't update the array)

      // FIRST PART: calculate the current spacing - should be expected_spacing!
      // fprintf(stdout, "calculate the sum_spacing...\n");
      int count_spacings = 0;
      double sum_spacing = 0.0;
      // we don't need to calculate the spacing anymore... we fix the x, z location of the grid and only allow the y direction to change
      /*       for (int w = 1; w < CGRID_W; w++) {   // this is a z-direction?
              for (int h = 1; h < CGRID_H; h++) { // this is the y-direction?
                // y is always regularized by bending energy and can be ignored here
                int i = w * CGRID_H + h;
                double ds = sqrt((control_points[i].x - control_points[i - 1].x) * (control_points[i].x - control_points[i - 1].x) +
                                 // (control_points[i].y - control_points[i - 1].y) * (control_points[i].y - control_points[i - 1].y) +
                                 (control_points[i].z - control_points[i - 1].z) * (control_points[i].z - control_points[i - 1].z));
                // fprintf(stdout, "current distance to add is: %f sum is: %f expected spacing: %f %f\n", ds, sum_spacing, expected_spacing[0],
         expected_spacing[1]); if (std::isfinite(ds)) { sum_spacing += (ds - expected_spacing[0]) * (ds - expected_spacing[0]); count_spacings++;
                }
                int j = (w - 1) * CGRID_H + h;
                ds = sqrt((control_points[i].x - control_points[j].x) * (control_points[i].x - control_points[j].x) +
                          //(control_points[i].y - control_points[j].y) * (control_points[i].y - control_points[j].y) +
                          (control_points[i].z - control_points[j].z) * (control_points[i].z - control_points[j].z));
                if (std::isfinite(ds)) {
                  sum_spacing += (ds - expected_spacing[1]) * (ds - expected_spacing[1]);
                  count_spacings++;
                }
                // fprintf(stdout, "current distance to add is: %f %f: divided %f\n", ds, sum_spacing, sum_spacing / count_spacings);
              }
            } */
      // fprintf(stdout, "done...\n");
      // THIRD PART: calculate for a given spacing the normal of our thin-plate spline
      // we can look at all the grid points now and calculate the normals
      int GRID_W = 60; // the resolution at which we compute the normals
      int GRID_H = 60;
      PointType point;
      ImageType::IndexType pixelIndex;

      // fprintf(stdout, "compute normal distances...\n");
      // we need to get voxel values from a position
      double sum_dist = 0.0;               // should be 0 in a perfect world - if location of spline is in each geometric mean of normal trace
      int num_sum_dist = 0;                // could now many samples we do and use the average to make the number smaller
      for (int w = 1; w < GRID_W; w++) {   // this is a z-direction?
        for (int h = 1; h < GRID_H; h++) { // this is the y-direction?
          // compute the normal at this point
          double ax = (*grid_pos)(w - 0, h - 0).x;
          double az = (*grid_pos)(w - 0, h - 0).z;
          double a[] = {ax, tps->interpolate_height(ax, az), az};
          ax = (*grid_pos)(w - 1, h - 0).x;
          az = (*grid_pos)(w - 1, h - 0).z;
          double b[] = {ax, tps->interpolate_height(ax, az), az};
          ax = (*grid_pos)(w - 0, h - 1).x;
          az = (*grid_pos)(w - 0, h - 1).z;
          double c[] = {ax, tps->interpolate_height(ax, az), az};
          double ab[] = V_MINUS(a, b);
          double cb[] = V_MINUS(c, b);
          double nn[] = V_CROSS(cb, ab); // this is the normal
                                         // now sample in the +- normal direction starting from a and report back the center of mass of the
                                         // foreground object is n of length 1?
          double n[] = V_NORMALIZE(nn);  // normal on a[]
          // fprintf(stdout, "normal here is : %f %f %f for a: %f %f %f b: %f %f %f c: %f %f %f ab: %f %f %f cb: %f %f %f\n", n[0], n[1], n[2], a[0], a[1],
          // a[2],
          //        b[0], b[1], b[2], c[0], c[1], c[2], ab[0], ab[1], ab[2], cb[0], cb[1], cb[2]);

          // problem with the approach below is that we don't get a propper gradient with this approach
          // only every pixel shift will we be getting a different weight, not inside the voxel itself,
          // this make it impossible for the linsearch to find the initial gradient, better to get a
          // true sub-pixel distance here

          // sample the line - this should be done by linmin instead, would make this faster
          double sum_weights = 0.0;
          double sum_weights_pos[3] = {0, 0, 0}; // calculate the center of mass in the space direction
          int num_samples = 0;                   // how many samples we use in the calculation (all inside the region of interest)
          for (float weight = 0.0; weight < (bb[3] - bb[0]); weight += 1.0) {
            // here we want to sample (h, x, z) - image space instead of (x, h, z) - tps-space
            point[0] = a[1] + weight * n[1];
            point[1] = a[0] + weight * n[0];
            point[2] = a[2] + weight * n[2];

            bool ok = labelField->TransformPhysicalPointToIndex(point, pixelIndex);
            if (!ok) // if not inside stop here
              break;
            unsigned int val = labelField->GetPixel(pixelIndex);
            if (val == labelFromLabelField) { // we have a voxel inside the label
              // dist to current point in double resolution, important for the gradient calculation
              PointType point2;
              labelField->TransformIndexToPhysicalPoint(pixelIndex, point2);
              sum_weights_pos[0] += point2[0];
              sum_weights_pos[1] += point2[1];
              sum_weights_pos[2] += point2[2];
              // double square_distance_to_center =
              //    (point2[0] - a[1]) * (point2[0] - a[1]) + (point2[1] - a[0]) * (point2[1] - a[0]) + (point2[2] - a[2]) * (point2[2] - a[2]);

              // sum_weights += square_distance_to_center; // add the weight works if n has length 1
              num_samples++;
            }
          }
          for (float weight = 0.0; weight < (bb[3] - bb[0]);
               weight += 1.0) { // instead of 1.0 this could be smallest voxel size, or the voxel size in the space dimension (linear interpolation between
                                // voxel_size[x] based on dot-product)
            point[0] = a[1] - weight * n[1];
            point[1] = a[0] - weight * n[0];
            point[2] = a[2] - weight * n[2];

            bool ok = labelField->TransformPhysicalPointToIndex(point, pixelIndex);
            if (!ok) // outside the volume, don't sample there
              break;
            unsigned int val = labelField->GetPixel(pixelIndex);
            if (val == labelFromLabelField) {
              PointType point2;
              labelField->TransformIndexToPhysicalPoint(pixelIndex, point2);
              sum_weights_pos[0] += point2[0];
              sum_weights_pos[1] += point2[1];
              sum_weights_pos[2] += point2[2];
              // double square_distance_to_center =
              //    (point2[0] - a[1]) * (point2[0] - a[1]) + (point2[1] - a[0]) * (point2[1] - a[0]) + (point2[2] - a[2]) * (point2[2] - a[2]);

              // sum_weights -= square_distance_to_center;
              num_samples++;
            }
          }
          // ok we got now d1 and d2 as the distances at which we are outside of the slice, but what if we are outside
          // to start with?, better to use the distance to the geometric mean on the line
          if (num_samples == 0) {
            // no point on this line, ignore this point
            // fprintf(stdout, "no points found.. ignore\n");
            continue; // do the next sample point
          }
          // center of mass along this spacial direction (n)
          sum_weights_pos[0] /= num_samples;
          sum_weights_pos[1] /= num_samples;
          sum_weights_pos[2] /= num_samples;
          double square_distance_to_center = (sum_weights_pos[0] - a[1]) * (sum_weights_pos[0] - a[1]) +
                                             (sum_weights_pos[1] - a[0]) * (sum_weights_pos[1] - a[0]) +
                                             (sum_weights_pos[2] - a[2]) * (sum_weights_pos[2] - a[2]);
          sum_dist += square_distance_to_center;
          // sum_dist += (sum_weights / num_samples); // in the direction of n at a, could be 0 when samples balance in +-
          num_sum_dist++;
        }
      }
      // fprintf(stdout, "add a sum: %f: %d...\n", sum_dist, num_sum_dist);
      // we are only interested in the absolute distance to the center
      if (sum_dist < 0)
        sum_dist *= -1.0;

      Doub metric = 0.0;
      // what is the scale of these three parts? They should have equal size to make the metric use all three to optimize
      // FIRST PART
      double t = 0.6; // weighting between data fit and bending energy
      metric += (sum_dist / num_sum_dist) * (1.0 - t);
      // SECOND PART
      metric += t * tps->bending_energy;
      // THIRD PART
      // metric += sum_spacing / count_spacings;
      fprintf(stdout, "metric: center-line = %f, bending energy = %f\n", (sum_dist / num_sum_dist) * (1.0 - t), (t * tps->bending_energy));

      return metric;
    }
  };

  Ftor f;
  f.tps = tps;
  f.CGRID_W = CGRID_W;
  f.CGRID_H = CGRID_H;
  f.grid_pos = &grid_pos;
  f.labelField = labelField;
  f.bb = bb;
  f.orig_control_points = &control_points; // keep these fixed (the coordinates that are not adjusted)
  f.labelFromLabelField = labelFromLabelField;
  f.expected_spacing[0] = expected_spacing[0];
  f.expected_spacing[1] = expected_spacing[1];

  Funcd<Ftor> aa(f);
  // we can only change the height of control points
  for (int i = 0; i < control_points.size(); i++) {
    p[i] = control_points[i].y;
  }
  fprintf(stdout, "start dfpmin with %f %f %f %f...\n", p[0], p[1], p[2], p[3]);
  dfpmin(p, GTOL, iter, fret, aa);
  // copy the updated control points back to the spline - we can use them later to sample the volume
  for (int i = 0; i < control_points.size(); i++) {
    //    control_points[i].x = p[i * 3 + 0];
    control_points[i].y = p[i];
    //    control_points[i].z = p[i * 3 + 2];
  }
  fprintf(stdout, "call calc_tps again...\n");
  tps->reset(control_points);
  tps->calc_tps();
  // now the final spline can be used to resample the data
  fprintf(stdout, "set pixel in final to some large value...\n");
  for (int w = 0; w < GRID_W - 1; w++) {   // this is a z-direction?
    for (int h = 0; h < GRID_H - 1; h++) { // this is the y-direction?
      PointType point;
      ImageType::IndexType pixelIndex;
      point[1] = grid_pos(w, h).x;
      point[2] = grid_pos(w, h).z;
      point[0] = tps->interpolate_height(point[1], point[2]);
      bool ok = final->TransformPhysicalPointToIndex(point, pixelIndex);
      // now set this point in the data to some value
      if (ok)
        final->SetPixel(pixelIndex, 3000);
    }
  }

  if (1) { // debug output
    fprintf(stdout, "and save to temp file...\n");
    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();

    // std::string a(labelfieldfilename + "trachea.nii");
    writer->SetFileName("/tmp/debug.nii");
    writer->SetInput(final);

    std::cout << "Writing the final mask as " << std::endl;
    std::cout << "/tmp/debug.nii" << std::endl;

    try {
      writer->Update();
    } catch (itk::ExceptionObject &ex) {
      std::cout << ex << std::endl;
      return ImageType::Pointer();
    }
  }

  // ok, now we are ready to create a new Image object that will contain our sampled image information
  ImageType::Pointer resliced = ImageType::New();
  ImageType::RegionType::SizeType size;
  size[0] = 1024;
  size[1] = 512;
  size[2] = 512;
  ImageType::RegionType::IndexType index; // lower left corner 1,1,1
  index.Fill(0);
  ImageType::RegionType region(index, size);
  // we assume that our image data is centered on one side of the image
  // because we want to combine left and right lung in one image
  resliced->SetRegions(region);
  resliced->Allocate();
  // the image does not have a special origin
  ImageType::PointType newOrigin;
  newOrigin.Fill(0.0);
  resliced->SetOrigin(newOrigin);

  ImageType::DirectionType direction;
  direction.SetIdentity();
  resliced->SetDirection(direction);

  ImageType::SpacingType spacing;
  // Units (e.g., mm, inches, etc.) are defined by the application.
  spacing[0] = 1; // spacing along X
  spacing[1] = 1; // spacing along Y
  spacing[2] = 1; // spacing along Z
  resliced->SetSpacing(spacing);

  // ok we can sample the volume now by translating one slice in x direction
  // (better) we start with the voxel in the output and go backwards - would result in better quality images
  // what is our offset?
  GRID_W = 512;
  GRID_H = 512;
  double padding = 5; // a 5 voxel padding on the border
  // store the coordinates (x,z) but calculate the height value y later with interpolate_height from tps
  matrix<Vec> grid_pos2(GRID_W, GRID_H);
  for (int w = 0; w < GRID_W; w++) {
    for (int h = 0; h < GRID_H; h++) {
      grid_pos2(w, h).x =
          (bb[1] - padding) + h * (((bb[4] + padding) - (bb[1] - padding)) / (GRID_H - 1)); // initially we have a flat sheet here - hope this will curve
      grid_pos2(w, h).z = (bb[2] - padding) + w * (((bb[5] + padding) - (bb[2] - padding)) / (GRID_W - 1));
      grid_pos2(w, h).y = tps->interpolate_height(grid_pos2(w, h).x, grid_pos2(w, h).z);
    }
  }

  for (int slice = 0; slice < size[2]; slice++) {
    for (int y = 0; y < size[1]; y++) {
      for (int x = 0; x < size[0]; x++) { // 1024
        int offset_x = (x - 512);         // assumes 1mm slice distance
        if (x < 512)
          continue;
        int xx = x - 512;
        PointType point;
        ImageType::IndexType pixelIndex;
        point[1] = grid_pos2(xx, y).x;
        point[2] = grid_pos2(xx, y).z;
        point[0] = grid_pos2(xx, y).y + (-255 + slice);
        // tps->interpolate_height(point[1], point[2]);
        ImageType::IndexType idx;
        idx[0] = x;
        idx[1] = slice;
        idx[2] = y;
        bool ok = final->TransformPhysicalPointToIndex(point, pixelIndex);
        if (ok) {
          double val = final->GetPixel(pixelIndex);
          // set this pixel into the current location
          resliced->SetPixel(idx, val);
        } else {
          resliced->SetPixel(idx, -1024);
        }
      }
    }
  }

  if (1) { // debug output
    fprintf(stdout, "and save resliced to temp file...\n");
    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();

    // std::string a(labelfieldfilename + "trachea.nii");
    writer->SetFileName("/tmp/debug_resliced.nii");
    writer->SetInput(resliced);

    std::cout << "Writing the resliced image as " << std::endl;
    std::cout << "/tmp/debug_resliced.nii" << std::endl;

    try {
      writer->Update();
    } catch (itk::ExceptionObject &ex) {
      std::cout << ex << std::endl;
      return ImageType::Pointer();
    }
  }

  return resliced;
}

#endif // INCLUDE_CURVEDRESLICE_HPP_