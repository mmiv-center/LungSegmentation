// set threads with
// ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=6 ./MaskRegistration ...

// ./MaskRegistration longData/049300006_049300006_seg/1.2.840.113654.2.3.1995.3.0.6.2017062709435100771.6151251220170622.nii
// longData/0493_00006_0493_00006_seg/1.2.840.113654.2.3.1995.3.0.6.2017053107400100715.101.051251220170523.nii longData/reg/volume.nii

// ./MaskRegistration -w -q 4 -m longData/049300006_049300006_seg/labels.nii -i longData/reg/difference.nii -b longData/reg/difference_before.nii -d
// longData/reg/deformation_field.nii longData/049300006_049300006_seg/1.2.840.113654.2.3.1995.3.0.6.2017062709435100771.6151251220170622.nii
// longData/0493_00006_0493_00006_seg/1.2.840.113654.2.3.1995.3.0.6.2017053107400100715.101.051251220170523.nii longData/reg/volume.nii

/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
// Software Guide : BeginLatex
//
// This example illustrates a realistic pipeline for solving a full deformable registration problem.
//
// First the two images are roughly aligned by using a transform
// initialization, then they are registered using a rigid transform, that in
// turn, is used to initialize a registration with an affine transform. The
// transform resulting from the affine registration is used as the bulk
// transform of a BSplineTransform. The deformable registration is
// computed, and finally the resulting transform is used to resample the moving
// image.
//
// Software Guide : EndLatex
#include "itkImageRegistrationMethod.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkTimeProbesCollectorBase.h"
//  Software Guide : BeginLatex
//
//  The following are the most relevant headers to this example.
//
//  \index{itk::VersorRigid3DTransform!header}
//  \index{itk::AffineTransform!header}
//  \index{itk::BSplineTransform!header}
//  \index{itk::RegularStepGradientDescentOptimizer!header}
//
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
#include "itkAffineTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkVersorRigid3DTransform.h"
// this new version does not know about bulk transform - useless!
#include "itkBSplineTransform.h"
// from old version 3.20
#include "itkBSplineDeformableTransform.h"

#include "itkRegularStepGradientDescentOptimizer.h"
// Software Guide : EndCodeSnippet
#include "itkBSplineDecompositionImageFilter.h"
#include "itkBSplineResampleImageFunction.h"
#include "itkCastImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkSqrtImageFilter.h"
#include "itkSquaredDifferenceImageFilter.h"
#include "itkTransformFileWriter.h"
#include <itkImageMaskSpatialObject.h>

#include "itkCastImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "json.hpp"
#include "metaCommand.h"

//#include <cmath>
#include <boost/filesystem.hpp>
#include <boost/timer/timer.hpp>

#include "itkDisplacementFieldJacobianDeterminantFilter.h"

using json = nlohmann::json;
using namespace boost::filesystem;
using namespace boost::timer;
bool verbose = false;

//  The following section of code implements a Command observer
//  used to monitor the evolution of the registration process.
//
#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command {
public:
  using Self = CommandIterationUpdate;
  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<Self>;
  itkNewMacro(Self);

protected:
  CommandIterationUpdate() = default;

public:
  using OptimizerType = itk::RegularStepGradientDescentOptimizer;
  using OptimizerPointer = const OptimizerType *;
  void Execute(itk::Object *caller, const itk::EventObject &event) override { Execute((const itk::Object *)caller, event); }
  void Execute(const itk::Object *object, const itk::EventObject &event) override {
    auto optimizer = static_cast<OptimizerPointer>(object);
    if (!(itk::IterationEvent().CheckEvent(&event))) {
      return;
    }
    if (verbose) {
      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer->GetValue() << "   ";
      std::cout << std::endl;
    }
  }
};

json resultJSON;

int main(int argc, char *argv[]) {

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(4);

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  command.SetName("Registration");
  command.SetDescription("Registration of intensity images based on the ITK toolkit.");
  command.SetVersion("0.0.1");
  command.AddField("fixedImageFile", "The fixed image", MetaCommand::STRING, true);
  command.AddField("movingImageFile", "The moving image", MetaCommand::STRING, true);
  command.AddField("outputImageFile", "Output moving image after elastic registration with fixed image.", MetaCommand::STRING, true);

  command.SetOption("fixedImageMaskFile", "m", false, "Specify a mask for the fixed image.");
  command.AddOptionField("fixedImageMaskFile", "mask", MetaCommand::STRING, true);

  command.SetOption("differenceOutputFile", "i", false, "Export the difference image.");
  command.AddOptionField("differenceOutputFile", "differenceOutputFile", MetaCommand::STRING, true);

  command.SetOption("differenceBeforeRegistration", "b", false, "Export the difference image before elastic registration to this file.");
  command.AddOptionField("differenceBeforeRegistration", "differenceBeforeRegistration", MetaCommand::STRING, true);

  command.SetOption("filenameForFinalTransformParameter", "t", false, "Export the transform as a file.");
  command.AddOptionField("filenameForFinalTransformParameter", "filenameForFinalTransformParameter", MetaCommand::STRING, true);

  command.SetOption("useExplicitPDFderivatives", "e", false, "Use explicit PDF derivatives.");

  command.SetOption("useCachingBSplineWeights", "c", false, "Use caching of weights.");

  command.SetOption("deformationField", "d", false, "Export the deformation field.");
  command.AddOptionField("deformationField", "deformationField", MetaCommand::STRING, true);

  command.SetOption("numberOfGridNodesInsideImageInOneDimensionCoarse", "c", false, "Number of grid nodes in coarse sampling.");
  command.AddOptionField("numberOfGridNodesInsideImageInOneDimensionCoarse", "numberOfGridNodesInsideImageInOneDimensionCoarse", MetaCommand::INT, true, "5");

  command.SetOption("numberOfGridNodesInsideImageInOneDimensionFine", "f", false, "Number of grid nodes in fine sampling.");
  command.AddOptionField("numberOfGridNodesInsideImageInOneDimensionFine", "numberOfGridNodesInsideImageInOneDimensionFine", MetaCommand::INT, true, "20");

  command.SetOption("maximumStepLength", "s", false, "Maximum step length.");
  command.AddOptionField("maximumStepLength", "maximumStepLength", MetaCommand::FLOAT, true, "10.0");

  command.SetOption("maximumNumberOfIterations", "r", false, "Maximum number of iterations.");
  command.AddOptionField("maximumNumberOfIterations", "maximumNumberOfIterations", MetaCommand::INT, true, "200");

  command.SetOption("maxNumberOfThreads", "q", false, "Maximum number of threads.");
  command.AddOptionField("maxNumberOfThreads", "maxNumberOfThreads", MetaCommand::INT, true, "4");

  command.SetOption("verbose", "w", false, "Print more verbose output");

  if (!command.Parse(argc, argv)) {
    return 1;
  }

  std::string fixedImageFile = command.GetValueAsString("fixedImageFile");   // argv[1]
  std::string movingImageFile = command.GetValueAsString("movingImageFile"); // argv[2]
  std::string outputImageFile = command.GetValueAsString("outputImageFile"); // argv[3]

  int maximumNumberOfIterations = 200; // argv[13]
  if (command.GetOptionWasSet("maximumNumberOfIterations"))
    maximumNumberOfIterations = command.GetValueAsInt("maximumNumberOfIterations", "maximumNumberOfIterations");

  float maximumStepLength = 10.0f; // argv[12]
  if (command.GetOptionWasSet("maximumStepLength"))
    maximumStepLength = command.GetValueAsFloat("maximumStepLength", "maximumStepLength");

  int numberOfGridNodesInsideImageInOneDimensionFine = 20; // argv[11]
  if (command.GetOptionWasSet("numberOfGridNodesInsideImageInOneDimensionFine"))
    numberOfGridNodesInsideImageInOneDimensionFine =
        command.GetValueAsInt("numberOfGridNodesInsideImageInOneDimensionFine", "numberOfGridNodesInsideImageInOneDimensionFine");

  int numberOfGridNodesInsideImageInOneDimensionCoarse = 5; // argv[10]
  if (command.GetOptionWasSet("numberOfGridNodesInsideImageInOneDimensionCoarse"))
    numberOfGridNodesInsideImageInOneDimensionCoarse =
        command.GetValueAsInt("numberOfGridNodesInsideImageInOneDimensionCoarse", "numberOfGridNodesInsideImageInOneDimensionCoarse");

  std::string deformationField; // argv[9]
  if (command.GetOptionWasSet("deformationField"))
    deformationField = command.GetValueAsString("deformationField", "deformationField");

  std::string fixedImageMaskFile;
  if (command.GetOptionWasSet("fixedImageMaskFile"))
    fixedImageMaskFile = command.GetValueAsString("fixedImageMaskFile", "mask");

  std::string filenameForFinalTransformParameter; // argv[6]
  if (command.GetOptionWasSet("filenameForFinalTransformParameter"))
    filenameForFinalTransformParameter = command.GetValueAsString("filenameForFinalTransformParameter", "filenameForFinalTransformParameter");

  std::string differenceOutputFile; // argv[4]
  if (command.GetOptionWasSet("differenceOutputFile")) {
    differenceOutputFile = command.GetValueAsString("differenceOutputFile", "differenceOutputFile");
  }

  std::string differenceBeforeRegistration; // argv[5]
  if (command.GetOptionWasSet("differenceBeforeRegistration"))
    differenceBeforeRegistration = command.GetValueAsString("differenceBeforeRegistration", "differenceBeforeRegistration");

  bool useCachingBSplineWeights = false; // argv[8]
  if (command.GetOptionWasSet("useCachingBSplineWeights"))
    useCachingBSplineWeights = true;

  bool useExplicitPDFderivatives = false; // argv[7]
  if (command.GetOptionWasSet("useExplicitPDFderivatives"))
    useExplicitPDFderivatives = true;

  if (command.GetOptionWasSet("maxNumberOfThreads"))
    itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(command.GetValueAsInt("maxNumberOfThreads", "maxNumberOfThreads"));

  // bool verbose = false; // global variable
  if (command.GetOptionWasSet("verbose")) {
    fprintf(stdout, "verbose level 1\n");
    fprintf(stdout, "use %d threads...\n", itk::MultiThreaderBase::GetGlobalMaximumNumberOfThreads());
    resultJSON["num_threads"] = itk::MultiThreaderBase::GetGlobalMaximumNumberOfThreads();
    verbose = true;
  }

  resultJSON["command_line"] = json::array();
  for (int i = 0; i < argc; i++) {
    resultJSON["command_line"].push_back(std::string(argv[i]));
  }

  cpu_timer timer;

  constexpr unsigned int ImageDimension = 3;
  using PixelType = signed short;
  using FixedImageType = itk::Image<PixelType, ImageDimension>;
  using MovingImageType = itk::Image<PixelType, ImageDimension>;
  const unsigned int SpaceDimension = ImageDimension;
  constexpr unsigned int SplineOrder = 3;
  using InternalPixelType = float;
  using InternalImageType = itk::Image<InternalPixelType, ImageDimension>;
  using CoordinateRepType = double;
  using RigidTransformType = itk::VersorRigid3DTransform<double>;
  using AffineTransformType = itk::AffineTransform<double, SpaceDimension>;

  // old from 3.20
  typedef itk::BSplineDeformableTransform<CoordinateRepType, SpaceDimension, SplineOrder> DeformableTransformType;

  /* using DeformableTransformType = itk::BSplineTransform<
                            CoordinateRepType,
                            SpaceDimension,
                            SplineOrder >; */
  using TransformInitializerType = itk::CenteredTransformInitializer<RigidTransformType, InternalImageType, InternalImageType>;
  using OptimizerType = itk::RegularStepGradientDescentOptimizer;
  using MetricType = /*  itk::MeanSquaresImageToImageMetric< */ itk::MattesMutualInformationImageToImageMetric<InternalImageType, InternalImageType>;
  using InterpolatorType = itk::LinearInterpolateImageFunction<InternalImageType, double>;
  using RegistrationType = itk::ImageRegistrationMethod<
      /* FixedImageType */ InternalImageType,
      /* MovingImageType */ InternalImageType>;
  MetricType::Pointer metric = MetricType::New();
  OptimizerType::Pointer optimizer = OptimizerType::New();
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  RegistrationType::Pointer registration = RegistrationType::New();
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
  registration->SetInterpolator(interpolator);
  // Auxiliary identity transform.
  using IdentityTransformType = itk::IdentityTransform<double, SpaceDimension>;
  IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();
  //
  //   Read the Fixed and Moving images.
  //
  using FixedImageReaderType = itk::ImageFileReader<FixedImageType>;
  using MovingImageReaderType = itk::ImageFileReader<MovingImageType>;
  FixedImageReaderType::Pointer fixedImageReader = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
  fixedImageReader->SetFileName(fixedImageFile /*  argv[1] */);
  movingImageReader->SetFileName(movingImageFile /*  argv[2] */);
  resultJSON["series_identifier"] = fixedImageFile /* argv[1] */;
  try {
    fixedImageReader->Update();
    movingImageReader->Update();
  } catch (itk::ExceptionObject &err) {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

  using ImageCasterType = itk::CastImageFilter<FixedImageType, InternalImageType>;

  if (verbose) {
    fprintf(stdout, "fixed image reader voxel size is: %f %f %f\n", fixedImageReader->GetOutput()->GetSpacing()[0],
            fixedImageReader->GetOutput()->GetSpacing()[1], fixedImageReader->GetOutput()->GetSpacing()[2]);
    fprintf(stdout, "fixed image reader origin is: %f %f %f\n", fixedImageReader->GetOutput()->GetOrigin()[0], fixedImageReader->GetOutput()->GetOrigin()[1],
            fixedImageReader->GetOutput()->GetOrigin()[2]);

    fprintf(stdout, "moving image reader voxel size is: %f %f %f\n", movingImageReader->GetOutput()->GetSpacing()[0],
            movingImageReader->GetOutput()->GetSpacing()[1], movingImageReader->GetOutput()->GetSpacing()[2]);
    fprintf(stdout, "moving image reader origin is: %f %f %f\n", movingImageReader->GetOutput()->GetOrigin()[0], movingImageReader->GetOutput()->GetOrigin()[1],
            movingImageReader->GetOutput()->GetOrigin()[2]);
  }

  ImageCasterType::Pointer fixedImageCaster = ImageCasterType::New();
  ImageCasterType::Pointer movingImageCaster = ImageCasterType::New();
  fixedImageCaster->SetInput(fixedImageReader->GetOutput());
  movingImageCaster->SetInput(movingImageReader->GetOutput());
  // match the histograms between source and target
  using MatchingFilterType = itk::HistogramMatchingImageFilter<InternalImageType, InternalImageType>;
  MatchingFilterType::Pointer matcher = MatchingFilterType::New();
  matcher->SetInput(movingImageCaster->GetOutput());
  matcher->SetReferenceImage(fixedImageCaster->GetOutput());
  matcher->SetNumberOfHistogramLevels(1024);
  matcher->SetNumberOfMatchPoints(7);
  matcher->ThresholdAtMeanIntensityOn();

  fixedImageCaster->Update();
  InternalImageType::Pointer fixedImage = fixedImageCaster->GetOutput();
  fixedImage->SetOrigin(fixedImageReader->GetOutput()->GetOrigin());
  fixedImage->SetSpacing(fixedImageReader->GetOutput()->GetSpacing());
  fixedImage->SetDirection(fixedImageReader->GetOutput()->GetDirection());

  movingImageCaster->Update();
  InternalImageType::Pointer movingImage = movingImageCaster->GetOutput();
  movingImage->SetOrigin(movingImageReader->GetOutput()->GetOrigin());
  movingImage->SetSpacing(movingImageReader->GetOutput()->GetSpacing());
  movingImage->SetDirection(movingImageReader->GetOutput()->GetDirection());

  // read in the Mask (if it has been provided)
  // typedef itk::ImageMaskSpatialObject<ImageDimension> MaskType;
  // MaskType::Pointer spatialMask = MaskType::New();
  using ImageMaskSpatialObjectType = itk::ImageMaskSpatialObject<ImageDimension>;
  typename ImageMaskSpatialObjectType::Pointer spatialMask = ImageMaskSpatialObjectType::New();

  typedef itk::Image<unsigned char, ImageDimension> ImageMaskType;
  typedef itk::ImageFileReader<ImageMaskType> MaskReaderType;
  MaskReaderType::Pointer maskReader = MaskReaderType::New();
  if (command.GetOptionWasSet("fixedImageMaskFile")) {
    if (verbose)
      fprintf(stdout, "reading in the fixed image mask...\n");
    maskReader->SetFileName(fixedImageMaskFile);
    maskReader->Update();
    // the mask should have to information of the fixed image
    // copy to make sure that is the case? Or should we trust the user?
    ImageMaskType::Pointer mm = maskReader->GetOutput();

    // lets try and create a mask similar to the Test for MattesMutualInformationImageToImageMetricTest.css
    // in that example the mask mimics the fixed image
    typename ImageMaskType::Pointer imgFixedMask = ImageMaskType::New();
    imgFixedMask->CopyInformation(fixedImage);
    typename FixedImageType::RegionType region = fixedImage->GetLargestPossibleRegion();
    imgFixedMask->SetRegions(region);
    imgFixedMask->Allocate(true); // initialize
    // copy the pixel over from the imported mm
    itk::ImageRegionIterator<ImageMaskType> maskIterator(mm, region);
    itk::ImageRegionIterator<ImageMaskType> fixedIterator(imgFixedMask, region);
    while (!maskIterator.IsAtEnd() && !fixedIterator.IsAtEnd()) {
      fixedIterator.Set(maskIterator.Get());
      ++fixedIterator;
      ++maskIterator;
    }

    // FixedImageType::RegionType fR = fixedImageReader->GetOutput()->GetBufferedRegion();
    // mm->SetRegions(fR);
    // mm->CopyInformation(fixedImage);
    // mm->SetRegions(region);
    // mm->SetOrigin(fixedImage->GetOrigin());
    // mm->SetSpacing(fixedImage->GetSpacing());
    // mm->SetDirection(fixedImage->GetDirection());

    spatialMask->SetImage(imgFixedMask);
    spatialMask->Update();
  }

  /*   fprintf(stdout, "fixed image caster voxel size is: %f %f %f\n", fixedImageCaster->GetOutput()->GetSpacing()[0],
          fixedImageCaster->GetOutput()->GetSpacing()[1],
          fixedImageCaster->GetOutput()->GetSpacing()[2]);
    fprintf(stdout, "fixed image caster origin is: %f %f %f\n",
          fixedImageCaster->GetOutput()->GetOrigin()[0],
          fixedImageCaster->GetOutput()->GetOrigin()[1],
          fixedImageCaster->GetOutput()->GetOrigin()[2]);

    fprintf(stdout, "moving image caster voxel size is: %f %f %f\n", movingImageCaster->GetOutput()->GetSpacing()[0],
          movingImageCaster->GetOutput()->GetSpacing()[1],
          movingImageCaster->GetOutput()->GetSpacing()[2]);
    fprintf(stdout, "moving image caster origin is: %f %f %f\n", movingImageCaster->GetOutput()->GetOrigin()[0],
          movingImageCaster->GetOutput()->GetOrigin()[1],
          movingImageCaster->GetOutput()->GetOrigin()[2]); */

  /*        regionGrowingField->SetOrigin(inputImage->GetOrigin());
          regionGrowingField->SetSpacing(inputImage->GetSpacing());
          regionGrowingField->SetDirection(inputImage->GetDirection());
   */

  registration->SetFixedImage(fixedImage);
  registration->SetMovingImage(/*  movingImageReader->GetOutput() */ movingImage /*  matcher->GetOutput()  */);
  //
  // Add a time and memory probes collector for profiling the computation time
  // of every stage.
  //
  itk::TimeProbesCollectorBase chronometer;
  itk::MemoryProbesCollectorBase memorymeter;
  //
  // Setup the metric parameters
  //
  metric->SetNumberOfHistogramBins(50);
  resultJSON["number_of_histogram_bins"] = 50;

  if (command.GetOptionWasSet("fixedImageMaskFile")) {
    if (verbose) {
      // where is the spatialMask?
      fprintf(stdout, "spatial mask has voxel size: %f %f %f\n", maskReader->GetOutput()->GetSpacing()[0], maskReader->GetOutput()->GetSpacing()[1],
              maskReader->GetOutput()->GetSpacing()[2]);
      fprintf(stdout, "spatial mask is at origin: %f %f %f\n", maskReader->GetOutput()->GetOrigin()[0], maskReader->GetOutput()->GetOrigin()[1],
              maskReader->GetOutput()->GetOrigin()[2]);
    }
    metric->SetFixedImageMask(spatialMask);
    resultJSON["use_fixed_image_mask"] = fixedImageMaskFile;
  }
  FixedImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();
  const unsigned int numberOfPixels = fixedRegion.GetNumberOfPixels();
  metric->ReinitializeSeed(42);
  if (useExplicitPDFderivatives) {
    // Define whether to calculate the metric derivative by explicitly
    // computing the derivatives of the joint PDF with respect to the Transform
    // parameters, or doing it by progressively accumulating contributions from
    // each bin in the joint PDF.
    metric->SetUseExplicitPDFDerivatives(1 /*  std::stoi( argv[7] ) */);
  }
  if (useCachingBSplineWeights) {
    // Define whether to cache the BSpline weights and indexes corresponding to
    // each one of the samples used to compute the metric. Enabling caching will
    // make the algorithm run faster but it will have a cost on the amount of memory
    // that needs to be allocated. This option is only relevant when using the
    // BSplineTransform.
    metric->SetUseCachingOfBSplineWeights(1 /*  std::stoi(argv[8]) */);
  }
  //
  //  Initialize a rigid transform by using Image Intensity Moments
  //
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
  initializer->SetTransform(rigidTransform);
  initializer->SetFixedImage(fixedImageCaster->GetOutput());
  initializer->SetMovingImage(movingImageCaster->GetOutput());
  initializer->MomentsOn();
  if (verbose) {
    std::cout << "Starting Rigid Transform Initialization " << std::endl;
    std::cout << timer.format() << '\n';
  }
  memorymeter.Start("Rigid Initialization");
  chronometer.Start("Rigid Initialization");
  initializer->InitializeTransform();
  chronometer.Stop("Rigid Initialization");
  memorymeter.Stop("Rigid Initialization");
  if (verbose) {
    std::cout << "Rigid Transform Initialization completed" << std::endl;
    std::cout << rigidTransform->GetParameters() << std::endl;
    std::cout << std::endl;
    std::cout << timer.format() << '\n';
  }
  registration->SetFixedImageRegion(fixedRegion);
  registration->SetInitialTransformParameters(rigidTransform->GetParameters());
  registration->SetTransform(rigidTransform);
  //
  //  Define optimizer normalization to compensate for different dynamic range
  //  of rotations and translations.
  //
  using OptimizerScalesType = OptimizerType::ScalesType;
  OptimizerScalesType optimizerScales(rigidTransform->GetNumberOfParameters());
  const double translationScale = 1.0 / 1000.0;
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = translationScale;
  optimizerScales[4] = translationScale;
  optimizerScales[5] = translationScale;
  optimizer->SetScales(optimizerScales);
  optimizer->SetMaximumStepLength(0.2000);
  optimizer->SetMinimumStepLength(0.0001);
  optimizer->SetNumberOfIterations(maximumNumberOfIterations * 2); // fixed to high number == 400
  //
  // The rigid transform has 6 parameters we use therefore a few samples to run
  // this stage.
  //
  // Regulating the number of samples in the Metric is equivalent to performing
  // multi-resolution registration because it is indeed a sub-sampling of the
  // image.
  metric->SetNumberOfSpatialSamples(80000L);
  //
  // Create the Command observer and register it with the optimizer.
  //
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();

  optimizer->AddObserver(itk::IterationEvent(), observer);
  if (verbose) {
    std::cout << "Starting Rigid Registration " << std::endl;
    std::cout << timer.format() << '\n';
  }
  try {
    memorymeter.Start("Rigid Registration");
    chronometer.Start("Rigid Registration");
    registration->Update();
    chronometer.Stop("Rigid Registration");
    memorymeter.Stop("Rigid Registration");
    if (verbose) {
      std::cout << "Optimizer stop condition = " << registration->GetOptimizer()->GetStopConditionDescription() << registration->GetLastTransformParameters()
                << std::endl;
    }
    std::ostringstream o;
    o << registration->GetLastTransformParameters();
    resultJSON["rigid_stop_condition"] = registration->GetOptimizer()->GetStopConditionDescription();
    resultJSON["rigid_transform"] = o.str();

  } catch (itk::ExceptionObject &err) {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }
  if (verbose) {
    std::cout << "Rigid Registration completed" << std::endl;
    std::cout << std::endl;
    std::cout << timer.format() << '\n';
  }
  rigidTransform->SetParameters(registration->GetLastTransformParameters());
  //
  //  Perform Affine Registration
  //
  AffineTransformType::Pointer affineTransform = AffineTransformType::New();
  affineTransform->SetCenter(rigidTransform->GetCenter());
  affineTransform->SetTranslation(rigidTransform->GetTranslation());
  affineTransform->SetMatrix(rigidTransform->GetMatrix());
  registration->SetTransform(affineTransform);
  registration->SetInitialTransformParameters(affineTransform->GetParameters());
  optimizerScales = OptimizerScalesType(affineTransform->GetNumberOfParameters());
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = 1.0;
  optimizerScales[4] = 1.0;
  optimizerScales[5] = 1.0;
  optimizerScales[6] = 1.0;
  optimizerScales[7] = 1.0;
  optimizerScales[8] = 1.0;
  optimizerScales[9] = translationScale;
  optimizerScales[10] = translationScale;
  optimizerScales[11] = translationScale;
  optimizer->SetScales(optimizerScales);
  optimizer->SetMaximumStepLength(0.1000);
  optimizer->SetMinimumStepLength(0.00001);
  optimizer->SetNumberOfIterations(20);
  //
  // The Affine transform has 12 parameters we use therefore more samples to run
  // this stage.
  //
  // Regulating the number of samples in the Metric is equivalent to performing
  // multi-resolution registration because it is indeed a sub-sampling of the
  // image.
  metric->SetNumberOfSpatialSamples(80000L);
  if (verbose) {
    std::cout << "Starting Affine Registration " << std::endl;
    std::cout << timer.format() << '\n';
  }
  try {
    memorymeter.Start("Affine Registration");
    chronometer.Start("Affine Registration");
    registration->Update();
    chronometer.Stop("Affine Registration");
    memorymeter.Stop("Affine Registration");

    std::ostringstream o;
    o << registration->GetLastTransformParameters();
    resultJSON["affine_stop_condition"] = registration->GetOptimizer()->GetStopConditionDescription();
    resultJSON["affine_transform"] = o.str();

  } catch (itk::ExceptionObject &err) {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }
  if (verbose) {
    std::cout << "Affine Registration completed" << std::endl;
    std::cout << registration->GetLastTransformParameters() << std::endl;
    std::cout << std::endl;
    std::cout << timer.format() << '\n';
  }
  affineTransform->SetParameters(registration->GetLastTransformParameters());

  // code from 3.20 that does know about bulk transforms

  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  //
  //  Perform Deformable Registration
  //
  DeformableTransformType::Pointer bsplineTransformCoarse = DeformableTransformType::New();

  unsigned int numberOfGridNodesInOneDimensionCoarse = numberOfGridNodesInsideImageInOneDimensionCoarse; // 5 default

  resultJSON["coarse_grid_nodes"] = numberOfGridNodesInOneDimensionCoarse;

  typedef DeformableTransformType::RegionType RegionType;
  RegionType bsplineRegion;
  RegionType::SizeType gridSizeOnImage;
  RegionType::SizeType gridBorderSize;
  RegionType::SizeType totalGridSize;

  gridSizeOnImage.Fill(numberOfGridNodesInOneDimensionCoarse);
  gridBorderSize.Fill(SplineOrder); // Border for spline order = 3 ( 1 lower, 2 upper )
  totalGridSize = gridSizeOnImage + gridBorderSize;

  bsplineRegion.SetSize(totalGridSize);

  typedef DeformableTransformType::SpacingType SpacingType;
  SpacingType spacing = fixedImage->GetSpacing();

  typedef DeformableTransformType::OriginType OriginType;
  OriginType origin = fixedImage->GetOrigin();

  FixedImageType::SizeType fixedImageSize = fixedRegion.GetSize();

  for (unsigned int r = 0; r < ImageDimension; r++) {
    spacing[r] *= static_cast<double>(fixedImageSize[r] - 1) / static_cast<double>(gridSizeOnImage[r] - 1);
  }

  FixedImageType::DirectionType gridDirection = fixedImage->GetDirection();
  SpacingType gridOriginOffset = gridDirection * spacing;

  OriginType gridOrigin = origin - gridOriginOffset;

  bsplineTransformCoarse->SetGridSpacing(spacing);
  bsplineTransformCoarse->SetGridOrigin(gridOrigin);
  bsplineTransformCoarse->SetGridRegion(bsplineRegion);
  bsplineTransformCoarse->SetGridDirection(gridDirection);

  bsplineTransformCoarse->SetBulkTransform(affineTransform);

  typedef DeformableTransformType::ParametersType ParametersType;

  unsigned int numberOfBSplineParameters = bsplineTransformCoarse->GetNumberOfParameters();

  optimizerScales = OptimizerScalesType(numberOfBSplineParameters);
  optimizerScales.Fill(1.0);

  optimizer->SetScales(optimizerScales);

  ParametersType initialDeformableTransformParameters(numberOfBSplineParameters);
  initialDeformableTransformParameters.Fill(0.0);

  bsplineTransformCoarse->SetParameters(initialDeformableTransformParameters);

  registration->SetInitialTransformParameters(bsplineTransformCoarse->GetParameters());
  registration->SetTransform(bsplineTransformCoarse);

  // Software Guide : EndCodeSnippet

  //  Software Guide : BeginLatex
  //
  //  Next we set the parameters of the RegularStepGradientDescentOptimizer object.
  //
  //  Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  optimizer->SetMaximumStepLength(maximumStepLength /* 10.0 */);
  optimizer->SetMinimumStepLength(0.01);
  resultJSON["optimizer_max_step_length"] = maximumStepLength;

  optimizer->SetRelaxationFactor(0.7);
  optimizer->SetNumberOfIterations(maximumNumberOfIterations /*  200 */);
  resultJSON["optimizer_number_iterations"] = maximumNumberOfIterations;

  // Software Guide : EndCodeSnippet

  // Optionally, get the step length from the command line arguments
  // if (argc > 11) {
  //  optimizer->SetMaximumStepLength(atof(argv[12]));
  //}

  // Optionally, get the number of iterations from the command line arguments
  // if (argc > 12) {
  //  optimizer->SetNumberOfIterations(atoi(argv[13]));
  //}

  //
  // The BSpline transform has a large number of parameters, we use therefore a
  // much larger number of samples to run this stage.
  //
  // Regulating the number of samples in the Metric is equivalent to performing
  // multi-resolution registration because it is indeed a sub-sampling of the
  // image.
  metric->SetNumberOfSpatialSamples(numberOfBSplineParameters * 100);
  resultJSON["optimizer_number_of_samples"] = numberOfBSplineParameters * 100;

  if (verbose) {
    std::cout << std::endl << "Starting Deformable Registration Coarse Grid" << std::endl;
    std::cout << timer.format() << '\n';
  }
  try {
    // itkProbesStart( "Deformable Registration Coarse" );
    // registration->StartRegistration();
    // itkProbesStop( "Deformable Registration Coarse" );

    memorymeter.Start("Deformable Registration Coarse");
    chronometer.Start("Deformable Registration Coarse");
    registration->Update();
    chronometer.Stop("Deformable Registration Coarse");
    memorymeter.Stop("Deformable Registration Coarse");

    std::ostringstream o;
    o << registration->GetLastTransformParameters();
    resultJSON["elastic_coarse_stop_condition"] = registration->GetOptimizer()->GetStopConditionDescription();
    if (verbose)
      resultJSON["elastic_coarse_transform"] = o.str();
  } catch (itk::ExceptionObject &err) {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

  if (verbose) {
    std::cout << "Deformable Registration Coarse Grid completed" << std::endl;
    std::cout << std::endl;
    std::cout << timer.format() << '\n';
  }

  OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();

  bsplineTransformCoarse->SetParameters(finalParameters);

  //  Software Guide : BeginLatex
  //
  //  Once the registration has finished with the low resolution grid, we
  //  proceed to instantiate a higher resolution
  //  \code{BSplineDeformableTransform}.
  //
  //  Software Guide : EndLatex

  DeformableTransformType::Pointer bsplineTransformFine = DeformableTransformType::New();

  unsigned int numberOfGridNodesInOneDimensionFine = numberOfGridNodesInsideImageInOneDimensionFine;

  // if (argc > 11) {
  //  numberOfGridNodesInOneDimensionFine = atoi(argv[11]);
  //}

  resultJSON["fine_grid_nodes"] = numberOfGridNodesInOneDimensionFine;

  RegionType::SizeType gridHighSizeOnImage;
  gridHighSizeOnImage.Fill(numberOfGridNodesInOneDimensionFine);
  totalGridSize = gridHighSizeOnImage + gridBorderSize;

  bsplineRegion.SetSize(totalGridSize);

  SpacingType spacingHigh = fixedImage->GetSpacing();
  OriginType originHigh = fixedImage->GetOrigin();

  for (unsigned int rh = 0; rh < ImageDimension; rh++) {
    spacingHigh[rh] *= static_cast<double>(fixedImageSize[rh] - 1) / static_cast<double>(gridHighSizeOnImage[rh] - 1);
    originHigh[rh] -= spacingHigh[rh];
  }

  SpacingType gridOriginOffsetHigh = gridDirection * spacingHigh;

  OriginType gridOriginHigh = origin - gridOriginOffsetHigh;

  bsplineTransformFine->SetGridSpacing(spacingHigh);
  bsplineTransformFine->SetGridOrigin(gridOriginHigh);
  bsplineTransformFine->SetGridRegion(bsplineRegion);
  bsplineTransformFine->SetGridDirection(gridDirection);

  bsplineTransformFine->SetBulkTransform(affineTransform);

  numberOfBSplineParameters = bsplineTransformFine->GetNumberOfParameters();

  ParametersType parametersHigh(numberOfBSplineParameters);
  parametersHigh.Fill(0.0);

  //  Software Guide : BeginLatex
  //
  //  Now we need to initialize the BSpline coefficients of the higher resolution
  //  transform. This is done by first computing the actual deformation field
  //  at the higher resolution from the lower resolution BSpline coefficients.
  //  Then a BSpline decomposition is done to obtain the BSpline coefficient of
  //  the higher resolution transform.
  //
  //  Software Guide : EndLatex

  unsigned int counter = 0;

  for (unsigned int k = 0; k < SpaceDimension; k++) {
    typedef DeformableTransformType::ImageType ParametersImageType;
    typedef itk::ResampleImageFilter<ParametersImageType, ParametersImageType> ResamplerType;
    ResamplerType::Pointer upsampler = ResamplerType::New();

    typedef itk::BSplineResampleImageFunction<ParametersImageType, double> FunctionType;
    FunctionType::Pointer function = FunctionType::New();

    upsampler->SetInput(bsplineTransformCoarse->GetCoefficientImages()[k]);
    upsampler->SetInterpolator(function);
    upsampler->SetTransform(identityTransform);
    upsampler->SetSize(bsplineTransformFine->GetGridRegion().GetSize());
    upsampler->SetOutputSpacing(bsplineTransformFine->GetGridSpacing());
    upsampler->SetOutputOrigin(bsplineTransformFine->GetGridOrigin());

    typedef itk::BSplineDecompositionImageFilter<ParametersImageType, ParametersImageType> DecompositionType;
    DecompositionType::Pointer decomposition = DecompositionType::New();

    decomposition->SetSplineOrder(SplineOrder);
    decomposition->SetInput(upsampler->GetOutput());
    decomposition->Update();

    ParametersImageType::Pointer newCoefficients = decomposition->GetOutput();

    // copy the coefficients into the parameter array
    typedef itk::ImageRegionIterator<ParametersImageType> Iterator;
    Iterator it(newCoefficients, bsplineTransformFine->GetGridRegion());
    while (!it.IsAtEnd()) {
      parametersHigh[counter++] = it.Get();
      ++it;
    }
  }

  optimizerScales = OptimizerScalesType(numberOfBSplineParameters);
  optimizerScales.Fill(1.0);

  optimizer->SetScales(optimizerScales);

  bsplineTransformFine->SetParameters(parametersHigh);

  //  Software Guide : BeginLatex
  //
  //  We now pass the parameters of the high resolution transform as the initial
  //  parameters to be used in a second stage of the registration process.
  //
  //  Software Guide : EndLatex

  if (verbose) {
    std::cout << "Starting Registration with high resolution transform" << std::endl;
    std::cout << timer.format() << '\n';
  }
  // Software Guide : BeginCodeSnippet
  registration->SetInitialTransformParameters(bsplineTransformFine->GetParameters());
  registration->SetTransform(bsplineTransformFine);

  //
  // The BSpline transform at fine scale has a very large number of parameters,
  // we use therefore a much larger number of samples to run this stage. In this
  // case, however, the number of transform parameters is closer to the number
  // of pixels in the image. Therefore we use the geometric mean of the two numbers
  // to ensure that the number of samples is larger than the number of transform
  // parameters and smaller than the number of samples.
  //
  // Regulating the number of samples in the Metric is equivalent to performing
  // multi-resolution registration because it is indeed a sub-sampling of the
  // image.
  const unsigned long numberOfSamples =
      static_cast<unsigned long>(std::sqrt(static_cast<double>(numberOfBSplineParameters) * static_cast<double>(numberOfPixels)));

  metric->SetNumberOfSpatialSamples(numberOfSamples);
  resultJSON["registation_fine_number_of_samples"] = numberOfSamples;

  try {
    memorymeter.Start("Deformable Registration Fine");
    chronometer.Start("Deformable Registration Fine");
    registration->Update();
    chronometer.Stop("Deformable Registration Fine");
    memorymeter.Stop("Deformable Registration Fine");

    std::ostringstream o;
    o << registration->GetLastTransformParameters();
    resultJSON["elastic_fine_stop_condition"] = registration->GetOptimizer()->GetStopConditionDescription();
    if (verbose) {
      std::cout << registration->GetOptimizer()->GetStopConditionDescription() << std::endl;
      resultJSON["elastic_fine_transform"] = o.str();
    }

    //    itkProbesStart( "Deformable Registration Fine" );
    //    registration->StartRegistration();
    //    itkProbesStop( "Deformable Registration Fine" );
  } catch (itk::ExceptionObject &err) {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }
  // Software Guide : EndCodeSnippet
  if (verbose) {
    std::cout << "Deformable Registration Fine Grid completed" << std::endl;
    std::cout << timer.format() << '\n';
  }
  // Report the time and memory taken by the registration
  // itkProbesReport( std::cout );

  finalParameters = registration->GetLastTransformParameters();

  bsplineTransformFine->SetParameters(finalParameters);

  // here is the newer code that does not know about bulktransform ///
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  /*
    //
    //  Perform Deformable Registration
    //
    DeformableTransformType::Pointer  bsplineTransformCoarse = DeformableTransformType::New();
    unsigned int numberOfGridNodesInOneDimensionCoarse = 5;
    DeformableTransformType::PhysicalDimensionsType   fixedPhysicalDimensions;
    DeformableTransformType::MeshSizeType             meshSize;
    DeformableTransformType::OriginType               fixedOrigin;
    for( unsigned int i=0; i< SpaceDimension; i++ )
      {
      fixedOrigin[i] = fixedImage->GetOrigin()[i];
      fixedPhysicalDimensions[i] = fixedImage->GetSpacing()[i] *
        static_cast<double>(
        fixedImage->GetLargestPossibleRegion().GetSize()[i] - 1 );
      }
    meshSize.Fill( numberOfGridNodesInOneDimensionCoarse - SplineOrder );
    bsplineTransformCoarse->SetTransformDomainOrigin( fixedOrigin );
    bsplineTransformCoarse->SetTransformDomainPhysicalDimensions(
      fixedPhysicalDimensions );
    bsplineTransformCoarse->SetTransformDomainMeshSize( meshSize );
    bsplineTransformCoarse->SetTransformDomainDirection(
      fixedImage->GetDirection() );

    // is this missing???
    //bsplineTransformCoarse->SetBulkTransformMatrix(affineTransform);

    using ParametersType = DeformableTransformType::ParametersType;
    unsigned int numberOfBSplineParameters = bsplineTransformCoarse->GetNumberOfParameters();
    optimizerScales = OptimizerScalesType( numberOfBSplineParameters );
    optimizerScales.Fill( 1.0 );
    optimizer->SetScales( optimizerScales );
    ParametersType initialDeformableTransformParameters( numberOfBSplineParameters );
    initialDeformableTransformParameters.Fill( 0.0 );
    bsplineTransformCoarse->SetParameters( initialDeformableTransformParameters );
    registration->SetInitialTransformParameters( bsplineTransformCoarse->GetParameters() );
    registration->SetTransform( bsplineTransformCoarse );
    // Software Guide : EndCodeSnippet
    //  Software Guide : BeginLatex
    //
    //  Next we set the parameters of the RegularStepGradientDescentOptimizer object.
    //
    //  Software Guide : EndLatex
    // Software Guide : BeginCodeSnippet
    optimizer->SetMaximumStepLength( 0.1 );
    optimizer->SetMinimumStepLength(  0.001 );
    optimizer->SetRelaxationFactor( 0.6 );
    optimizer->SetNumberOfIterations( 80 );
    // Software Guide : EndCodeSnippet
    // Optionally, get the step length from the command line arguments
    if( argc > 11 )
      {
      optimizer->SetMaximumStepLength( std::stod( argv[12] ) );
      }
    // Optionally, get the number of iterations from the command line arguments
    if( argc > 12 )
      {
      optimizer->SetNumberOfIterations( std::stoi( argv[13] ) );
      }
    //
    // The BSpline transform has a large number of parameters, we use therefore a
    // much larger number of samples to run this stage.
    //
    // Regulating the number of samples in the Metric is equivalent to performing
    // multi-resolution registration because it is indeed a sub-sampling of the
    // image.
    metric->SetNumberOfSpatialSamples( numberOfBSplineParameters * 100 );
    std::cout << std::endl << "Starting Deformable Registration Coarse Grid" << std::endl;
    try
      {
      memorymeter.Start( "Deformable Registration Coarse" );
      chronometer.Start( "Deformable Registration Coarse" );
      registration->Update();
      chronometer.Stop( "Deformable Registration Coarse" );
      memorymeter.Stop( "Deformable Registration Coarse" );
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
      }
    std::cout << "Deformable Registration Coarse Grid completed" << std::endl;
    std::cout << std::endl;
    OptimizerType::ParametersType finalParameters =
                      registration->GetLastTransformParameters();
    bsplineTransformCoarse->SetParameters( finalParameters );
    //  Software Guide : BeginLatex
    //
    //  Once the registration has finished with the low resolution grid, we
    //  proceed to instantiate a higher resolution
    //  \code{BSplineTransform}.
    //
    //  Software Guide : EndLatex
    DeformableTransformType::Pointer  bsplineTransformFine = DeformableTransformType::New();
    unsigned int numberOfGridNodesInOneDimensionFine = 20;
    meshSize.Fill( numberOfGridNodesInOneDimensionFine - SplineOrder );
    bsplineTransformFine->SetTransformDomainOrigin( fixedOrigin );
    bsplineTransformFine->SetTransformDomainPhysicalDimensions(
      fixedPhysicalDimensions );
    bsplineTransformFine->SetTransformDomainMeshSize( meshSize );
    bsplineTransformFine->SetTransformDomainDirection(
      fixedImage->GetDirection() );
    numberOfBSplineParameters = bsplineTransformFine->GetNumberOfParameters();
    ParametersType parametersHigh( numberOfBSplineParameters );
    parametersHigh.Fill( 0.0 );
    //  Software Guide : BeginLatex
    //
    //  Now we need to initialize the BSpline coefficients of the higher resolution
    //  transform. This is done by first computing the actual deformation field
    //  at the higher resolution from the lower resolution BSpline coefficients.
    //  Then a BSpline decomposition is done to obtain the BSpline coefficient of
    //  the higher resolution transform.
    //
    //  Software Guide : EndLatex
    unsigned int counter = 0;
    for ( unsigned int k = 0; k < SpaceDimension; k++ )
      {
      using ParametersImageType = DeformableTransformType::ImageType;
      using ResamplerType = itk::ResampleImageFilter<ParametersImageType,ParametersImageType>;
      ResamplerType::Pointer upsampler = ResamplerType::New();
      using FunctionType = itk::BSplineResampleImageFunction<ParametersImageType,double>;
      FunctionType::Pointer function = FunctionType::New();
      upsampler->SetInput( bsplineTransformCoarse->GetCoefficientImages()[k] );
      upsampler->SetInterpolator( function );
      upsampler->SetTransform( identityTransform );
      upsampler->SetSize( bsplineTransformFine->GetCoefficientImages()[k]->
        GetLargestPossibleRegion().GetSize() );
      upsampler->SetOutputSpacing( bsplineTransformFine->GetCoefficientImages()[k]->
        GetSpacing() );
      upsampler->SetOutputOrigin( bsplineTransformFine->GetCoefficientImages()[k]->
        GetOrigin() );
      using DecompositionType =
          itk::BSplineDecompositionImageFilter<ParametersImageType,ParametersImageType>;
      DecompositionType::Pointer decomposition = DecompositionType::New();
      decomposition->SetSplineOrder( SplineOrder );
      decomposition->SetInput( upsampler->GetOutput() );
      decomposition->Update();
      ParametersImageType::Pointer newCoefficients = decomposition->GetOutput();
      // copy the coefficients into the parameter array
      using Iterator = itk::ImageRegionIterator<ParametersImageType>;
      Iterator it( newCoefficients, bsplineTransformFine->GetCoefficientImages()[k]->
        GetLargestPossibleRegion() );
      while ( !it.IsAtEnd() )
        {
        parametersHigh[ counter++ ] = it.Get();
        ++it;
        }
      }
    optimizerScales = OptimizerScalesType( numberOfBSplineParameters );
    optimizerScales.Fill( 1.0 );
    optimizer->SetScales( optimizerScales );
    bsplineTransformFine->SetParameters( parametersHigh );
    //  Software Guide : BeginLatex
    //
    //  We now pass the parameters of the high resolution transform as the initial
    //  parameters to be used in a second stage of the registration process.
    //
    //  Software Guide : EndLatex
    std::cout << "Starting Registration with high resolution transform" << std::endl;
    // Software Guide : BeginCodeSnippet
    registration->SetInitialTransformParameters(
                                        bsplineTransformFine->GetParameters() );
    registration->SetTransform( bsplineTransformFine );
    //
    // The BSpline transform at fine scale has a very large number of parameters,
    // we use therefore a much larger number of samples to run this stage. In
    // this case, however, the number of transform parameters is closer to the
    // number of pixels in the image. Therefore we use the geometric mean of the
    // two numbers to ensure that the number of samples is larger than the number
    // of transform parameters and smaller than the number of samples.
    //
    // Regulating the number of samples in the Metric is equivalent to performing
    // multi-resolution registration because it is indeed a sub-sampling of the
    // image.
    const auto numberOfSamples = static_cast<unsigned long>(
         std::sqrt( static_cast<double>( numberOfBSplineParameters ) *
                   static_cast<double>( numberOfPixels ) ) );
    metric->SetNumberOfSpatialSamples( numberOfSamples );
    try
      {
      memorymeter.Start( "Deformable Registration Fine" );
      chronometer.Start( "Deformable Registration Fine" );
      registration->Update();
      chronometer.Stop( "Deformable Registration Fine" );
      memorymeter.Stop( "Deformable Registration Fine" );
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
      }
    // Software Guide : EndCodeSnippet
    std::cout << "Deformable Registration Fine Grid completed" << std::endl;
    std::cout << std::endl;
    // Report the time and memory taken by the registration
    chronometer.Report( std::cout );
    memorymeter.Report( std::cout );
    finalParameters = registration->GetLastTransformParameters();
    bsplineTransformFine->SetParameters( finalParameters );


    */

  using ResampleFilterType = itk::ResampleImageFilter<InternalImageType, InternalImageType>;
  ResampleFilterType::Pointer resample = ResampleFilterType::New();
  resample->SetTransform(bsplineTransformFine);
  resample->SetInput(movingImageCaster->GetOutput());
  resample->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
  resample->SetOutputOrigin(fixedImage->GetOrigin());
  resample->SetOutputSpacing(fixedImage->GetSpacing());
  resample->SetOutputDirection(fixedImage->GetDirection());
  // This value is set to zero in order to make easier to perform
  // regression testing in this example. However, for didactic
  // exercise it will be better to set it to a medium gray value
  // such as 100 or 128.
  resample->SetDefaultPixelValue(-1024);
  using OutputPixelType = signed short;
  using OutputImageType = itk::Image<OutputPixelType, ImageDimension>;
  using CastFilterType = itk::CastImageFilter<InternalImageType, OutputImageType>;
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  WriterType::Pointer writer = WriterType::New();
  CastFilterType::Pointer caster = CastFilterType::New();
  writer->SetFileName(outputImageFile /* argv[3] */);
  // create the path if it does not exist
  path p(outputImageFile /* argv[3] */);
  create_directories(p.parent_path());
  resultJSON["resampled_moving_image"] = outputImageFile /* argv[3] */;

  caster->SetInput(resample->GetOutput());
  writer->SetInput(caster->GetOutput());
  if (verbose) {
    std::cout << "Writing resampled moving image...";
    std::cout << timer.format() << '\n';
  }
  try {
    writer->Update();
  } catch (itk::ExceptionObject &err) {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }
  if (verbose) {
    std::cout << " Done!" << std::endl;
    std::cout << timer.format() << '\n';
  }
  resultJSON["resampled_moving_image"] = outputImageFile /* argv[3] */;
  using DifferenceFilterType = itk::SquaredDifferenceImageFilter<InternalImageType, InternalImageType, InternalImageType>;
  DifferenceFilterType::Pointer difference = DifferenceFilterType::New();
  using SqrtFilterType = itk::SqrtImageFilter<InternalImageType, InternalImageType>;
  SqrtFilterType::Pointer sqrtFilter = SqrtFilterType::New();
  sqrtFilter->SetInput(difference->GetOutput());
  using DifferenceImageWriterType = itk::ImageFileWriter<InternalImageType>;
  DifferenceImageWriterType::Pointer writer2 = DifferenceImageWriterType::New();
  InternalImageType::Pointer sqrtF = sqrtFilter->GetOutput();
  sqrtF->SetOrigin(fixedImage->GetOrigin());
  sqrtF->SetSpacing(fixedImage->GetSpacing());
  sqrtF->SetDirection(fixedImage->GetDirection());
  writer2->SetInput(sqrtF);
  // Compute the difference image between the
  // fixed and resampled moving image.
  if (command.GetOptionWasSet("differenceOutputFile")) {
    difference->SetInput1(fixedImageCaster->GetOutput());
    difference->SetInput2(resample->GetOutput());
    writer2->SetFileName(differenceOutputFile /* argv[4] */);

    // create the path if it does not exist
    path p(differenceOutputFile /* argv[4] */);
    create_directories(p.parent_path());
    resultJSON["difference_image_after_elastic"] = differenceOutputFile /* argv[4] */;

    if (verbose) {
      std::cout << "Writing difference image after registration...";
      std::cout << timer.format() << '\n';
    }
    try {
      writer2->Update();
    } catch (itk::ExceptionObject &err) {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
    if (verbose) {
      std::cout << " Done!" << std::endl;
      std::cout << timer.format() << '\n';
    }
  }
  // Compute the difference image between the
  // affine registered and elastic registered image.
  if (command.GetOptionWasSet("differenceBeforeRegistration")) {
    writer2->SetFileName(differenceBeforeRegistration /* argv[5] */);

    // create the path if it does not exist
    path p(differenceBeforeRegistration /* argv[5] */);
    create_directories(p.parent_path());

    InternalImageType::Pointer sqrtF = fixedImageCaster->GetOutput();
    sqrtF->SetOrigin(fixedImage->GetOrigin());
    sqrtF->SetSpacing(fixedImage->GetSpacing());
    sqrtF->SetDirection(fixedImage->GetDirection());

    difference->SetInput1(sqrtF);
    resample->SetTransform(affineTransform);
    if (verbose) {
      std::cout << "Writing difference image before elastic registration...";
      std::cout << timer.format() << '\n';
    }
    resultJSON["difference_image_before_elastic"] = differenceBeforeRegistration /* argv[5] */;

    try {
      writer2->Update();
    } catch (itk::ExceptionObject &err) {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
    if (verbose) {
      std::cout << " Done!" << std::endl;
      std::cout << timer.format() << '\n';
    }
  }
  // Generate the explicit deformation field resulting from
  // the registration.
  if (command.GetOptionWasSet("deformationField")) {
    if (verbose) {
      std::cout << "Create displacement field for export... ";
      std::cout << timer.format() << '\n';
    }
    using VectorType = itk::Vector<float, ImageDimension>;
    using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;
    DisplacementFieldType::Pointer field = DisplacementFieldType::New();
    field->SetRegions(fixedRegion);
    field->SetOrigin(fixedImage->GetOrigin());
    field->SetSpacing(fixedImage->GetSpacing());
    field->SetDirection(fixedImage->GetDirection());
    field->Allocate();
    using FieldIterator = itk::ImageRegionIterator<DisplacementFieldType>;
    FieldIterator fi(field, fixedRegion);
    fi.GoToBegin();
    DeformableTransformType::InputPointType fixedPoint;
    DeformableTransformType::OutputPointType movingPoint;
    DisplacementFieldType::IndexType index;
    VectorType displacement;
    while (!fi.IsAtEnd()) {
      index = fi.GetIndex();
      field->TransformIndexToPhysicalPoint(index, fixedPoint);
      movingPoint = bsplineTransformFine->TransformPoint(fixedPoint);
      displacement = movingPoint - fixedPoint;
      fi.Set(displacement);
      ++fi;
    }
    if (verbose) {
      std::cout << "Done!" << std::endl;
      std::cout << timer.format() << '\n';
    }
    using FieldWriterType = itk::ImageFileWriter<DisplacementFieldType>;
    FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
    fieldWriter->SetInput(field);
    fieldWriter->SetFileName(deformationField /* argv[9] */);
    resultJSON["deformation_field"] = deformationField /* argv[9] */;

    // create the path if it does not exist
    path p(deformationField /* argv[9] */);
    create_directories(p.parent_path());

    if (verbose) {
      std::cout << "Writing deformation field ...";
      std::cout << timer.format() << '\n';
    }
    try {
      fieldWriter->Update();
    } catch (itk::ExceptionObject &excp) {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    if (verbose) {
      std::cout << " Done!" << std::endl;
      std::cout << timer.format() << '\n';
    }

    // we should write out the Jacobian of the deformation field as well
    // we can use this filter: itkDisplacementFieldJacobianDeterminantFilter.h
    // deformation field is: field
    if (1) {
      using JacobianFilterType = itk::DisplacementFieldJacobianDeterminantFilter<DisplacementFieldType, float, InternalImageType>;
      JacobianFilterType::Pointer computeJacobian = JacobianFilterType::New();
      computeJacobian->SetInput(field);
      typedef itk::ImageFileWriter<InternalImageType> WriterType;
      WriterType::Pointer writer = WriterType::New();

      std::string ofn = std::string(deformationField /* argv[9] */);
      std::string ofn2 = ofn;
      size_t lastdot = ofn.find_last_of(".");
      if (lastdot == std::string::npos)
        ofn2 = ofn + "_walls.nii";
      else
        ofn2 = ofn.substr(0, lastdot) + "_jacobian.nii";
      resultJSON["deformation_field_jacobian"] = ofn2;
      std::string volFileName = ofn2;
      path p(volFileName);
      create_directories(p.parent_path());

      // std::string a(labelfieldfilename + "trachea.nii");
      writer->SetFileName(volFileName);
      InternalImageType::Pointer cJ = computeJacobian->GetOutput();
      cJ->SetSpacing(fixedImage->GetSpacing());
      cJ->SetOrigin(fixedImage->GetOrigin());
      cJ->SetDirection(fixedImage->GetDirection());
      writer->SetInput(cJ);

      if (verbose) {
        std::cout << "Writing the jacobian scalar image as " << std::endl;
        std::cout << volFileName << std::endl << std::endl;
        std::cout << timer.format() << '\n';
      }
      resultJSON["jacobian"] = volFileName;

      try {
        writer->Update();
      } catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
      }
      if (verbose) {
        std::cout << "Done!" << std::endl;
        std::cout << timer.format() << '\n';
      }
    }
  }
  // Optionally, save the transform parameters in a file
  if (command.GetOptionWasSet("filenameForFinalTransformParameter")) {
    if (verbose) {
      std::cout << "Writing transform parameter file ...";
      std::cout << timer.format() << '\n';
    }
    using TransformWriterType = itk::TransformFileWriter;
    TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter->AddTransform(bsplineTransformFine);
    transformWriter->SetFileName(filenameForFinalTransformParameter /* argv[6] */);
    resultJSON["transform_parameter_file"] = filenameForFinalTransformParameter /*  argv[6] */;
    transformWriter->Update();
    if (verbose) {
      std::cout << " Done!" << std::endl;
      std::cout << timer.format() << '\n';
    }
  }
  resultJSON["runtime"] = timer.format();

  std::string res = resultJSON.dump(4) + "\n";
  std::ostringstream o;
  std::string si(resultJSON["series_identifier"]);
  si.erase(std::remove(si.begin(), si.end(), '\"'), si.end());
  boost::filesystem::path p2(outputImageFile /* argv[3] */);
  // std::string output = p.parent_path();
  o << p2.parent_path() << "/" << si << ".json";
  std::ofstream out(o.str());
  out << res;
  out.close();

  fprintf(stdout, "%s", res.c_str());

  return EXIT_SUCCESS;
}
