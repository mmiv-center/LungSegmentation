// todo: define attenuation coefficients for the diffent tissues types
// todo: use a low frequency dataset as a mask for high frequency structures (cells in organs)

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include "itkMetaDataObject.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
//#include "itkExtractImageFilter.h"
//#include "itkPasteImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkImageAdaptor.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkRGBPixel.h"
#include "itkSliceBySliceImageFilter.h"
#include "itkSymmetricEigenAnalysisImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
//#include <itkPixelAccessor.h>
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include "itkBSplineInterpolateImageFunction.h"
#include "itkExtractImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkWindowedSincInterpolateImageFunction.h"

#include "gdcmAnonymizer.h"
#include "gdcmAttribute.h"
#include "gdcmDataSetHelper.h"
#include "gdcmFileDerivation.h"
#include "gdcmFileExplicitFilter.h"
#include "gdcmGlobal.h"
#include "gdcmImageApplyLookupTable.h"
#include "gdcmImageChangePlanarConfiguration.h"
#include "gdcmImageChangeTransferSyntax.h"
#include "gdcmImageHelper.h"
#include "gdcmImageReader.h"
#include "gdcmImageWriter.h"
#include "gdcmMediaStorage.h"
#include "gdcmRescaler.h"
#include "gdcmStringFilter.h"
#include "gdcmUIDGenerator.h"
#include "itkConstantPadImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkGDCMImageIO.h"

#include "itkMetaDataDictionary.h"
#include "json.hpp"
#include "metaCommand.h"
#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <map>

#include <stdlib.h>
#include <time.h>

using json = nlohmann::json;
using namespace boost::filesystem;

template <typename TFilter> class CommandIterationUpdate : public itk::Command {
public:
  typedef CommandIterationUpdate Self;
  typedef itk::Command Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro(Self);

protected:
  CommandIterationUpdate() {}

public:
  virtual void Execute(itk::Object *caller, const itk::EventObject &event) ITK_OVERRIDE { Execute((const itk::Object *)caller, event); }

  virtual void Execute(const itk::Object *object, const itk::EventObject &event) ITK_OVERRIDE {
    const TFilter *filter = dynamic_cast<const TFilter *>(object);

    if (typeid(event) != typeid(itk::IterationEvent)) {
      return;
    }
    if (filter->GetElapsedIterations() == 1) {
      std::cout << "Current level = " << filter->GetCurrentLevel() + 1 << std::endl;
    }
    std::cout << "  Iteration " << filter->GetElapsedIterations() << " (of " << filter->GetMaximumNumberOfIterations()[filter->GetCurrentLevel()] << ").  ";
    std::cout << " Current convergence value = " << filter->GetCurrentConvergenceMeasurement() << " (threshold = " << filter->GetConvergenceThreshold() << ")"
              << std::endl;
  }
};

json resultJSON;

int main(int argc, char *argv[]) {
  boost::posix_time::ptime timeLocal = boost::posix_time::second_clock::local_time();
  resultJSON["run_date_time"] = to_simple_string(timeLocal);

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(4);

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  command.SetDescription("FakeLungVolumes simulating CT contrast volumes at high resolution. Exports volume files like nrrd or nifti based on the "
                         "provided file extension for the output image. The algorithm calculates the discrete intersection of two iso-surfaces "
                         "from band-pass filtered white noise volumes. Based on the amount of band-pass filtering and the threshold for the "
                         "intersection detection blood-vessel like pattern can be generated that are densely packed around alveoli like structures.");
  command.AddField("outfile", "Exported file name.", MetaCommand::STRING, true);

  command.SetOption("Resolution", "r", false, "Specify the resolution of the volume to be generated (in pixel as in 64x64x64).");
  command.AddOptionField("Resolution", "resolution", MetaCommand::STRING, false);

  command.SetOption("SmoothingKernelSize", "k", false, "Specify the kernel size for the Gaussian in pixel (7).");
  command.AddOptionField("SmoothingKernelSize", "kernelSize", MetaCommand::INT, false);

  command.SetOption("SmoothingIterations", "i", false, "Specify the number of times the Gaussian kernels are applied (2).");
  command.AddOptionField("SmoothingIterations", "iterations", MetaCommand::INT, false);

  command.SetOption("Threshold", "t", false, "Specify the threshold for zero-crossing (0.0001).");
  command.AddOptionField("Threshold", "threshold", MetaCommand::FLOAT, false);

  // we could use a pair of values x,y for the threshold at each of the two fields, that would rotate the
  // vessels - r, theta
  command.SetOption("Zero", "z", false, "Specify at what value the intersection should be calculated (0).");
  command.AddOptionField("Zero", "zero", MetaCommand::FLOAT, false);

  command.SetOption("finalSmooth", "f", false, "Specify the kernel size of a smoothing with a Gaussian at the end of the process (0).");
  command.AddOptionField("finalSmooth", "finalsmooth", MetaCommand::FLOAT, false);

  command.SetOption("VoidSpaces", "w", false,
                    "Create void spaces with a given distance away from the lines. Default is that this option is not used.\nIn the resulting volume 0 will be "
                    "the gap space right next to each vessel (label 4095) with 1, 2, 3, 4 the values of voxel that are in void space.");
  command.AddOptionField("VoidSpaces", "voidspaces", MetaCommand::FLOAT, false);

  command.SetOption("addLesion", "l", false, "Specify a lesion of a specific size (5). Requires the option VoidSpaces.");
  command.AddOptionField("addLesion", "addlesion", MetaCommand::INT, false);

  command.SetOption("Mask", "m", false, "Specify a mask file (assumption is that the mask fits in resolution with the volume created).");
  command.AddOptionField("Mask", "mask", MetaCommand::STRING, false);

  command.SetOption("randomSeed", "s", false,
                    "Specify the value used for initialization of the random numbers (time based). The same value should produce the same fields.");
  command.AddOptionField("randomSeed", "randomseed", MetaCommand::INT, false);

  command.SetOption("Force", "f", false, "Ignore existing files and force overwrite.");

  command.SetOption("Verbose", "V", false, "Print more verbose output");

  if (!command.Parse(argc, argv)) {
    return 1;
  }

  int randomSeed = 0;
  if (command.GetOptionWasSet("randomSeed")) {
    randomSeed = command.GetValueAsInt("randomSeed", "randomseed");
    fprintf(stdout, "random seed value is: %d\n", randomSeed);
  }

  float finalSmooth = 0;
  if (command.GetOptionWasSet("finalSmooth")) {
    finalSmooth = command.GetValueAsFloat("finalSmooth", "finalsmooth");
    fprintf(stdout, "final smoothing kernel is: %f\n", finalSmooth);
  }

  int iterations = 2;
  if (command.GetOptionWasSet("Iterations")) {
    iterations = command.GetValueAsInt("Iterations", "iterations");
    fprintf(stdout, "iterations is now: %d\n", iterations);
  }

  float threshold = 0.0001;
  if (command.GetOptionWasSet("Threshold")) {
    threshold = command.GetValueAsFloat("Threshold", "threshold");
    fprintf(stdout, "threshold is now: %f\n", threshold);
  }

  std::string maskName = ""; // we should make sure that we have a way to specify what label to use in the mask as well
  // we could use a <filename>,label1,label2,... way to specify options. Or another argument or list of arguments. For
  // now just use all labels.
  if (command.GetOptionWasSet("Mask")) {
    maskName = command.GetValueAsString("Mask", "mask");
    fprintf(stdout, "Mask volume requested: %s\n", maskName.c_str());
  }

  float zero = 0;
  if (command.GetOptionWasSet("Zero")) {
    zero = command.GetValueAsFloat("Zero", "zero");
    fprintf(stdout, "zero-crossing at: %f\n", zero);
  }

  float voidSpaces = 0.0001;
  if (command.GetOptionWasSet("VoidSpaces")) {
    voidSpaces = command.GetValueAsFloat("VoidSpaces", "voidspaces");
    fprintf(stdout, "void spaces at a distance of: %f\n", voidSpaces);
  }

  int smoothingKernelSize = 7;
  if (command.GetOptionWasSet("SmoothingKernelSize")) {
    smoothingKernelSize = command.GetValueAsInt("SmoothingKernelSize", "kernelSize");
    fprintf(stdout, "kernel size is now: %d\n", smoothingKernelSize);
  }
  std::string output = command.GetValueAsString("outfile");

  bool verbose = false;
  if (command.GetOptionWasSet("Verbose"))
    verbose = true;
  bool force = false;
  if (command.GetOptionWasSet("Force"))
    force = true;

  std::string resolution = "64x64x64";
  if (command.GetOptionWasSet("Resolution")) {
    resolution = command.GetValueAsString("Resolution", "resolution");
  }

  // store information in the result json file
  resultJSON["command_line"] = json::array();
  for (int i = 0; i < argc; i++) {
    resultJSON["command_line"].push_back(std::string(argv[i]));
  }

  typedef unsigned short OutputPixelType;
  typedef float FloatPixelType;
  const unsigned int Dimension = 3;

  //  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef itk::Image<FloatPixelType, Dimension> ImageType;
  typedef itk::Image<OutputPixelType, Dimension> OutputImageType;

  // mask stuff
  OutputImageType::Pointer mask;
  if (command.GetOptionWasSet("Mask")) {
    typedef itk::ImageFileReader<OutputImageType> MaskReaderType;
    MaskReaderType::Pointer maskreader = MaskReaderType::New();
    maskreader->SetFileName(maskName);
    try {
      maskreader->Update();
      mask = maskreader->GetOutput();
      mask->DisconnectPipeline();
    } catch (...) {
      mask = ITK_NULLPTR;
    }
  }

  ImageType::Pointer imageA = ImageType::New();
  ImageType::Pointer imageB = ImageType::New();
  ImageType::IndexType start;
  start[0] = 0; // first index on X
  start[1] = 0; // first index on Y
  start[2] = 0; // first index on Z

  ImageType::SizeType size;
  sscanf(resolution.c_str(), "%lux%lux%lu", &(size[0]), &(size[1]), &(size[2]));
  fprintf(stdout, "generate volume with: %lu %lu %lu voxel\n", size[0], size[1], size[2]);
  ImageType::RegionType region;

  region.SetSize(size);
  region.SetIndex(start);
  imageA->SetRegions(region);
  imageA->Allocate();
  imageB->SetRegions(region);
  imageB->Allocate();

  // set voxel values to random between -0.5 and 0.5
  using IteratorType = itk::ImageRegionIterator<ImageType>;
  IteratorType IteratorA(imageA, imageA->GetLargestPossibleRegion());
  IteratorType IteratorB(imageB, imageB->GetLargestPossibleRegion());
  if (command.GetOptionWasSet("randomSeed")) {
    srand(randomSeed);
  } else {
    srand(time(NULL));
  }
  for (IteratorA.GoToBegin(), IteratorB.GoToBegin(); !IteratorA.IsAtEnd() && !IteratorB.IsAtEnd(); ++IteratorA, ++IteratorB) {
    IteratorA.Set(((float)rand() / (float)RAND_MAX) - 0.5f);
    IteratorB.Set(((float)rand() / (float)RAND_MAX) - 0.5f);
  }

  // do Gaussian Smoothing for N iterations
  using FilterType = itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>;
  // FilterType::Pointer smoothingFilterA = FilterType::New();
  // FilterType::Pointer smoothingFilterB = FilterType::New();
  ImageType::Pointer tmpA = imageA;
  ImageType::Pointer tmpB = imageB;

  int numSmoothingSteps = iterations;
  for (int i = 0; i < numSmoothingSteps; i++) {
    FilterType::Pointer sFA = FilterType::New();
    sFA->SetSigma(smoothingKernelSize);
    sFA->SetInput(tmpA);
    sFA->Update();
    tmpA = sFA->GetOutput();

    // if we have a mask after every smoothing step we should tweak the volume
    if (mask != ITK_NULLPTR) {
      // assumption is that we have the same size, we can therefore iterate and adjust the border locations
      // of the label in the mask to be a zero-crossing. This has to be done in 3D coordinates... slow?
      // we should have a list of lists of voxel locations that are a border in the image. Based on how many
      // elements are in the list we can change shift the values towards zero. But what if two orthogonal
      // elements disagree? Can they? Because they share one voxel the cannot - that voxel has a fixed sign.
      // But the other two voxel can be larger than the seed voxel or smaller, we can average the combined effect
      // they have and adjust the nudge the intensities in that direction. - This is not geometry.... Isn't it
      // better to adjust the frequency spectrum of the data to start with? There we could filter in frequency
      // space and get two well defined frequencies.
    }
  }
  for (int i = 0; i < numSmoothingSteps; i++) {
    FilterType::Pointer sFB = FilterType::New();
    sFB->SetSigma(smoothingKernelSize);
    sFB->SetInput(tmpB);
    sFB->Update();
    tmpB = sFB->GetOutput();
  }

  // now compute the zero-crossing between the two volumes
  ImageType::Pointer erg = ImageType::New();
  erg->SetRegions(region);
  erg->Allocate();

  IteratorType IteratorE(erg, erg->GetLargestPossibleRegion());
  IteratorType itA(tmpA, tmpA->GetLargestPossibleRegion());
  IteratorType itB(tmpB, tmpB->GetLargestPossibleRegion());
  for (IteratorE.GoToBegin(), itA.GoToBegin(), itB.GoToBegin(); !itA.IsAtEnd() && !itB.IsAtEnd() && !IteratorE.IsAtEnd(); ++itA, ++itB, ++IteratorE) {
    if ((itA.Get() < (zero + threshold) && (itA.Get() > (zero - threshold))) && (itB.Get() < (zero + threshold) && (itB.Get() > (zero - threshold))))
      IteratorE.Set(4095.0);
    else
      IteratorE.Set(0.0);
  }

  // if we want to have void spaces we can create them here
  if (command.GetOptionWasSet("VoidSpaces")) {
    // use voidSpaces distance away and signs for placing void materials at intensity 1, 2, 3 and 4
    for (IteratorE.GoToBegin(), itA.GoToBegin(), itB.GoToBegin(); !itA.IsAtEnd() && !itB.IsAtEnd() && !IteratorE.IsAtEnd(); ++itA, ++itB, ++IteratorE) {
      if (fabs(itA.Get()) >= (threshold + voidSpaces) && fabs(itB.Get()) >= (threshold + voidSpaces)) {
        float testA = itA.Get();
        float testB = itB.Get();
        int type = 1; // both are negative
        if (testA > 0 && testB > 0)
          type = 2;
        else if (testA > 0 && testB < 0)
          type = 3;
        else if (testA < 0 && testB > 0)
          type = 4;
        IteratorE.Set(type);
      }
    }
    // add a lesion if we have to
    if (command.GetOptionWasSet("addLesion")) {
      int lesion_size = command.GetValueAsInt("addLesion", "addlesion");
      // size of the lesion should be in lesion_size
      if (lesion_size < 1) {
        fprintf(stderr, "ERROR: lesion size must be greater than 0. Set to 1 and continue.\n");
        lesion_size = 1;
      }

      // make a copy of erg with only the void space
      ImageType::Pointer ergVoidSpace = ImageType::New();
      ergVoidSpace->SetRegions(region);
      ergVoidSpace->Allocate();

      //
      IteratorType ierg(erg, erg->GetLargestPossibleRegion());
      IteratorType iergVoidSpace(ergVoidSpace, ergVoidSpace->GetLargestPossibleRegion());
      for (ierg.GoToBegin(), iergVoidSpace.GoToBegin(); !iergVoidSpace.IsAtEnd() && !ierg.IsAtEnd(); ++ierg, ++iergVoidSpace) {
        if (ierg.Get() == 1 || ierg.Get() == 2 || ierg.Get() == 3 || ierg.Get() == 4) {
          iergVoidSpace.Set(1);
        } else {
          iergVoidSpace.Set(0);
        }
      }
      // for this 0/1 volume we want to shrink it (to be able to ensure a lesion that is spherical)
      using StructuringElementType = itk::BinaryBallStructuringElement<OutputPixelType, 3>;
      using ErodeFilterType = itk::BinaryErodeImageFilter<ImageType, ImageType, StructuringElementType>;
      StructuringElementType structuringElement;
      structuringElement.SetRadius(1); // 3x3 structuring element
      structuringElement.CreateStructuringElement();
      ErodeFilterType::Pointer binaryErode = ErodeFilterType::New();
      binaryErode->SetKernel(structuringElement);
      binaryErode->SetInput(ergVoidSpace);
      binaryErode->SetErodeValue(lesion_size); // size of the lesion
      binaryErode->Update();
      ImageType::Pointer placeForLesion = binaryErode->GetOutput();

      // now we have to look for a random point in that volume (are there any voxel we can use?)
      IteratorType iplaceForLesion(placeForLesion, placeForLesion->GetLargestPossibleRegion());
      int validVoxel = 0;
      for (iplaceForLesion.GoToBegin(); !iplaceForLesion.IsAtEnd(); ++iplaceForLesion) {
        if (iplaceForLesion.Get() == 1) {
          validVoxel++;
        }
      }
      if (validVoxel < 1) {
        fprintf(stderr, "ERROR: no valid voxel found for a lesion.\n");
        return EXIT_FAILURE;
      }

      // create a random ellipsoid shape
      {
        using EllipseType = itk::EllipseSpatialObject<3>;
        ImageType::SizeType esize;
        esize[0] = lesion_size;
        esize[1] = lesion_size;
        esize[2] = lesion_size;
        SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
        imageFilter->SetSize(esize);
        ImageType::SpacingType spacing;
        spacing[0] = 1.0;
        spacing[1] = 1.0;
        spacing[2] = 1.0;
        imageFilter->SetSpacing(spacing);

        EllipseType::Pointer ellipse = EllipseType::New();
        EllipseType::ArrayType radiusArray;
        radiusArray[0] = lesion_size;
        radiusArray[1] = lesion_size / 2;
        radiusArray[1] = lesion_size / 2;
        // ellipse->SetRadiusInObjectSpace(  size[0] * 0.2 * spacing[0] );
        ellipse->SetRadiusInObjectSpace(radiusArray);

        ellipse->SetDefaultInsideValue(2048);
        ellipse->SetDefaultOutsideValue(0);

        using TransformType = EllipseType::TransformType;

        TransformType::Pointer transform = TransformType::New();

        transform->SetIdentity();

        TransformType::OutputVectorType translation;
        translation[0] = lesion_size / 2.0;
        translation[1] = lesion_size / 2.0;
        translation[2] = lesion_size / 2.0;
        transform->Translate(translation, false);

        ellipse->SetObjectToParentTransform(transform);
        imageFilter->SetInput(ellipse);
        imageFilter->SetUseObjectValue(true);
        imageFilter->SetOutsideValue(0);
      }
      // pick a location for the lesion
      unsigned int randomVoxel = rand() % validVoxel;

      // in the original volume put an object at this location
      int validVoxel = 0;
      for (iplaceForLesion.GoToBegin(), ierg.GoToBegin(); !iplaceForLesion.IsAtEnd() && !ierg.IsAtEnd(); ++iplaceForLesion, ++ierg) {
        if (iplaceForLesion.Get() == 1) {
          if (randomVoxel == validVoxel) {
            ierg.Set(2048);
            break;
          }
          validVoxel++;
        }
      }
    }
  }

  // we should blur the result?
  ImageType::Pointer result = ImageType::New();
  if (finalSmooth == 0) { // does this work for a float?
    result = erg;
  } else {
    FilterType::Pointer s = FilterType::New();
    s->SetSigma(finalSmooth);
    s->SetInput(erg);
    s->Update();
    result = s->GetOutput();
    // cleanup smoothing mess
    IteratorType iter(result, result->GetLargestPossibleRegion());
    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
      if (iter.Get() > 4068 || iter.Get() < 0)
        iter.Set(0);
    }
  }

  // cast to the output pixel type
  using OutputFilterType = itk::CastImageFilter<ImageType, OutputImageType>;
  OutputFilterType::Pointer filter = OutputFilterType::New();
  filter->SetInput(result);

  if (1) { // don't save this because its not with lungs separated
    typedef itk::ImageFileWriter<OutputImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();

    std::string volFileName = output;
    path p(volFileName);
    create_directories(p.parent_path());

    writer->SetFileName(volFileName);
    writer->SetInput(filter->GetOutput()); // not the right field

    std::cout << "Writing the Fake Lung Volume " << std::endl << std::endl;
    std::cout << volFileName << std::endl << std::endl;
    resultJSON["output"] = std::string(volFileName);

    try {
      writer->Update();
    } catch (itk::ExceptionObject &ex) {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::string res = resultJSON.dump(4) + "\n";
  std::ostringstream o;
  std::string si(output);
  si.erase(std::remove(si.begin(), si.end(), '\"'), si.end());
  o << output << "/" << si << ".json";
  std::ofstream out(o.str());
  out << res;
  out.close();

  fprintf(stdout, "%s", res.c_str());

  return EXIT_SUCCESS;
}