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
  command.SetDescription("FakeLungVolumes simulating CT contrast volumes. Exports volume files like nrrd or nifti based on the "
                         "provided file extension for the output image. The algorithm calculates the discrete intersection of two iso-surfaces "
                         "from band-pass filtered white noise volumes. Based on the amount of band-pass filtering and the threshold for the "
                         "intersection detection blood-vessel like pattern can be generated.");
  command.AddField("outfile", "Exported file name.", MetaCommand::STRING, true);

  command.SetOption("Resolution", "r", false, "Specify the resolution of the volume to be generated (in pixel as in 64x64x64).");
  command.AddOptionField("Resolution", "resolution", MetaCommand::STRING, false);

  command.SetOption("SmoothingKernelSize", "k", false, "Specify the kernel size for the Gaussian in pixel (7).");
  command.AddOptionField("SmoothingKernelSize", "kernelSize", MetaCommand::INT, false);

  command.SetOption("SmoothingIterations", "i", false, "Specify the number of times the Gaussian kernels are applied (2).");
  command.AddOptionField("SmoothingIterations", "iterations", MetaCommand::INT, false);

  command.SetOption("Threshold", "t", false, "Specify the threshold for zero-crossing (0.0001).");
  command.AddOptionField("Threshold", "threshold", MetaCommand::FLOAT, false);

  command.SetOption("finalSmooth", "f", false, "Specify the kernel size of a smoothing with a Gaussian at the end of the process (0).");
  command.AddOptionField("finalSmooth", "finalsmooth", MetaCommand::FLOAT, false);

  command.SetOption("Force", "f", false, "Ignore existing files and force overwrite.");

  command.SetOption("Verbose", "V", false, "Print more verbose output");

  if (!command.Parse(argc, argv)) {
    return 1;
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

  ImageType::Pointer imageA = ImageType::New();
  ImageType::Pointer imageB = ImageType::New();
  ImageType::IndexType start;
  start[0] =   0;  // first index on X 
  start[1] =   0;  // first index on Y
  start[2] =   0;  // first index on Z

  ImageType::SizeType  size;
  sscanf(resolution.c_str(), "%lux%lux%lu", &(size[0]), &(size[1]), &(size[2]));
  fprintf(stdout, "generate volume with: %lu %lu %lu voxel\n", size[0], size[1], size[2]);
  ImageType::RegionType region;
  
  region.SetSize( size );
  region.SetIndex( start );
  imageA->SetRegions( region );
  imageA->Allocate();
  imageB->SetRegions( region );
  imageB->Allocate();

  // set voxel values to random between -0.5 and 0.5
  using IteratorType = itk::ImageRegionIterator< ImageType >;
  IteratorType IteratorA( imageA, imageA->GetLargestPossibleRegion() ); 
  IteratorType IteratorB( imageB, imageB->GetLargestPossibleRegion() ); 
  srand(time(NULL));
  for ( IteratorA.GoToBegin(), IteratorB.GoToBegin(); !IteratorA.IsAtEnd() && !IteratorB.IsAtEnd(); ++IteratorA, ++IteratorB) {
    IteratorA.Set( ((float)rand() / (float)RAND_MAX) -0.5f );
    IteratorB.Set( ((float)rand() / (float)RAND_MAX) -0.5f );
  }

  // do Gaussian Smoothing for N iterations
  using FilterType = itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>;
  //FilterType::Pointer smoothingFilterA = FilterType::New();
  //FilterType::Pointer smoothingFilterB = FilterType::New();
  ImageType::Pointer tmpA = imageA;
  ImageType::Pointer tmpB = imageB;

  int numSmoothingSteps = iterations;
  for (int i = 0; i < numSmoothingSteps; i++) {
    FilterType::Pointer sFA = FilterType::New();
    sFA->SetSigma(smoothingKernelSize);
    sFA->SetInput(tmpA);
    sFA->Update();
    tmpA = sFA->GetOutput();
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
  erg->SetRegions( region );
  erg->Allocate();

  IteratorType IteratorE( erg, erg->GetLargestPossibleRegion() ); 
  IteratorType itA( tmpA, tmpA->GetLargestPossibleRegion() ); 
  IteratorType itB( tmpB, tmpB->GetLargestPossibleRegion() ); 
  for ( IteratorE.GoToBegin(), itA.GoToBegin(), itB.GoToBegin(); !itA.IsAtEnd() && !itB.IsAtEnd() && !IteratorE.IsAtEnd(); ++itA, ++itB, ++IteratorE) {
    if (fabs(itA.Get()) < threshold && fabs(itB.Get()) < threshold)
      IteratorE.Set( 4095.0 );
    else 
      IteratorE.Set( 0.0 );
  }

  // we should blur the result?
  ImageType::Pointer result = ImageType::New();
  if ( finalSmooth == 0 ) { // does this work for a float?
    result = erg;
  } else {
    FilterType::Pointer s = FilterType::New();
    s->SetSigma(finalSmooth);
    s->SetInput(erg);
    s->Update();
    result = s->GetOutput();
    // cleanup smoothing mess
    IteratorType iter( result, result->GetLargestPossibleRegion() ); 
    for ( iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
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