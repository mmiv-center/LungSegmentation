// Nice map of lung segments:
//    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5489231/

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

#include "itkGDCMImageIO.h"

#include "itkMetaDataDictionary.h"
#include "json.hpp"
#include "metaCommand.h"
#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <map>

#include "mytypes.h"

// this may need some types defined later, we should share those types in another header file
#include "curvedReslice/curvedReslice.hpp"

using json = nlohmann::json;
using namespace boost::filesystem;

// forward declaration
void CopyDictionary(itk::MetaDataDictionary &fromDict, itk::MetaDataDictionary &toDict);

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

template <typename TValue> TValue Convert(std::string optionString) {
  TValue value;
  std::istringstream iss(optionString);

  iss >> value;
  return value;
}

template <typename TValue> std::vector<TValue> ConvertVector(std::string optionString) {
  std::vector<TValue> values;
  std::string::size_type crosspos = optionString.find('x', 0);

  if (crosspos == std::string::npos) {
    values.push_back(Convert<TValue>(optionString));
  } else {
    std::string element = optionString.substr(0, crosspos);
    TValue value;
    std::istringstream iss(element);
    iss >> value;
    values.push_back(value);
    while (crosspos != std::string::npos) {
      std::string::size_type crossposfrom = crosspos;
      crosspos = optionString.find('x', crossposfrom + 1);
      if (crosspos == std::string::npos) {
        element = optionString.substr(crossposfrom + 1, optionString.length());
      } else {
        element = optionString.substr(crossposfrom + 1, crosspos);
      }
      std::istringstream iss2(element);
      iss2 >> value;
      values.push_back(value);
    }
  }
  return values;
}

json resultJSON;

int main(int argc, char *argv[]) {
  boost::posix_time::ptime timeLocal = boost::posix_time::second_clock::local_time();
  resultJSON["run_date_time"] = to_simple_string(timeLocal);

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(4);

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  command.SetDescription("LungSegmentation on CT DICOM images. Read DICOM image series and perform "
                         "a lung segmentation into left, right and trachea regions of interest. Exports a new DICOM series.");
  command.AddField("indir", "Directory with input DICOM image series.", MetaCommand::STRING, true);
  command.AddField("outdir", "Directory for output DICOM image series.", MetaCommand::STRING, true);

  // command.SetOption("Mask", "m", false, "Provide a mask image as an additional input. If not
  // supplied the mask will be calculated."); command.AddOptionField("Mask", "mask",
  // MetaCommand::STRING, false);

  // command.SetOption("Write", "s", false, "The shrink factor will make the problem easier to
  // handle (sub-sample data). The larger the value the faster."); command.AddOptionField("Write",
  // "shrinkFactor", MetaCommand::INT, false, "3");

  command.SetOption("SeriesName", "n", false, "Select series by series name (if more than one series is present).");
  command.AddOptionField("SeriesName", "seriesname", MetaCommand::STRING, false);

  command.SetOption("SaveLabelfield", "b", false, "Save the label field as a nifty file in the current directory");
  command.AddOptionField("SaveLabelfield", "labelfieldfilename", MetaCommand::STRING, true);

  command.SetOption("SaveNifty", "u", false, "Save the corrected dataset as a nifty image to the current directory");
  command.AddOptionField("SaveNifty", "niftyfilename", MetaCommand::STRING, true);

  command.SetOption("SaveReslice", "r", false, "Save a resliced version of the data as DICOM");
  // command.AddOptionField("SaveReslice", "reslicefilename", MetaCommand::STRING, true);

  command.SetOption("Force", "f", false, "Ignore existing directories and force reprocessing. Default is to stop processing if directory already exists.");

  command.SetOption("Verbose", "V", false, "Print more verbose output");

  if (!command.Parse(argc, argv)) {
    return 1;
  }

  std::string input = command.GetValueAsString("indir");
  std::string output = command.GetValueAsString("outdir");

  bool verbose = false;
  if (command.GetOptionWasSet("Verbose"))
    verbose = true;
  bool force = false;
  if (command.GetOptionWasSet("Force"))
    force = true;

  bool saveLabelField = false;
  bool saveNifty = false;
  bool seriesIdentifierFlag = false;

  if (command.GetOptionWasSet("SaveLabelfield"))
    saveLabelField = true;
  if (command.GetOptionWasSet("SaveNifty")) {
    saveNifty = true;
    fprintf(stdout, "will save nifty\n");
  }
  // std::string reslicefilename("");
  // if (command.GetOptionWasSet("SaveReslice")) {
  //  reslicefilename = command.GetValueAsString("SaveReslice", "reslicefilename");
  //}
  if (command.GetOptionWasSet("SeriesName"))
    seriesIdentifierFlag = true;
  std::string labelfieldfilename = command.GetValueAsString("SaveLabelfield", "labelfieldfilename");
  // todo: the argument could not be there, in this case the labelfieldfilename might be empty
  if (!saveLabelField) {
    labelfieldfilename = output + "/label_field.nii";
  }
  std::string niftyfilename = command.GetValueAsString("SaveNifty", "niftyfilename");
  std::string niftyfilename2 = niftyfilename + "_walls.nii";
  size_t lastdot = niftyfilename.find_last_of(".");
  if (lastdot == std::string::npos)
    niftyfilename2 = niftyfilename + "_walls.nii";
  else
    niftyfilename2 = niftyfilename.substr(0, lastdot) + "_walls.nii";

  std::string seriesName = command.GetValueAsString("SeriesName", "seriesname");

  // store information in the result json file
  resultJSON["command_line"] = json::array();
  for (int i = 0; i < argc; i++) {
    resultJSON["command_line"].push_back(std::string(argv[i]));
  }

  //  typedef signed short PixelType;
  typedef float FloatPixelType;
  const unsigned int Dimension = 3;

  //  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef itk::Image<FloatPixelType, Dimension> FloatImageType;

  typedef itk::ImageSeriesReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  typedef itk::Image<unsigned char, Dimension> MaskImageType;
  typedef itk::Image<unsigned char, 2> MaskSliceImageType;

  using StructuringElementType = itk::BinaryBallStructuringElement<PixelType, Dimension>;
  using ErodeFilterType = itk::BinaryErodeImageFilter<MaskImageType, MaskImageType, StructuringElementType>;
  using DilateFilterType = itk::BinaryDilateImageFilter<MaskImageType, MaskImageType, StructuringElementType>;

  typedef itk::ConnectedComponentImageFilter<MaskImageType, ImageType> ConnectedComponentImageFilterType;

  typedef itk::GDCMImageIO ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  dicomIO->LoadPrivateTagsOn();

  reader->SetImageIO(dicomIO);

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails(true);
  nameGenerator->AddSeriesRestriction("0008|0021");
  nameGenerator->SetRecursive(true);
  nameGenerator->SetDirectory(input);

  try {
    typedef std::vector<std::string> SeriesIdContainer;

    const SeriesIdContainer &seriesUID = nameGenerator->GetSeriesUIDs();

    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while (seriesItr != seriesEnd) {
      std::cout << "Found DICOM Series: ";
      std::cout << std::endl;
      std::cout << "  " << seriesItr->c_str() << std::endl;
      ++seriesItr;
    }

    std::string seriesIdentifier;

    SeriesIdContainer runThese;
    if (seriesIdentifierFlag) { // If no optional series identifier
      // seriesIdentifier = seriesName;
      runThese.push_back(seriesName);
    } else {
      // Todo: here we select only the first series. We should run
      // N3 on all series.

      seriesItr = seriesUID.begin();
      seriesEnd = seriesUID.end();
      // seriesIdentifier = seriesUID.begin()->c_str();
      while (seriesItr != seriesEnd) {
        runThese.push_back(seriesItr->c_str());
        ++seriesItr;
      }
      // Todo: if we have multiple phases they will all be in the same series.
      // It does not make sense to handle them here as one big volume, we should
      // look for the slice locations (consecutive) identify the first volume
      // and run N3 on that one. The resulting bias field should be applied to
      // all phases of the series.
    }

    seriesItr = runThese.begin();
    seriesEnd = runThese.end();
    while (seriesItr != seriesEnd) {
      seriesIdentifier = seriesItr->c_str();
      ++seriesItr;

      std::cout << "Processing series: " << std::endl;
      std::cout << "  " << seriesIdentifier << std::endl;

      std::string outputSeries = output + "/" + seriesIdentifier;
      if (!force && itksys::SystemTools::FileIsDirectory(outputSeries.c_str())) {
        fprintf(stdout, "Skip this series %s, output directory exists already...\n", outputSeries.c_str());
        exit(0); // this is no skip, that is giving up...
      }

      typedef std::vector<std::string> FileNamesContainer;
      FileNamesContainer fileNames;

      fileNames = nameGenerator->GetFileNames(seriesIdentifier);

      if (fileNames.size() < 10) {
        std::cout << "skip processing, not enough images in this series..." << std::endl;
        continue;
      }
      fprintf(stdout, "sufficient number of images (%lu) in this series\n", fileNames.size());
      resultJSON["series_identifier"] = seriesIdentifier;
      // for (int i = 0; i < fileNames.size(); i++) {
      //  resultJSON["file_names"].push_back(fileNames[i]);
      //}

      // here we read in all the slices as a single volume for processing
      // if we want to write them back out we have to read them slice by
      // slice and get a copy of the meta data for each slice
      reader->SetFileNames(fileNames);
      reader->ForceOrthogonalDirectionOff(); // do we need this?

      try {
        reader->Update();
      } catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
      }

      // read the data dictionary
      ImageType::Pointer inputImage = reader->GetOutput();
      typedef itk::MetaDataDictionary DictionaryType;
      DictionaryType &dictionary = dicomIO->GetMetaDataDictionary();
      fprintf(stdout, "pixel spacing of input is: %f %f %f\n", inputImage->GetSpacing()[0], inputImage->GetSpacing()[1], inputImage->GetSpacing()[2]);
      // overwrite a value in the dicom dictionary
      // std::string entryId( "0010|0010" );
      // std::string value( "MYNAME" );
      // itk::EncapsulateMetaData<std::string>( dictionary, entryId, value );

      std::string studyDescription;
      std::string seriesDescription;
      std::string patientID;
      std::string patientName;
      std::string sopClassUID;
      std::string seriesDate;
      std::string seriesTime;
      std::string studyDate;
      std::string studyTime;
      std::string patientSex;
      std::string convolutionKernelGroup;
      std::string modality;
      std::string manufacturer;
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|1030", studyDescription))
        resultJSON["SeriesDescription"] = boost::algorithm::trim_copy(seriesDescription);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|103e", seriesDescription))
        resultJSON["StudyDescription"] = boost::algorithm::trim_copy(studyDescription);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0016", sopClassUID))
        resultJSON["SOPClassUID"] = boost::algorithm::trim_copy(sopClassUID);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0021", seriesDate))
        resultJSON["StudyDate"] = boost::algorithm::trim_copy(studyDate);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0031", seriesTime))
        resultJSON["SeriesTime"] = boost::algorithm::trim_copy(seriesTime);
      if (itk::ExposeMetaData<std::string>(dictionary, "0010|0020", patientID))
        resultJSON["PatientID"] = boost::algorithm::trim_copy(patientID);
      if (itk::ExposeMetaData<std::string>(dictionary, "0010|0010", patientName))
        resultJSON["PatientName"] = boost::algorithm::trim_copy(patientName);
      if (itk::ExposeMetaData<std::string>(dictionary, "0010|0040", patientSex))
        resultJSON["PatientSex"] = boost::algorithm::trim_copy(patientSex);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0030", studyTime))
        resultJSON["StudyTime"] = boost::algorithm::trim_copy(studyTime);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0020", studyDate))
        resultJSON["SeriesDate"] = boost::algorithm::trim_copy(seriesDate);
      if (itk::ExposeMetaData<std::string>(dictionary, "0018|9316", convolutionKernelGroup))
        resultJSON["CTConvolutionKernelGroup"] = boost::algorithm::trim_copy(convolutionKernelGroup); // LUNG, BRAIN, BONE, SOFT_TISSUE
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0060", modality))
        resultJSON["Modality"] = boost::algorithm::trim_copy(modality);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0080", manufacturer))
        resultJSON["Manufacturer"] = boost::algorithm::trim_copy(manufacturer);

      MaskImageType::Pointer maskImage;

      // if we have a mask on the command line
      if (command.GetOptionWasSet("Mask")) {
        std::string maskName = command.GetValueAsString("Mask", "mask");
        typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
        MaskReaderType::Pointer maskreader = MaskReaderType::New();
        maskreader->SetFileName(maskName);
        try {
          maskreader->Update();
          maskImage = maskreader->GetOutput();
          maskImage->DisconnectPipeline();
        } catch (...) {
          maskImage = ITK_NULLPTR;
        }
      }

      if (!maskImage) { // segment tissue (body will be 1, air 0)
        typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> gaussianFilterType;
        gaussianFilterType::Pointer gaussianFilter = gaussianFilterType::New();
        gaussianFilter->SetInput(inputImage);
        gaussianFilter->SetVariance(2.0f);

        typedef itk::BinaryThresholdImageFilter<ImageType, MaskImageType> ThresholderType;
        ThresholderType::Pointer thresh = ThresholderType::New();
        thresh->SetInput(gaussianFilter->GetOutput());
        thresh->SetLowerThreshold(-20000);
        thresh->SetUpperThreshold(-600); // not too bright because otherwise the airways will connect with the lungs
        thresh->SetOutsideValue(1);      // this marks our lungs as empty space and the body as 1!
        thresh->SetInsideValue(0);

        thresh->Update();
        maskImage = thresh->GetOutput();
        maskImage->DisconnectPipeline();
        fprintf(stdout, "pixel spacing of mask image is: %f %f %f\n", maskImage->GetSpacing()[0], maskImage->GetSpacing()[1], maskImage->GetSpacing()[2]);
      }

      //
      // improve the mask for the body by hole filling and smoothing (grow + shrink)
      //
      DilateFilterType::Pointer binaryDilate = DilateFilterType::New(); // grows inside the tissue
      DilateFilterType::Pointer binaryDilate2 = DilateFilterType::New();
      ErodeFilterType::Pointer binaryErode = ErodeFilterType::New(); // grows inside the lungs
      ErodeFilterType::Pointer binaryErode2 = ErodeFilterType::New();
      ErodeFilterType::Pointer binaryErode3 = ErodeFilterType::New();
      ErodeFilterType::Pointer binaryErode4 = ErodeFilterType::New();

      StructuringElementType structuringElement;
      structuringElement.SetRadius(1); // 3x3 structuring element
      structuringElement.CreateStructuringElement();
      binaryDilate->SetKernel(structuringElement);
      binaryDilate2->SetKernel(structuringElement);
      binaryErode->SetKernel(structuringElement);
      binaryErode2->SetKernel(structuringElement);
      binaryErode3->SetKernel(structuringElement);
      binaryErode4->SetKernel(structuringElement);

      binaryDilate->SetInput(maskImage);
      binaryDilate2->SetInput(binaryDilate->GetOutput());
      binaryErode->SetInput(binaryDilate2->GetOutput());
      // binaryErode2->SetInput( binaryErode->GetOutput() );
      // binaryErode3->SetInput( binaryErode->GetOutput() );
      binaryErode4->SetInput(binaryErode->GetOutput());

      binaryDilate->SetDilateValue(1);
      binaryDilate2->SetDilateValue(1);
      binaryErode->SetErodeValue(1);
      binaryErode2->SetErodeValue(1);
      binaryErode3->SetErodeValue(1);
      binaryErode4->SetErodeValue(1);
      //        binaryDilate->Update();
      //        binaryDilate2->Update();
      //        binaryErode->Update();
      // binaryErode2->Update();
      // binaryErode3->Update();
      //        binaryErode4->Update();

      //
      // fill holes in the binary image slice by slice (fills in air inside the body)
      //
      typedef itk::SliceBySliceImageFilter<MaskImageType, MaskImageType> SliceFilterType;
      SliceFilterType::Pointer sliceFilter = SliceFilterType::New();
      sliceFilter->SetInput(maskImage /* binaryErode4->GetOutput() */); // smoothing should get us a better body mask
                                                                        // (touching bed problem)

      if (0) { // debug, save the body mask
        typedef itk::ImageFileWriter<MaskImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();
        // check if that directory exists, create before writing
        path p("/tmp/");
        writer->SetFileName("/tmp/body_mask.nii");
        writer->SetInput(/* binaryErode4->GetOutput() */ maskImage);

        std::cout << "Writing the initial mask image as " << std::endl;
        std::cout << "/tmp/body_mask.nii" << std::endl << std::endl;

        try {
          writer->Update();
        } catch (itk::ExceptionObject &ex) {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }
      }

      typedef itk::BinaryFillholeImageFilter<SliceFilterType::InternalInputImageType> HoleFillFilterType;
      HoleFillFilterType::Pointer holefillfilter = HoleFillFilterType::New();
      holefillfilter->SetForegroundValue(1);

      sliceFilter->SetFilter(holefillfilter);
      sliceFilter->Update();

      ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
      connected->SetBackgroundValue(0);
      connected->SetInput(sliceFilter->GetOutput());
      connected->Update();

      // keep only the large connected component
      typedef itk::LabelShapeKeepNObjectsImageFilter<ImageType> LabelShapeKeepNObjectsImageFilterType;
      LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter = LabelShapeKeepNObjectsImageFilterType::New();
      labelShapeKeepNObjectsImageFilter->SetInput(connected->GetOutput());
      labelShapeKeepNObjectsImageFilter->SetBackgroundValue(0);
      labelShapeKeepNObjectsImageFilter->SetNumberOfObjects(1);
      labelShapeKeepNObjectsImageFilter->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
      labelShapeKeepNObjectsImageFilter->Update();
      // the above image labelShapeKeepNObjectsImageFilter->GetOutput() is now a good mask for the
      // body of the participant next we need to segment air inside that mask (do we need to make
      // the mask tighter?, do we need to smooth the input? do we need to calculate this on a
      // smaller scale?) use binaryErode2->GetOutput(),
      // labelShapeKeepNObjectsImageFilter->GetOutput(), and inputImage to get lung mask

      MaskImageType::Pointer mask1 = maskImage;                                  // binaryErode4->GetOutput();
      ImageType::Pointer mask2 = labelShapeKeepNObjectsImageFilter->GetOutput(); // body mask
      // make mask2 appear over mask1
      mask2->SetOrigin(mask1->GetOrigin());
      mask2->SetSpacing(mask1->GetSpacing());
      mask2->SetDirection(mask1->GetDirection());

      MaskImageType::Pointer mask = MaskImageType::New();
      MaskImageType::RegionType maskRegion = inputImage->GetLargestPossibleRegion();
      MaskImageType::RegionType mask1Region = mask1->GetLargestPossibleRegion();
      ImageType::RegionType mask2Region = mask2->GetLargestPossibleRegion();

      mask->SetRegions(maskRegion);
      mask->Allocate();
      mask->FillBuffer(itk::NumericTraits<PixelType>::Zero);
      mask->SetOrigin(inputImage->GetOrigin());
      mask->SetSpacing(inputImage->GetSpacing());
      mask->SetDirection(inputImage->GetDirection());
      MaskImageType::SizeType regionSize = maskRegion.GetSize();
      itk::ImageRegionIterator<MaskImageType> maskIterator(mask, maskRegion);
      itk::ImageRegionIterator<MaskImageType> inputMask1Iterator(mask1, mask1Region);
      itk::ImageRegionIterator<ImageType> inputMask2Iterator(mask2, mask2Region);
      // everything that is 1 in mask2 and 0 in mask1
      // problem is that mask2 can also be a value larger than 1! example is 0263 where the value is
      // == 2
      while (!maskIterator.IsAtEnd() && !inputMask1Iterator.IsAtEnd() && !inputMask2Iterator.IsAtEnd()) {
        if (inputMask2Iterator.Get() > 0 && inputMask1Iterator.Get() == 0) {
          // if (maskIterator.GetIndex()[0] > static_cast<ImageType::IndexValueType>(regionSize[0])
          // / 2) {
          maskIterator.Set(1);
        }
        ++maskIterator;
        ++inputMask1Iterator;
        ++inputMask2Iterator;
      }
      fprintf(stdout, "pixel spacing of grow/shrunk image is: %f %f %f\n", mask->GetSpacing()[0], mask->GetSpacing()[1], mask->GetSpacing()[2]);
      if (0) { // debug, save the body mask
        typedef itk::ImageFileWriter<MaskImageType> WriterType;
        WriterType::Pointer writer1 = WriterType::New();
        // check if that directory exists, create before writing
        writer1->SetFileName("/tmp/bodyLung_mask01.nii");
        writer1->SetInput(mask1 /* mask */);

        std::cout << "Writing the lung mask as " << std::endl;
        std::cout << "/tmp/bodyLung_mask01.nii" << std::endl << std::endl;

        try {
          writer1->Update();
        } catch (itk::ExceptionObject &ex) {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }
        typedef itk::ImageFileWriter<ImageType> WriterType2;
        WriterType2::Pointer writer2 = WriterType2::New();
        // check if that directory exists, create before writing
        writer2->SetFileName("/tmp/bodyLung_mask02.nii");
        writer2->SetInput(mask2 /* mask */);

        std::cout << "Writing the lung mask as " << std::endl;
        std::cout << "/tmp/bodyLung_mask02.nii" << std::endl << std::endl;

        try {
          writer2->Update();
        } catch (itk::ExceptionObject &ex) {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }

        // connected->GetOutput()
        WriterType2::Pointer writer3 = WriterType2::New();
        // check if that directory exists, create before writing
        writer3->SetFileName("/tmp/bodyLung_mask03.nii");
        writer3->SetInput(connected->GetOutput());

        std::cout << "Writing the lung mask as " << std::endl;
        std::cout << "/tmp/bodyLung_mask03.nii" << std::endl << std::endl;

        try {
          writer3->Update();
        } catch (itk::ExceptionObject &ex) {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }

        WriterType::Pointer writer4 = WriterType::New();
        // check if that directory exists, create before writing
        writer4->SetFileName("/tmp/bodyLung_mask04.nii");
        writer4->SetInput(mask);

        std::cout << "Writing the lung mask as " << std::endl;
        std::cout << "/tmp/bodyLung_mask04.nii" << std::endl << std::endl;

        try {
          writer4->Update();
        } catch (itk::ExceptionObject &ex) {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }
      }

      // now in mask we have every air filled cavity inside the body, look for the largest one,
      // assume its the lungs this will remove the outside the body parts that are created because
      // the slice filling adds them in 2d
      ConnectedComponentImageFilterType::Pointer connected2 = ConnectedComponentImageFilterType::New();
      connected2->SetBackgroundValue(0);
      connected2->SetInput(mask);
      connected2->Update();
      // std::cout << "Number of connected components: " << connected->GetObjectCount() <<
      // std::endl;

      // std::cout << "Number of connected components: " << connected->GetObjectCount() <<
      // std::endl; we sometimes get the two lungs separated here, we should look for the size of
      // connected components we do have a volume range we are looking for and getting two lungs
      // with different volumes would be possible (see 0183)

      using LabelType = unsigned short;
      using ShapeLabelObjectType = itk::ShapeLabelObject<LabelType, 3>;
      using LabelMapType = itk::LabelMap<ShapeLabelObjectType>;
      using labelType = itk::LabelImageToShapeLabelMapFilter<ImageType, LabelMapType>;
      labelType::Pointer label = labelType::New();
      label->SetInput(connected2->GetOutput());
      label->SetComputePerimeter(true);
      label->Update();

      int useNObjects = 0;
      LabelMapType *labelMap = label->GetOutput();
      if (labelMap->GetNumberOfLabelObjects() == 0) {
        // error case
        fprintf(stderr, "Could not find any object in the data using the current set of "
                        "thresholds, we can try to start again by lowering the threshold?\n");
      }

      for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n) {
        ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
        // fprintf(stdout, "object %d has %0.4f liters, %lu voxel\n", n,
        // labelObject->GetPhysicalSize() / 1000000, labelObject->GetNumberOfPixels());

        if (labelObject->GetNumberOfPixels() > 150000) // magick number using sample of 1 - should be selected based on volume instead
                                                       // of number of pixel
          useNObjects++;
      }
      if (useNObjects < 1) {
        fprintf(stdout, "useNObjects is: %d. Set manually to 1.\n", useNObjects);
        useNObjects = 1;
      }
      fprintf(stdout, "found %d object%s with number of voxel large enough to be of interest...\n", useNObjects, useNObjects == 1 ? "" : "s");

      // keep only the large connected component
      typedef itk::LabelShapeKeepNObjectsImageFilter<ImageType> LabelShapeKeepNObjectsImageMaskFilterType;
      LabelShapeKeepNObjectsImageMaskFilterType::Pointer labelShapeKeepNObjectsImageFilter2 = LabelShapeKeepNObjectsImageMaskFilterType::New();
      labelShapeKeepNObjectsImageFilter2->SetInput(connected2->GetOutput());
      labelShapeKeepNObjectsImageFilter2->SetBackgroundValue(0);
      labelShapeKeepNObjectsImageFilter2->SetNumberOfObjects(useNObjects);
      labelShapeKeepNObjectsImageFilter2->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
      labelShapeKeepNObjectsImageFilter2->Update();

      // the label filter will keep the id of the label, its not 1, what is it?
      /*typedef itk::MinimumMaximumImageCalculator <ImageType> ImageCalculatorFilterType;
        ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New
        (); imageCalculatorFilter->SetImage( labelShapeKeepNObjectsImageFilter2->GetOutput() );
        imageCalculatorFilter->ComputeMaximum();
        int maxValue = imageCalculatorFilter->GetMaximum();*/
      int maxValue = 1;

      // instead of keeping the labels we should get a single label here with the value 1
      ImageType::Pointer lungs = labelShapeKeepNObjectsImageFilter2->GetOutput();
      ImageType::RegionType lungRegion = lungs->GetLargestPossibleRegion();
      itk::ImageRegionIterator<ImageType> lungIterator(lungs, lungRegion);
      while (!lungIterator.IsAtEnd()) {
        if (lungIterator.Get() > 0) {
          lungIterator.Set(1);
        }
        ++lungIterator;
      }

      //
      // and fill holes in the final segmentation of the lung
      //
      typedef itk::SliceBySliceImageFilter<ImageType, ImageType> SliceFilterImageType;
      SliceFilterImageType::Pointer sliceFilter2 = SliceFilterImageType::New();
      sliceFilter2->SetInput(labelShapeKeepNObjectsImageFilter2->GetOutput());

      typedef itk::BinaryFillholeImageFilter<SliceFilterImageType::InternalInputImageType> HoleFillFilterType2;
      HoleFillFilterType2::Pointer holefillfilter2 = HoleFillFilterType2::New();
      holefillfilter2->SetForegroundValue(maxValue);

      sliceFilter2->SetFilter(holefillfilter2);
      sliceFilter2->Update();

      // now apply the lung mask to the input image and export that instead
      ImageType::Pointer final = sliceFilter2->GetOutput();
      ImageType::RegionType imRegion = inputImage->GetLargestPossibleRegion();
      ImageType::RegionType maskRegion2 = final->GetLargestPossibleRegion();
      itk::ImageRegionIterator<ImageType> imIterator(inputImage, imRegion);
      itk::ImageRegionIterator<ImageType> maskIterator2(final, maskRegion2);

      MaskImageType::Pointer labelField = MaskImageType::New();
      MaskImageType::RegionType labelFieldRegion = inputImage->GetLargestPossibleRegion();
      labelField->SetRegions(labelFieldRegion);
      labelField->Allocate();
      labelField->SetOrigin(inputImage->GetOrigin());
      labelField->SetSpacing(inputImage->GetSpacing());
      labelField->SetDirection(inputImage->GetDirection());
      itk::ImageRegionIterator<MaskImageType> labelFieldIterator(labelField, labelFieldRegion);

      // as the background voxel we can use the first voxel overall
      int firstVoxel;
      bool gotFirstVoxel = false;
      size_t numberOfPixel = 0;
      while (!imIterator.IsAtEnd() && !maskIterator2.IsAtEnd()) {
        if (!gotFirstVoxel) {
          gotFirstVoxel = true;
          firstVoxel = imIterator.Get();
        }

        if (maskIterator2.Get() != 0) {
          // if (maskIterator.GetIndex()[0] > static_cast<ImageType::IndexValueType>(regionSize[0])
          // / 2) {
          maskIterator2.Set(imIterator.Get());
          labelFieldIterator.Set(1);
          numberOfPixel++;
        } else {
          maskIterator2.Set(firstVoxel); // this depends on the pixel representation its ok
          labelFieldIterator.Set(0);
        }
        ++maskIterator2;
        ++imIterator;
        ++labelFieldIterator;
      }
      resultJSON["lung_number_of_voxel"] = numberOfPixel;
      // calculate the volume of the lungs
      double spacingx = inputImage->GetSpacing()[0];
      double spacingy = inputImage->GetSpacing()[1];
      double spacingz = inputImage->GetSpacing()[2];
      double originx = inputImage->GetOrigin()[0];
      double originy = inputImage->GetOrigin()[1];
      double originz = inputImage->GetOrigin()[2];
      resultJSON["voxel_size"] = json::array();
      resultJSON["voxel_size"].push_back(spacingx);
      resultJSON["voxel_size"].push_back(spacingy);
      resultJSON["voxel_size"].push_back(spacingz);
      resultJSON["data_dims"] = json::array();
      resultJSON["data_dims"].push_back(inputImage->GetLargestPossibleRegion().GetSize()[0]);
      resultJSON["data_dims"].push_back(inputImage->GetLargestPossibleRegion().GetSize()[1]);
      resultJSON["data_dims"].push_back(inputImage->GetLargestPossibleRegion().GetSize()[2]);

      // in liters
      resultJSON["lung_volume"] = numberOfPixel * spacingx * spacingy * spacingz * 1e-6;

      // save the label field (one lung)
      if (saveLabelField) {
        typedef itk::ImageFileWriter<MaskImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();
        // check if that directory exists, create before writing
        path p(labelfieldfilename);
        create_directories(p.parent_path());
        writer->SetFileName(labelfieldfilename);
        writer->SetInput(labelField);

        std::cout << "Writing the label field as " << std::endl;
        std::cout << labelfieldfilename << std::endl << std::endl;
        resultJSON["label_field_filename"] = std::string(labelfieldfilename);
        fprintf(stdout, "label field voxel size is: %f %f %f\n", labelField->GetSpacing()[0], labelField->GetSpacing()[1], labelField->GetSpacing()[2]);
        // std::string res2 = resultJSON.dump(4) + "\n";
        // fprintf(stdout, "%s", res2.c_str());

        try {
          writer->Update();
        } catch (itk::ExceptionObject &ex) {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }
      }

      // in some slices we have more information, for example if there are only two objects in a
      // slice we can assume that we have the left and the right lung. If there are three objects
      // the middle one would be the trachae.

      // we can also just give up and start computing a lung average image after elastic
      // registration that one could be segmented and we take the label from there (if air and over
      // average, use label from average).

      ////////////////////////////////////////////////////////////////////////////////////////////////
      // we should now split the lungs into the airways and the lungs using region growing
      // for this we need to find out if we have a label in the top image... we should have 3
      // regions of interest if we start from the top, we should set a seed point into the middle
      // one
      using ExtractFilterType = itk::ExtractImageFilter<MaskImageType, MaskImageType>;
      ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
      extractFilter->SetDirectionCollapseToSubmatrix();
      const MaskImageType *inImage = labelField;
      ImageType::RegionType inputRegion = inImage->GetLargestPossibleRegion();
      // GetBufferedRegion();
      ImageType::SizeType size = inputRegion.GetSize();
      unsigned int lastSlice = size[2];
      unsigned int sliceNumber = size[2] - 1;
      int minRegion = 0;
      float minSeedx = 0;
      float minSeedy = 0;
      int minSeedIdx = 0;
      int minSeedIdy = 0;
      int minDist = 0;
      int sliceSeedStart = 0;
      int minNumObjects = 0;
      for (; sliceNumber >= 0; sliceNumber--) {
        size[2] = 1; // we extract along z direction
        ImageType::IndexType start = inputRegion.GetIndex();
        start[2] = sliceNumber;
        ImageType::RegionType desiredRegion;
        desiredRegion.SetSize(size);
        desiredRegion.SetIndex(start);
        extractFilter->SetExtractionRegion(desiredRegion);
        extractFilter->SetInput(inImage);

        // Now we can run connected components on the 2D image of the first slice
        // we expect 3 distinct regions
        ConnectedComponentImageFilterType::Pointer connected3 = ConnectedComponentImageFilterType::New();
        connected3->SetBackgroundValue(0);
        connected3->SetInput(extractFilter->GetOutput());
        connected3->Update();

        // using LabelType = unsigned short;
        // using ShapeLabelObjectType = itk::ShapeLabelObject< LabelType, 3 >;
        // using LabelMapType = itk::LabelMap< ShapeLabelObjectType >;
        // using labelType = itk::LabelImageToShapeLabelMapFilter< ImageType, LabelMapType>;
        labelType::Pointer labelTmp = labelType::New();
        labelTmp->SetInput(connected3->GetOutput());
        labelTmp->SetComputePerimeter(true);
        labelTmp->Update();

        // int useNObjects = 0;
        LabelMapType *labelMapTmp = labelTmp->GetOutput();
        if (labelMapTmp->GetNumberOfLabelObjects() == 0) {
          // error case
          fprintf(stderr,
                  "Look at slice %d, Could not find any object in the data using the current set "
                  "of thresholds, we can try to start again by "
                  "lowering the threshold?\n",
                  sliceNumber);
        }

        for (unsigned int n = 0; n < labelMapTmp->GetNumberOfLabelObjects(); ++n) {
          ShapeLabelObjectType *labelObject = labelMapTmp->GetNthLabelObject(n);
          fprintf(stdout, "top slice: object %d has %0.4f liters, %lu voxel\n", n, labelObject->GetPhysicalSize() / 1000000, labelObject->GetNumberOfPixels());

          // if (labelObject->GetNumberOfPixels() > 150000) // magick number using sample of 1 -
          // should be selected based on volume instead of number of pixel
          //  useNObjects++;
        }
        // it might be better to be able to stop early, there are some series with only a single
        // lung for those we have to stop at 2. What we could do is stop as well if we are too far
        // down in the image stack
        if (labelMapTmp->GetNumberOfLabelObjects() == 3 || sliceNumber < 2 * lastSlice / 3) {
          sliceSeedStart = sliceNumber;
          fprintf(stdout, "found %lu objects in slice %d of %d\n", labelMapTmp->GetNumberOfLabelObjects(), sliceNumber + 1, lastSlice);
          // seed point would be in the object that is closest to the center of mass in the image
          // The above assumption is sometimes wrong. One of the lungs might be closer to the center of the image.
          // Instead we could assume that the lungs are oriented in the same way in all images. That allows us to
          // look at a single coordinate and to select the object that is in the middle of the other two objects.

          // we could go to the last slice that has 3 objects as well...
          // what if we have only one lung? We should detect the diameter of the airways instead..
          // something like at least two objects

          double spacingx = inputImage->GetSpacing()[0];
          double spacingy = inputImage->GetSpacing()[1];
          double spacingz = inputImage->GetSpacing()[2];
          double originx = inputImage->GetOrigin()[0];
          double originy = inputImage->GetOrigin()[1];
          double originz = inputImage->GetOrigin()[2];

          // what is in the middle in the x-direction? (assumes we have three objects!)
          double x1 = labelMapTmp->GetNthLabelObject(0)->GetCentroid()[0];
          double x2 = labelMapTmp->GetNthLabelObject(1)->GetCentroid()[0];
          double x3 = labelMapTmp->GetNthLabelObject(2)->GetCentroid()[0];
          if ((x1 > x2 && x1 < x3) || (x1 > x3 && x1 < x2)) {
            // point 0 is the one we are looking for
            minDist = 0;
            minRegion = 0;
            double x = labelMapTmp->GetNthLabelObject(0)->GetCentroid()[0];
            double y = labelMapTmp->GetNthLabelObject(0)->GetCentroid()[1];
            int seedx = roundf((x - originx) / spacingx);
            int seedy = roundf((y - originy) / spacingy);
            minSeedx = seedx;
            minSeedy = seedy;
            minNumObjects = labelMapTmp->GetNumberOfLabelObjects();
          } else if ((x2 > x1 && x2 < x3) || (x2 > x3 && x2 < x1)) {
            // point 1 is the one we are looking for
            minDist = 0;
            minRegion = 0;
            double x = labelMapTmp->GetNthLabelObject(1)->GetCentroid()[0];
            double y = labelMapTmp->GetNthLabelObject(1)->GetCentroid()[1];
            int seedx = roundf((x - originx) / spacingx);
            int seedy = roundf((y - originy) / spacingy);
            minSeedx = seedx;
            minSeedy = seedy;
            minNumObjects = labelMapTmp->GetNumberOfLabelObjects();
          } else {
            // point 3 is the one we are looking for
            minDist = 0;
            minRegion = 0;
            double x = labelMapTmp->GetNthLabelObject(2)->GetCentroid()[0];
            double y = labelMapTmp->GetNthLabelObject(2)->GetCentroid()[1];
            int seedx = roundf((x - originx) / spacingx);
            int seedy = roundf((y - originy) / spacingy);
            minSeedx = seedx;
            minSeedy = seedy;
            minNumObjects = labelMapTmp->GetNumberOfLabelObjects();
          }

          // which one is closest to the center (center of mass)
          /*for (unsigned int n = 0; n < labelMapTmp->GetNumberOfLabelObjects(); n++) {
            ShapeLabelObjectType *labelObject = labelMapTmp->GetNthLabelObject(n);
            double x = labelObject->GetCentroid()[0];
            double y = labelObject->GetCentroid()[1];
            // fprintf(stdout, "location of object in pixel %d is %f x %f\n", n, x, y);
            MaskImageType::RegionType bb = inputRegion;
            double spacingx = inputImage->GetSpacing()[0];
            double spacingy = inputImage->GetSpacing()[1];
            double spacingz = inputImage->GetSpacing()[2];
            double originx = inputImage->GetOrigin()[0];
            double originy = inputImage->GetOrigin()[1];
            double originz = inputImage->GetOrigin()[2];
            x = x;
            // originx + x *(spacingx);
            y = y;
            // originy + y *(spacingy);
            fprintf(stdout, "location of object in patient space %d is %f x %f\n", n, x, y);

            double COMx = originx + (spacingx) * ((size[0] - 1) / 2.0);
            double COMy = originy + (spacingy) * ((size[1] - 1) / 2.0);
            int seedx = roundf((x - originx) / spacingx);
            int seedy = roundf((y - originy) / spacingy);

            fprintf(stdout, "%lu %lu %lu, %f %f %f, %f %f %f, center %f %f\n", size[0], size[1], size[2], originx, originy, originz, spacingx, spacingy,
                    spacingz, COMx, COMy);

            double dist = sqrtf((COMx - x) * (COMx - x) + (COMy - y) * (COMy - y));
            if (n == 0 || dist < minDist) {
              minDist = dist;
              minRegion = n;
              minSeedx = seedx;
              minSeedy = seedy;
              minNumObjects = labelMapTmp->GetNumberOfLabelObjects();
            }
            fprintf(stdout, "object %d is %f far away from center, centroid is at: %d %d\n", n, dist, seedx, seedy);

            // inputRegion
          } */
          break;
        }
      }
      fprintf(stdout, "min object is: %d, seed is at %f %f slice %d\n", minRegion, minSeedx, minSeedy, sliceSeedStart);
      resultJSON["trachea_slice_location_pixel"] = json::array();
      resultJSON["trachea_slice_location_pixel"].push_back(minSeedx);
      resultJSON["trachea_slice_location_pixel"].push_back(minSeedy);
      resultJSON["trachea_slice_location_pixel"].push_back(sliceSeedStart);
      resultJSON["trachea_slice_number_objects"] = minNumObjects;
      resultJSON["trachea_slice_min_region"] = minRegion;
      resultJSON["trachea_slice_min_distance"] = minDist;

      // before we try to separate the two lungs we should shrink our label again,
      // that should help in cases where the trachea touch the lung wall
      // if we do this we also have to grow the region again before we do an assignment to tissue
      // type

      ///////////////////////////////////////////////////////////////////////////////
      // ok now do a region growing in labelField starting at minSeedx maxSeedy for all voxel that
      // have the value 1 we should have a second label field here how do we find neighbors?
      MaskImageType::Pointer regionGrowingField = MaskImageType::New();
      MaskImageType::RegionType regionGrowingFieldRegion = inputImage->GetLargestPossibleRegion();
      regionGrowingField->SetRegions(regionGrowingFieldRegion);
      regionGrowingField->Allocate();
      regionGrowingField->FillBuffer(itk::NumericTraits<PixelType>::Zero);
      regionGrowingField->SetOrigin(inputImage->GetOrigin());
      regionGrowingField->SetSpacing(inputImage->GetSpacing());
      regionGrowingField->SetDirection(inputImage->GetDirection());
      // itk::ImageRegionIterator<MaskImageType> regionGrowingFieldIterator(regionGrowingField,
      // regionGrowingFieldRegion); how to we find neighboring voxel from current location minSeedx
      // and minSeedy and sliceSeedStart?
      using PointType = itk::Point<int, MaskImageType::ImageDimension>;
      std::vector<PointType> front;
      PointType p1;
      p1[0] = minSeedx;
      p1[1] = minSeedy;
      p1[2] = sliceSeedStart;
      front.push_back(p1);
      // we should have a 1 at this place
      ImageType::IndexType pixelIndex = {{p1[0], p1[1], p1[2]}};
      regionGrowingField->SetPixel(pixelIndex, 1);
      unsigned int count = 0;
      unsigned int frontSize = 0;
      size = regionGrowingFieldRegion.GetSize();
      int burnInSize = 50000 * (1.0 / (inputImage->GetSpacing()[0] * inputImage->GetSpacing()[1] * inputImage->GetSpacing()[2])) / (1.0 / (0.84 * 0.84 * 0.8));
      resultJSON["trachea_burn_in_size"] = burnInSize;
      fprintf(stdout, "size of output: %lu %lu %lu\n", size[0], size[1], size[2]);
      while (1) {
        if (front.size() == 0)
          break;

        count++;
        // fprintf(stdout, "frontSize is %d\n", front.size());
        // how do we check the first value from front in labelField? We want to know if the value is
        // 0 fprintf(stderr, "size of front before step %d: %lu\n", count, front.size());
        PointType p;
        p[0] = front[0][0];
        p[1] = front[0][1];
        p[2] = front[0][2];
        front.erase(front.begin());
        // fprintf(stderr, "remove one front pixel in step %d. Size of front is now: %lu\n", count,
        // front.size()); pixelIndex = {{p[0], p[1], p[2]}}; fprintf(stdout, "read pixel value at
        // location %d %d %d as : %d\n", p[0], p[1], p[2], inImage->GetPixel(pixelIndex));

        // get the 6 neighbors
        PointType pp1;
        pp1[0] = p[0] + 1;
        pp1[1] = p[1];
        pp1[2] = p[2];
        // inside boundaries?
        if (pp1[0] < size[0]) {
          // are we still in the region of interest?
          pixelIndex = {{pp1[0], pp1[1], pp1[2]}};
          if (inImage->GetPixel(pixelIndex) == 1 && regionGrowingField->GetPixel(pixelIndex) == 0) {
            // add to front as a new pixel
            regionGrowingField->SetPixel(pixelIndex, 1);
            front.push_back(pp1);
            // fprintf(stderr, "add pp1\n");
          }
        }

        PointType pp2;
        pp2[0] = p[0] - 1;
        pp2[1] = p[1];
        pp2[2] = p[2];
        // inside boundaries?
        if (pp2[0] >= 0) {
          // are we still in the region of interest?
          pixelIndex = {{pp2[0], pp2[1], pp2[2]}};
          if (inImage->GetPixel(pixelIndex) == 1 && regionGrowingField->GetPixel(pixelIndex) == 0) {
            // add to front as a new pixel
            regionGrowingField->SetPixel(pixelIndex, 1);
            front.push_back(pp2);
            // fprintf(stderr, "add pp2\n");
          }
        }

        PointType pp3;
        pp3[0] = p[0];
        pp3[1] = p[1] + 1;
        pp3[2] = p[2];
        // inside boundaries?
        if (pp3[1] < size[1]) {
          // are we still in the region of interest?
          pixelIndex = {{pp3[0], pp3[1], pp3[2]}};
          if (inImage->GetPixel(pixelIndex) == 1 && regionGrowingField->GetPixel(pixelIndex) == 0) {
            // add to front as a new pixel
            regionGrowingField->SetPixel(pixelIndex, 1);
            front.push_back(pp3);
            // fprintf(stderr, "add pp3\n");
          }
        }

        PointType pp4;
        pp4[0] = p[0];
        pp4[1] = p[1] - 1;
        pp4[2] = p[2];
        // inside boundaries?
        if (pp4[1] >= 0) {
          // are we still in the region of interest?
          pixelIndex = {{pp4[0], pp4[1], pp4[2]}};
          if (inImage->GetPixel(pixelIndex) == 1 && regionGrowingField->GetPixel(pixelIndex) == 0) {
            // add to front as a new pixel
            regionGrowingField->SetPixel(pixelIndex, 1);
            front.push_back(pp4);
            // fprintf(stderr, "add pp4\n");
          }
        }

        PointType pp5;
        pp5[0] = p[0];
        pp5[1] = p[1];
        pp5[2] = p[2] + 1;
        // inside boundaries?
        if (pp5[2] < size[2]) {
          // are we still in the region of interest?
          pixelIndex = {{pp5[0], pp5[1], pp5[2]}};
          if (inImage->GetPixel(pixelIndex) == 1 && regionGrowingField->GetPixel(pixelIndex) == 0) {
            // add to front as a new pixel
            regionGrowingField->SetPixel(pixelIndex, 1);
            front.push_back(pp5);
            // fprintf(stderr, "add pp5\n");
          }
        }

        PointType pp6;
        pp6[0] = p[0];
        pp6[1] = p[1];
        pp6[2] = p[2] - 1;
        // inside boundaries?
        if (pp6[2] >= 0) {
          // are we still in the region of interest?
          pixelIndex = {{pp6[0], pp6[1], pp6[2]}};
          if (inImage->GetPixel(pixelIndex) == 1 && regionGrowingField->GetPixel(pixelIndex) == 0) {
            // add to front as a new pixel
            regionGrowingField->SetPixel(pixelIndex, 1);
            front.push_back(pp6);
            // fprintf(stderr, "add pp6\n");
          }
        }
        // how long is the list in front? This burn in phase might not be required,
        // if we can add all pixel of the initial trachea into the region of interest
        // problem here is that burn in should be different is we have larger voxel sizes
        // we need a shorter burn in if we have voxel that are three times larger
        if (count == burnInSize) {
          frontSize = front.size();
          fprintf(stdout, "frontSize after burn in size of %d is %d\n", burnInSize, frontSize);
        }
        if (count > burnInSize) {
          // fprintf(stdout, "iteration on frontSize: %lu\n", front.size());
          if (front.size() > frontSize * 2) {
            // found the last entry
            fprintf(stdout, "stop at count == %d\n", count);
            resultJSON["trachea_stop_voxel_count"] = count;
            resultJSON["trachea_stop_voxel_front_size"] = front.size();
            break; // don't continue with region growing, show the output
          }
        }
      }
      // can we save the air ways right now?
      if (saveLabelField) {
        typedef itk::ImageFileWriter<MaskImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();
        std::string a;
        size_t lastdot = labelfieldfilename.find_last_of(".");
        if (lastdot == std::string::npos)
          a = labelfieldfilename + "_trachea.nii";
        else
          a = labelfieldfilename.substr(0, lastdot) + "_trachea.nii";

        // make sure this directory exists
        path p(labelfieldfilename);
        create_directories(p.parent_path());

        // std::string a(labelfieldfilename + "trachea.nii");
        writer->SetFileName(a);
        writer->SetInput(regionGrowingField);

        std::cout << "Writing the trachea field as " << std::endl;
        std::cout << a << std::endl << std::endl;
        resultJSON["trachea_filename"] = a;

        try {
          writer->Update();
        } catch (itk::ExceptionObject &ex) {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }
      }

      // after separating the airways out we should check again if we have two large regions of
      // interest left those would be left and right lung respectively. Each one should get its own
      // label. new Mask for the lung (not in trachea)
      MaskImageType::Pointer no_trachea = MaskImageType::New();
      MaskImageType::RegionType no_tracheaRegion = inputImage->GetLargestPossibleRegion();
      no_trachea->SetRegions(no_tracheaRegion);
      no_trachea->Allocate();
      no_trachea->SetOrigin(inputImage->GetOrigin());
      no_trachea->SetSpacing(inputImage->GetSpacing());
      no_trachea->SetDirection(inputImage->GetDirection());
      itk::ImageRegionIterator<MaskImageType> no_tracheaIterator(no_trachea, no_tracheaRegion);
      itk::ImageRegionIterator<MaskImageType> tracheaIterator(regionGrowingField, no_tracheaRegion);
      itk::ImageRegionIterator<MaskImageType> lungMaskIterator(labelField, no_tracheaRegion);
      while (!lungMaskIterator.IsAtEnd() && !tracheaIterator.IsAtEnd() && !no_tracheaIterator.IsAtEnd()) {
        if (lungMaskIterator.Get() > 0 && tracheaIterator.Get() != 1) {
          no_tracheaIterator.Set(1);
        } else {
          no_tracheaIterator.Set(0);
        }
        ++lungMaskIterator;
        ++tracheaIterator;
        ++no_tracheaIterator;
      }
      // can we save the lungs without trachea now?
      if (saveLabelField) {
        typedef itk::ImageFileWriter<MaskImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();
        std::string a;
        size_t lastdot = labelfieldfilename.find_last_of(".");
        if (lastdot == std::string::npos)
          a = labelfieldfilename + "_lung_no_trachea.nii";
        else
          a = labelfieldfilename.substr(0, lastdot) + "_lung_no_trachea.nii";

        path p(labelfieldfilename);
        create_directories(p.parent_path());

        // std::string a(labelfieldfilename + "trachea.nii");
        writer->SetFileName(a);
        writer->SetInput(no_trachea);

        std::cout << "Writing the no trachea lung mask as " << std::endl;
        std::cout << a << std::endl << std::endl;
        resultJSON["no_trachea_lung_mask_filename"] = a;

        try {
          writer->Update();
        } catch (itk::ExceptionObject &ex) {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }
      }

      // ok special sauce to separate the two lung regions from each other (might touch at the
      // front) we can start at the level that has three objects and look at each slice going down.
      // If we ignore the trachea we can calculate an area where the lungs touch (first slice after
      // we have one object instead of two). If we make the two objects larger (dilate) and subtract
      // them from each other we should get the touching voxel (in both) that we can remove - to get
      // two separate objects again, we can repeat this procedure for as long as we have our objects
      // (2) in the current slice. That should separate them from each other. Additionally we could
      // remove high intensity voxel...

      // lung mask without trachea: no_trachea

      // trachea are in: regionGrowingField

      // slice where we have three objects (trachea and two lungs): sliceSeedStart
      int slice = sliceSeedStart;
      ImageType::Pointer lastTwoLungRegions;
      MaskImageType::Pointer label1;
      MaskImageType::Pointer label2;
      for (; slice >= 0; slice--) { // correct orientation to go through the volume (top to bottom of body)
        // how many objects do we have in this slice if we ignore the trachea?
        size[2] = 1; // we extract along z direction
        ImageType::IndexType start = inputRegion.GetIndex();
        start[2] = slice;
        ImageType::RegionType desiredRegion;
        desiredRegion.SetSize(size);
        desiredRegion.SetIndex(start);
        extractFilter->SetExtractionRegion(desiredRegion);
        extractFilter->SetInput(no_trachea);

        // Now we can run connected components on the 2D image of the slice
        // we expect 2 distinct regions (the two lungs)
        ConnectedComponentImageFilterType::Pointer connected4 = ConnectedComponentImageFilterType::New();
        connected4->SetBackgroundValue(0);
        connected4->SetInput(extractFilter->GetOutput());
        // connected4->Update();

        labelType::Pointer labelTmp = labelType::New();
        labelTmp->SetInput(connected4->GetOutput());
        labelTmp->SetComputePerimeter(true);
        labelTmp->Update();

        LabelMapType *labelMapTmp = labelTmp->GetOutput();
        // labelMapTmp->DisconnectPipeline(); // we want to use extractFilter again, without running
        // the next part
        ImageType::Pointer twoLungRegions;
        // problem is that the number is not a good indicator of separation, instead it should be
        // the number after filtering out the small islands per slice, so we need to filter all
        // objects by area with a minimum area per slice
        int lungAreas = 0; // areas large enough to be a single lung
        for (unsigned int n = 0; n < labelMapTmp->GetNumberOfLabelObjects(); ++n) {
          ShapeLabelObjectType *labelObject = labelMapTmp->GetNthLabelObject(n);
          if (labelObject->GetNumberOfPixels() > 10000) { // per slice, todo: convert to volume, problem could be that we are too far
                                                          // up in the stack of slices with small volumes
            fprintf(stdout,
                    "could be a lung because its large (>0.2l) object: %d has %0.4f liters, %lu "
                    "voxel\n",
                    n, labelObject->GetPhysicalSize() / 1000000, labelObject->GetNumberOfPixels());
            lungAreas++;
          }
          // fprintf(stdout, "is this a lung? (>0.2l) object: %d has %0.4f liters, %lu voxel\n", n,
          // labelObject->GetPhysicalSize() / 1000000, labelObject->GetNumberOfPixels());
        }
        // fprintf(stdout, "FOUND %d lung area%s\t", lungAreas, lungAreas != 1 ? "s" : "");
        if (lungAreas == 1) {
          // todo: would be better to do this twice, once from below and once from above (limit the
          // effect of direction) would also be good to look for bright pixel in the separating
          // region, using image information would be better
          fprintf(stdout,
                  "SEPARATE THEM (using information from previous slice)!\n"); // using the labels from
                                                                               // the slice above
        }

        // if we have more than 2 objects, lets only use the largest 2, this ignores islands in
        // single slices if (labelMapTmp->GetNumberOfLabelObjects() >= 2) { keep only the large
        // connected component
        typedef itk::LabelShapeKeepNObjectsImageFilter<ImageType> LabelShapeKeepNObjectsImageMaskFilterType;
        LabelShapeKeepNObjectsImageMaskFilterType::Pointer labelShapeKeepNObjectsImageFilter3 = LabelShapeKeepNObjectsImageMaskFilterType::New();
        labelShapeKeepNObjectsImageFilter3->SetInput(connected4->GetOutput());
        labelShapeKeepNObjectsImageFilter3->SetBackgroundValue(0);
        labelShapeKeepNObjectsImageFilter3->SetNumberOfObjects(2);
        labelShapeKeepNObjectsImageFilter3->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
        labelShapeKeepNObjectsImageFilter3->Update();

        twoLungRegions = labelShapeKeepNObjectsImageFilter3->GetOutput();

        if (slice == sliceSeedStart) {
          lastTwoLungRegions = twoLungRegions;
        }
        fprintf(stderr, "slice %d, #objects: %lu...\n", slice, labelMapTmp->GetNumberOfLabelObjects());
        if (labelMapTmp->GetNumberOfLabelObjects() == 0)
          continue; // we could end up here at the bottom of the stack
        //}
        if (lungAreas == 1) { // we have 1 (or none) region, here we assume that we have a previous
                              // ok slice in label1 and label2
          fprintf(stderr, "slice %d, single object, separate them now...\n", slice);
          // the two objects are in lastTwoLungRegions and marked there with two different values
          // (from connected4)

          // do a region growing on label1 and label2
          DilateFilterType::Pointer dilate = DilateFilterType::New();
          StructuringElementType structElement;
          structElement.SetRadius(3); // 3x3 structuring element, larger radius will make the
                                      // separation larger (depends on voxel size)
          structElement.CreateStructuringElement();
          dilate->SetKernel(structElement);
          dilate->SetDilateValue(1);

          dilate->SetInput(label1);
          dilate->Update();
          MaskImageType::Pointer tmpLabel1 = dilate->GetOutput();
          tmpLabel1->DisconnectPipeline();
          dilate->SetInput(label2);
          dilate->Update();
          MaskImageType::Pointer tmpLabel2 = dilate->GetOutput();
          tmpLabel2->DisconnectPipeline();
          label1 = tmpLabel1;
          label2 = tmpLabel2;
          // now remove the overlap from the 'slice' mask (this should work if the images are
          // sufficiently similar)
          ImageType::RegionType re = label1->GetLargestPossibleRegion();
          itk::ImageRegionIterator<MaskImageType> imIterator1(label1, re);
          itk::ImageRegionIterator<MaskImageType> imIterator2(label2, re);
          itk::ImageRegionIterator<MaskImageType> dR(no_trachea,
                                                     desiredRegion); // should be the location of the current slice in the whole volume
          itk::ImageRegionIterator<ImageType> iR(inputImage, desiredRegion);
          int thresholdAir = -800;
          while (!imIterator1.IsAtEnd() && !imIterator2.IsAtEnd() && !dR.IsAtEnd()) {
            // fprintf(stdout, "Start iterating after dilation...\n");
            if (imIterator1.Get() > 0 && imIterator2.Get() > 0 && iR.Get() > thresholdAir) { // this voxel is an overlap between the two slices,
                                                                                             // remove from our mask
              // fprintf(stdout, "Found overlap voxel in no_trachea with label %d\n", dR.Get());
              dR.Set(0); // remove from no_trachea
              imIterator1.Set(0);
              imIterator2.Set(0);
            }
            ++imIterator1;
            ++imIterator2;
            ++dR;
            ++iR;
          }
          // now we have changed the slice (hopefully they are separated now)
          // we need to update the label1 and label2 regions again - because we have a new one for
          // the next slice
          extractFilter->SetExtractionRegion(desiredRegion);
          extractFilter->SetInput(no_trachea);

          //////////////////////////////////////////////////////////////////////////
          // we expect 2 distinct regions (the two lungs)
          // as a sanity check we should make sure that we really have two large enough objects now
          // - that the separation was successful
          ConnectedComponentImageFilterType::Pointer connected5 = ConnectedComponentImageFilterType::New();
          connected5->SetBackgroundValue(0);
          connected5->SetInput(extractFilter->GetOutput());
          // connected5->Update();
          labelType::Pointer labelTmp2 = labelType::New();
          labelTmp2->SetInput(connected5->GetOutput());
          labelTmp2->SetComputePerimeter(true);
          labelTmp2->Update();

          LabelMapType *labelMapTmp2 = labelTmp2->GetOutput();
          fprintf(stdout, "after separation found %lu objects in slice, use the largest two only\n", labelMapTmp2->GetNumberOfLabelObjects());

          // we expect to have two regions now - because the set(0) from above hopefully separated
          // them, if not we should increase the structuring element and try again
          typedef itk::LabelShapeKeepNObjectsImageFilter<ImageType> LabelShapeKeepNObjectsImageMaskFilterType;
          LabelShapeKeepNObjectsImageMaskFilterType::Pointer labelShapeKeepNObjectsImageFilter4 = LabelShapeKeepNObjectsImageMaskFilterType::New();
          labelShapeKeepNObjectsImageFilter4->SetInput(connected5->GetOutput());
          labelShapeKeepNObjectsImageFilter4->SetBackgroundValue(0);
          labelShapeKeepNObjectsImageFilter4->SetNumberOfObjects(2);
          labelShapeKeepNObjectsImageFilter4->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
          labelShapeKeepNObjectsImageFilter4->Update();

          twoLungRegions = labelShapeKeepNObjectsImageFilter4->GetOutput();
          connected4 = connected5; // we are looking at these in the next section
                                   // this should be sufficient, the next part below is copying the
                                   // values into label1 and label2 again for the next iteration
        }

        // remember the last slices regions (exactly 2 regions in here)
        lastTwoLungRegions = twoLungRegions;
        // store the two regions as separate slices region1 and region2
        // if (lungAreas == 1 /*  labelMapTmp->GetNumberOfLabelObjects() >= 2 */) {
        label1 = MaskImageType::New();
        label2 = MaskImageType::New();
        ImageType::RegionType re = twoLungRegions->GetLargestPossibleRegion();
        label1->SetRegions(re);
        label2->SetRegions(re);
        label1->Allocate();
        label2->Allocate();
        label1->FillBuffer(itk::NumericTraits<PixelType>::Zero);
        label2->FillBuffer(itk::NumericTraits<PixelType>::Zero);
        itk::ImageRegionIterator<ImageType> imIt(connected4->GetOutput(), re);
        itk::ImageRegionIterator<ImageType> imIt2(twoLungRegions, re);
        itk::ImageRegionIterator<MaskImageType> imIterator1(label1, re);
        itk::ImageRegionIterator<MaskImageType> imIterator2(label2, re);
        std::vector<int> labelValues;
        while (!imIt.IsAtEnd()) {
          if (imIt.Get() != 0) {
            // fprintf(stdout, "found non-zero voxel: %d\n", imIt.Get());
            if (!(std::find(labelValues.begin(), labelValues.end(), imIt.Get()) != labelValues.end()))
              labelValues.push_back(imIt.Get());
          }
          ++imIt;
        }
        imIt.GoToBegin();
        if (labelValues.size() > 1)
          fprintf(stdout, "found %lu label values: %d %d\n", labelValues.size(), labelValues[0], labelValues[1]);
        else
          fprintf(stdout, "found a single %lu label value: %d\n", labelValues.size(), labelValues[0]);
        // now copy values into label1 and label2
        if (labelValues.size() > 1) {
          while (!imIt.IsAtEnd() && !imIt2.IsAtEnd() && !imIterator1.IsAtEnd() && !imIterator2.IsAtEnd()) {
            // ok we copy into label1 if we find labelValue[0] in connected4 and >0 in imIt2
            if (imIt2.Get() > 0 && imIt.Get() == labelValues[0]) {
              imIterator1.Set(1);
            }
            if (imIt2.Get() > 0 && imIt.Get() == labelValues[1]) {
              imIterator2.Set(1);
            }
            ++imIt;
            ++imIt2;
            ++imIterator1;
            ++imIterator2;
          }
        } else {
          fprintf(stdout, "WE HAVE GIVEN UP HERE, single object even after separation!!! Increase "
                          "separating size and try again?\n");
        }
        // we need to make sure that label1 and label2 look ok now

        //}
      }

      // ok, if we did this right we should have now 2 separate regions of interest (left and right
      // lung) input is no_trachea
      ConnectedComponentImageFilterType::Pointer connected_final_mask = ConnectedComponentImageFilterType::New();
      connected_final_mask->SetBackgroundValue(0);
      connected_final_mask->SetInput(no_trachea);
      connected_final_mask->Update();

      // ToDo: we should check here which side is left/right to be able to name the lungs
      // appropriately
      ImageType::Pointer finalLabelField = connected_final_mask->GetOutput();

      ////////////////////////////////////////////////////////////////////////////////////////////////////
      // We did a 2D region growing before, we should do a 3D regions growing for each of the lungs next
      // This would fill in the blood vessels better. Best would be a rolling ball filter here...
      // We have to do this for each area in the lung seperately.
      // This will only close smaller vessels!
      if (1) {
        std::vector<int> labelIds;
        MaskImageType::RegionType reg = finalLabelField->GetLargestPossibleRegion();
        itk::ImageRegionIterator<ImageType> labelIter0(finalLabelField, reg);
        while (!labelIter0.IsAtEnd()) {
          int l = labelIter0.Get();
          if (std::find(labelIds.begin(), labelIds.end(), l) != labelIds.end()) {
            labelIds.push_back(l);
          }
          ++labelIter0;
        }
        if (verbose) {
          fprintf(stdout, "Without trachea there are %ld labels in the label field\n", labelIds.size());
        }
        for (int label = 0; label < labelIds.size(); label++) {
          if (labelIds[label] == 0) {
            // ignore background
            continue;
          }
          // otherwise extract this label to a new volume for region growing
          ImageType::Pointer tlabel = ImageType::New();
          tlabel->SetRegions(reg);
          tlabel->Allocate();
          tlabel->SetOrigin(finalLabelField->GetOrigin());
          tlabel->SetSpacing(finalLabelField->GetSpacing());
          tlabel->SetDirection(finalLabelField->GetDirection());
          itk::ImageRegionIterator<ImageType> labelIter01(finalLabelField, reg);
          itk::ImageRegionIterator<ImageType> labelIter02(tlabel, reg);
          while (!labelIter01.IsAtEnd() && !labelIter02.IsAtEnd()) {
            if (labelIter01.Get() == labelIds[label]) {
              labelIter02.Set(1);
            } else {
              labelIter02.Set(0);
            }

            ++labelIter02;
            ++labelIter01;
          }
          // now that we have the values inside the volume, lets do region growing (2x)
          // and shrinking 2x

          using ErodeFilterType3D = itk::BinaryErodeImageFilter<ImageType, ImageType, StructuringElementType>;
          using DilateFilterType3D = itk::BinaryDilateImageFilter<ImageType, ImageType, StructuringElementType>;

          DilateFilterType3D::Pointer dilate01 = DilateFilterType3D::New(); // grows inside the tissue
          DilateFilterType3D::Pointer dilate02 = DilateFilterType3D::New(); // grows inside the tissue
          DilateFilterType3D::Pointer dilate03 = DilateFilterType3D::New(); // grows inside the tissue
          DilateFilterType3D::Pointer dilate04 = DilateFilterType3D::New(); // grows inside the tissue
          ErodeFilterType3D::Pointer erode01 = ErodeFilterType3D::New();
          ErodeFilterType3D::Pointer erode02 = ErodeFilterType3D::New();
          ErodeFilterType3D::Pointer erode03 = ErodeFilterType3D::New();
          ErodeFilterType3D::Pointer erode04 = ErodeFilterType3D::New();
          StructuringElementType structuringElement2;
          structuringElement2.SetRadius(1); // 3x3 structuring element
          structuringElement2.CreateStructuringElement();
          dilate01->SetKernel(structuringElement2);
          dilate02->SetKernel(structuringElement2);
          dilate03->SetKernel(structuringElement2);
          dilate04->SetKernel(structuringElement2);
          erode01->SetKernel(structuringElement2);
          erode02->SetKernel(structuringElement2);
          erode03->SetKernel(structuringElement2);
          erode04->SetKernel(structuringElement2);
          dilate01->SetDilateValue(1);
          dilate02->SetDilateValue(1);
          dilate03->SetDilateValue(1);
          dilate04->SetDilateValue(1);
          erode01->SetErodeValue(1);
          erode02->SetErodeValue(1);
          erode03->SetErodeValue(1);
          erode04->SetErodeValue(1);
          dilate01->SetInput(tlabel);
          dilate02->SetInput(dilate01->GetOutput());
          dilate03->SetInput(dilate02->GetOutput());
          dilate04->SetInput(dilate03->GetOutput());
          erode01->SetInput(dilate04->GetOutput());
          erode02->SetInput(erode01->GetOutput());
          erode03->SetInput(erode02->GetOutput());
          erode04->SetInput(erode03->GetOutput());
          erode04->Update();

          // ok, now copy this back to the image, overwrite the previous value labelIds[label]
          itk::ImageRegionIterator<ImageType> labelIter03(finalLabelField, reg);
          itk::ImageRegionIterator<ImageType> labelIter04(erode04->GetOutput(), reg);
          while (!labelIter03.IsAtEnd() && !labelIter04.IsAtEnd()) {
            if (labelIter04.Get() == 1) {
              labelIter03.Set(labelIds[label]);
            }
            // we only add voxel but never remove
            ++labelIter03;
            ++labelIter04;
          }
        }
      }
      // copy back the trachea to get a full label field
      // the data is in regionGrowingField
      MaskImageType::RegionType tracheaRegion = inputImage->GetLargestPossibleRegion();
      itk::ImageRegionIterator<MaskImageType> tracheaIter(regionGrowingField, tracheaRegion);
      itk::ImageRegionIterator<ImageType> labelIter(finalLabelField, tracheaRegion);
      // we should also compute the size here (number of voxel per region)
      size_t count_label1 = 0;
      size_t count_label2 = 0;
      size_t count_trachea = 0;
      while (!tracheaIter.IsAtEnd() && !labelIter.IsAtEnd()) {
        if (tracheaIter.Get() > 0) {
          labelIter.Set(3);
          count_trachea++;
        }
        if (labelIter.Get() == 1)
          count_label1++;
        if (labelIter.Get() == 2)
          count_label2++;
        ++tracheaIter;
        ++labelIter;
      }
      resultJSON["lung_separated_count"] = json::array();
      json l1;
      l1["id"] = 1;
      l1["label"] = "Lung1";
      l1["count"] = count_label1;
      l1["volume"] = count_label1 * spacingx * spacingy * spacingz * 1e-6;
      resultJSON["lung_separated_count"].push_back(l1);
      json l2;
      l2["id"] = 2;
      l2["label"] = "Lung2";
      l2["count"] = count_label2;
      l2["volume"] = count_label2 * spacingx * spacingy * spacingz * 1e-6;
      resultJSON["lung_separated_count"].push_back(l2);
      json l3;
      l3["id"] = 3;
      l3["trachea"] = "Trachea";
      l3["count"] = count_trachea;
      l3["volume"] = count_trachea * spacingx * spacingy * spacingz * 1e-6;
      resultJSON["lung_separated_count"].push_back(l3);

      if (1) { // debug output
        typedef itk::ImageFileWriter<ImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();

        std::string a;
        size_t lastdot = labelfieldfilename.find_last_of(".");
        if (lastdot == std::string::npos)
          a = labelfieldfilename + "_separated.nii";
        else
          a = labelfieldfilename.substr(0, lastdot) + "_separated.nii";

        // make sure we can create that path
        path p(a);
        create_directories(p.parent_path());

        // std::string a(labelfieldfilename + "trachea.nii");
        writer->SetFileName(a);
        writer->SetInput(finalLabelField /* no_trachea */);

        std::cout << "Writing the final mask as " << std::endl;
        std::cout << a << std::endl << std::endl;
        resultJSON["label_field_separated"] = std::string(a);

        try {
          writer->Update();
        } catch (itk::ExceptionObject &ex) {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }
      }

      // we need a DICOM version of the label field as well, just a copy over the input
      if (1) { // save the label field - without image information
        gdcm::UIDGenerator suid;
        std::string newSeriesUID = suid.Generate();
        int seriesNumber = 0;
        // frameOfReference as been defined below
        // std::string frameOfReferenceUID = fuid.Generate();

        // create the output directory for the DICOM data
        itksys::SystemTools::MakeDirectory(output);
        outputSeries = output + "/" + seriesIdentifier + "_separated";
        itksys::SystemTools::MakeDirectory(outputSeries);
        // fprintf(stdout, "save data to %s\n", output.c_str());
        std::cout << "Writing output images to " << outputSeries << std::endl;
        resultJSON["output_series_separated"] = std::string(outputSeries);
        // now read in each input file in a loop, copy the result data over and write out as DICOM
        for (int i = 0; i < fileNames.size(); i++) {
          // std::cout << "use slice: " << fileNames[i] << " as template for output" << std::endl;

          // this is 2D work
          typedef signed short InputPixelType;
          const unsigned int Dimension = 2;
          typedef itk::Image<InputPixelType, Dimension> InputImageType;

          typedef itk::ImageFileReader<InputImageType> ReaderType;
          ReaderType::Pointer reader = ReaderType::New();
          reader->SetFileName(fileNames[i]);

          typedef itk::GDCMImageIO ImageIOType;
          ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
          reader->SetImageIO(gdcmImageIO);

          try {
            reader->Update();
          } catch (itk::ExceptionObject &e) {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
            return EXIT_FAILURE;
          }
          // ReaderType::DictionaryRawPointer inputDict =
          // (*(reader->GetMetaDataDictionaryArray()))[0];

          InputImageType::Pointer inputImage = reader->GetOutput();
          InputImageType::RegionType region;
          region = inputImage->GetBufferedRegion();
          InputImageType::SizeType size = region.GetSize();
          // std::cout << "size is: " << size[0] << " " << size[1] << std::endl;

          InputImageType::PixelContainer *container;
          container = inputImage->GetPixelContainer();
          container->SetContainerManageMemory(false);
          unsigned int bla = sizeof(InputImageType::PixelType);
          InputImageType::PixelType *buffer2 = container->GetBufferPointer();

          ImageType::Pointer nImage = finalLabelField;
          InputImageType::PixelContainer *container2;
          container2 = nImage->GetPixelContainer();
          InputImageType::PixelType *buffer3 = container2->GetBufferPointer();

          // Here we copy all values over, that is 0, 1, 2, 3 but also additional labels
          // that have been selected before (air in intestines for example).
          memcpy(buffer2, &(buffer3[i * size[0] * size[1]]), size[0] * size[1] * bla);
          // We can clean the data (remove all other label).
          for (int k = 0; k < size[0] * size[1]; k++) {
            if (buffer2[k] > 3) {
              buffer2[k] = 0; // set to background
            }
          }

          typedef itk::MetaDataDictionary DictionaryType;
          DictionaryType &dictionary = dicomIO->GetMetaDataDictionary();

          std::string studyUID;
          std::string sopClassUID;
          std::string seriesNumber;
          std::string acquisitionNumber;
          std::string instanceNumber;
          std::string frameOfReferenceUID;
          itk::ExposeMetaData<std::string>(dictionary, "0020|000d", studyUID);
          itk::ExposeMetaData<std::string>(dictionary, "0008|0016", sopClassUID);
          itk::ExposeMetaData<std::string>(dictionary, "0020|0011", seriesNumber);
          itk::ExposeMetaData<std::string>(dictionary, "0020|0012", acquisitionNumber);
          itk::ExposeMetaData<std::string>(dictionary, "0020|0013", instanceNumber);
          itk::ExposeMetaData<std::string>(dictionary, "0020|0052", frameOfReferenceUID); // read out

          int newSeriesNumber = 1000 + atoi(seriesNumber.c_str()) + 2;
          int newAcquisitionNumber = 1000 + atoi(acquisitionNumber.c_str()) + 2;
          int newInstanceNumber = atoi(instanceNumber.c_str());

          // without this we get different study instance uids for each image in the series
          gdcmImageIO->KeepOriginalUIDOn();

          std::string sopInstanceUID = suid.Generate();

          // std::string entryId( "0008|103e" );
          // std::string value( "Intensity Corrected" );
          // itk::EncapsulateMetaData<std::string>( dictionary, entryId, value );

          // DictionaryType *dict = new DictionaryType();

          // Copy the dictionary from the first slice
          // CopyDictionary (dictionary, *dict);

          // Set the UID's for the study, series, SOP  and frame of reference
          itk::MetaDataDictionary &dictionarySlice = reader->GetOutput()->GetMetaDataDictionary();

          // fprintf(stdout, "Try to set this studyUID %s\n", studyUID.c_str());
          // fprintf(stdout, "Try to set this newSeriesUID %s\n", newSeriesUID.c_str());
          // fprintf(stdout, "Try to set this frame of reference UID: %s\n", frameOfReferenceUID.c_str());
          // fprintf(stdout, "Try to set this SOP instance UID: %s\n", sopInstanceUID.c_str());

          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|000d", studyUID);
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|000e", newSeriesUID);
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|0011", std::to_string(newSeriesNumber));
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|0012", std::to_string(newAcquisitionNumber));
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|0013", std::to_string(newInstanceNumber));
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|0052", frameOfReferenceUID); // apply
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0008|0018", sopInstanceUID);
          // these keys don't exist - results in error
          // itk::EncapsulateMetaData<std::string>(dictionary,"0020|0052", "0"); // Intercept
          // itk::EncapsulateMetaData<std::string>(dictionary,"0020|0053", "1"); // Slope

          std::string oldSeriesDesc;
          itk::ExposeMetaData<std::string>(dictionary, "0008|103e", oldSeriesDesc);

          std::ostringstream value;
          std::string extension = " (lung labels)";
          value.str("");
          value << oldSeriesDesc;
          // This is a long string and there is a 64 character limit in the
          // standard
          unsigned lengthDesc = value.str().length();

          // std::string seriesDesc(value.str(), 0, lengthDesc > 64 ? 64 : lengthDesc);
          // fprintf(stdout, "Try to set this series description: %s\n", seriesDesc.c_str());
          // itk::EncapsulateMetaData<std::string>(dictionarySlice, "0008|103e", seriesDesc);

          std::string seriesDesc(value.str(), 0, lengthDesc + extension.length() > 64 ? 64 - extension.length() : lengthDesc + extension.length());
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0008|103e", seriesDesc + extension);

          // itk::EncapsulateMetaData<std::string>(dictionary, "0008|103e", seriesDesc);

          // set a lung window -600 ... 1600
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0028|1051", std::to_string(4));
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0028|1050", std::to_string(1.5));

          // copy the values for this slice over
          // CopyDictionary (*dict, dictionary);

          // write out the result as a DICOM again
          typedef itk::ImageFileWriter<InputImageType> Writer1Type;
          Writer1Type::Pointer writer1 = Writer1Type::New();

          writer1->SetInput(inputImage);
          std::ostringstream o;
          o << outputSeries << "/dicom" << i << ".dcm";
          writer1->SetFileName(o.str());
          writer1->SetImageIO(gdcmImageIO);
          writer1->Update();
        }
      }

      // write the original intensity image again only in the area where we have a label (>0)
      // this should replace the previous extracted lung image because that one does not separate
      // the lungs
      ImageType::Pointer lungDensity = ImageType::New();
      ImageType::RegionType lungDensityRegion = inputImage->GetLargestPossibleRegion();
      lungDensity->SetRegions(lungDensityRegion);
      lungDensity->Allocate();
      lungDensity->FillBuffer(-1024); // density for air
      lungDensity->SetOrigin(inputImage->GetOrigin());
      lungDensity->SetSpacing(inputImage->GetSpacing());
      lungDensity->SetDirection(inputImage->GetDirection());
      itk::ImageRegionIterator<ImageType> labelIterator(finalLabelField, lungDensityRegion);
      itk::ImageRegionIterator<ImageType> lungDensityIterator(inputImage, lungDensityRegion);
      itk::ImageRegionIterator<ImageType> lungDensityOutIterator(lungDensity, lungDensityRegion);
      while (!labelIterator.IsAtEnd() && !lungDensityIterator.IsAtEnd()) {
        // is this a copy of do we really write the color here?
        PixelType value = labelIterator.Value();
        if (value > 0) {
          lungDensityOutIterator.Set(lungDensityIterator.Get());
        }
        ++labelIterator;
        ++lungDensityIterator;
        ++lungDensityOutIterator;
      }

      if (saveNifty) {
        typedef itk::ImageFileWriter<ImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();

        std::string volFileName = output + "/" + seriesIdentifier + ".nii";
        path p(volFileName);
        create_directories(p.parent_path());

        // std::string a(labelfieldfilename + "trachea.nii");
        writer->SetFileName(volFileName);
        writer->SetInput(lungDensity /* no_trachea */);

        std::cout << "Writing the lung density image as " << std::endl;
        std::cout << volFileName << std::endl << std::endl;
        resultJSON["lung_intensity"] = volFileName;

        try {
          writer->Update();
        } catch (itk::ExceptionObject &ex) {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }
      }

      gdcm::UIDGenerator fuid;
      std::string frameOfReferenceUID = fuid.Generate();

      /////////////////////////////////////////////////////////
      // lets write a DICOM dataset for the labels as well (color overlay)
      if (1) {
        gdcm::UIDGenerator suid;
        std::string newSeriesUID = suid.Generate();

        // https://itk.org/Doxygen/html/Examples_2IO_2RGBImageSeriesReadWrite_8cxx-example.html
        using CPixelType = itk::RGBPixel<unsigned char>;
        using CImageType = itk::Image<CPixelType, Dimension>;
        using CWriterType = itk::ImageFileWriter<CImageType>;
        CImageType::Pointer fused = CImageType::New();
        CImageType::RegionType fusedRegion = inputImage->GetLargestPossibleRegion();
        fused->SetRegions(fusedRegion);
        fused->Allocate();
        fused->FillBuffer(itk::NumericTraits<CPixelType>::Zero);
        fused->SetOrigin(inputImage->GetOrigin());
        fused->SetSpacing(inputImage->GetSpacing());
        fused->SetDirection(inputImage->GetDirection());

        itk::ImageRegionIterator<ImageType> fusedLabelIterator(finalLabelField, fusedRegion);
        itk::ImageRegionIterator<ImageType> inputIterator(inputImage, fusedRegion);
        itk::ImageRegionIterator<ImageType> finalIterator(final, fusedRegion);
        itk::ImageRegionIterator<CImageType> fusedIterator(fused, fusedRegion);
        // now compute the fused image
        // max value of input is?
        using ImageCalculatorFilterType = itk::MinimumMaximumImageCalculator<ImageType>;
        ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
        imageCalculatorFilter->SetImage(inputImage);
        imageCalculatorFilter->Compute();
        int minGray = imageCalculatorFilter->GetMinimum();
        int maxGray = imageCalculatorFilter->GetMaximum();
        // what are good contrast and brightness values?
        // compute a histogram first
        using HistogramGeneratorType = itk::Statistics::ScalarImageToHistogramGenerator<ImageType>;
        HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
        histogramGenerator->SetInput(inputImage);
        int histogramSize = 1024;
        histogramGenerator->SetNumberOfBins(histogramSize);
        histogramGenerator->SetHistogramMin(minGray);
        histogramGenerator->SetHistogramMax(maxGray);
        histogramGenerator->Compute();
        using HistogramType = HistogramGeneratorType::HistogramType;
        const HistogramType *histogram = histogramGenerator->GetOutput();
        double lowerT = 0.01;
        double upperT = 0.999;
        double t1 = -1;
        double t2 = -1;
        double sum = 0;
        double total = 0;
        for (unsigned int bin = 0; bin < histogramSize; bin++) {
          total += histogram->GetFrequency(bin, 0);
        }
        for (unsigned int bin = 0; bin < histogramSize; bin++) {
          double f = histogram->GetFrequency(bin, 0) / total;
          // fprintf(stdout, "bin %d, value is %f\n", bin, f);
          sum += f;
          if (t1 == -1 && sum > lowerT) {
            t1 = minGray + (maxGray - minGray) * (bin / (float)histogramSize);
          }
          if (t2 == -1 && sum > upperT) {
            t2 = minGray + (maxGray - minGray) * (bin / (float)histogramSize);
            break;
          }
        }
        fprintf(stdout, "calculated best threshold low: %f, high: %f\n", t1, t2);

        float f = 0.6;
        while (!fusedLabelIterator.IsAtEnd() && !fusedIterator.IsAtEnd() && !finalIterator.IsAtEnd()) {
          // is this a copy of do we really write the color here?
          CPixelType value = fusedIterator.Value();
          float scaledP = (inputIterator.Get() - t1) / (t2 - t1);
          float red = scaledP;
          if (fusedLabelIterator.Get() == 3)
            red = f * scaledP + (1 - f) * (1);
          float blue = scaledP;
          if (fusedLabelIterator.Get() == 1)
            blue = f * scaledP + (1 - f) * (1);
          float green = scaledP;
          if (fusedLabelIterator.Get() == 2)
            green = f * scaledP + (1 - f) * (1);

          red = std::min<float>(1, std::max<float>(0, red));
          green = std::min<float>(1, std::max<float>(0, green));
          blue = std::min<float>(1, std::max<float>(0, blue));
          // fprintf(stdout, "red: %f, green: %f, blue: %f\n", red, green, blue);

          value.SetRed((int)(red * 255));
          value.SetGreen((int)(green * 255));
          value.SetBlue((int)(blue * 255));
          fusedIterator.Set(value);
          // update the final image as well (with all gray-value voxel from the labels)
          if (fusedLabelIterator.Get() > 0) {
            finalIterator.Set(inputIterator.Get());
          } else {
            finalIterator.Set(-1024);
          }

          ++fusedIterator;
          ++fusedLabelIterator;
          ++inputIterator;
          ++finalIterator;
        }
        // fprintf(stdout, "PhotometricInterpretation is: %d\n",
        // fused->GetPhotometricInterpretation());

        CWriterType::Pointer cwriter = CWriterType::New();
        // https://github.com/malaterre/GDCM/blob/master/Examples/Cxx/CreateARGBImage.cxx

        // create the output directory for the DICOM data
        itksys::SystemTools::MakeDirectory(output);
        outputSeries = output + "/" + seriesIdentifier + "_fused";
        itksys::SystemTools::MakeDirectory(outputSeries);
        // fprintf(stdout, "save data to %s\n", output.c_str());
        std::cout << "Writing fused label images to " << outputSeries << std::endl;
        resultJSON["output_label_series"] = std::string(outputSeries);
        resultJSON["output_label_images"] = json::array();
        std::string newSeriesUID2 = suid.Generate();

        // now read in each input file in a loop, copy the result data over and write out as DICOM
        for (int i = 0; i < fileNames.size(); i++) {
          // std::cout << "use slice: " << fileNames[i] << " as template for output" << std::endl;
          // fprintf(stdout, "start with %d of %lu\n", i, fileNames.size());
          // this is 2D work
          typedef signed short InputPixelType;
          const unsigned int Dimension = 2;
          typedef itk::Image<InputPixelType, Dimension> InputImageType;

          typedef itk::ImageFileReader<InputImageType> ReaderType;
          ReaderType::Pointer reader = ReaderType::New();
          reader->SetFileName(fileNames[i]);

          typedef itk::GDCMImageIO ImageIOType;
          ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
          reader->SetImageIO(gdcmImageIO);

          try {
            reader->Update();
          } catch (itk::ExceptionObject &e) {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
            return EXIT_FAILURE;
          }
          // ReaderType::DictionaryRawPointer inputDict = (*(reader->GetMetaDataDictionaryArray()))[0];
          itk::MetaDataDictionary &dictionarySlice = reader->GetOutput()->GetMetaDataDictionary();

          InputImageType::Pointer inputImage = reader->GetOutput();
          InputImageType::RegionType region;
          region = inputImage->GetBufferedRegion();
          InputImageType::SizeType size = region.GetSize();
          // std::cout << "size is: " << size[0] << " " << size[1] << std::endl;

          InputImageType::PixelContainer *container;
          container = inputImage->GetPixelContainer();
          container->SetContainerManageMemory(false);
          unsigned int bla = sizeof(InputImageType::PixelType);
          InputImageType::PixelType *buffer2 = container->GetBufferPointer();

          gdcm::PixelFormat pf = gdcm::PixelFormat::UINT8;
          pf.SetSamplesPerPixel(3);
          gdcm::SmartPointer<gdcm::Image> simage = new gdcm::Image;
          gdcm::Image &image = *simage;
          image.SetNumberOfDimensions(2);
          typedef itk::Image<PixelType, 2> ImageType2D;
          ImageType2D::RegionType inRegion = inputImage->GetLargestPossibleRegion();
          image.SetDimension(0, static_cast<unsigned int>(inRegion.GetSize()[0]));
          image.SetDimension(1, static_cast<unsigned int>(inRegion.GetSize()[1]));
          // image.SetDimension(2, m_Dimensions[2] );
          image.SetSpacing(0, inputImage->GetSpacing()[0]);
          image.SetSpacing(1, inputImage->GetSpacing()[1]);

          image.SetPixelFormat(pf);
          gdcm::PhotometricInterpretation pi = gdcm::PhotometricInterpretation::RGB;
          image.SetPhotometricInterpretation(pi);
          image.SetTransferSyntax(gdcm::TransferSyntax::ExplicitVRLittleEndian);
          // copy the DICOM tags over from inputImage to image
          gdcm::DataElement pixeldata(gdcm::Tag(0x7fe0, 0x0010));
          // get Pixel buffer from fused
          CImageType::Pointer fusedNImage = fused;
          CImageType::PixelContainer *container22 = fusedNImage->GetPixelContainer();
          CImageType::PixelType *buffer22 = container22->GetBufferPointer();

          // for this image create a new image instance UID
          gdcm::UIDGenerator sopuid;
          std::string sopInstanceUID = sopuid.Generate();

          // now copy all the DICOM tags over
          using DictionaryType = itk::MetaDataDictionary;
          const DictionaryType &dictionaryIn = dicomIO->GetMetaDataDictionary();
          using MetaDataStringType = itk::MetaDataObject<std::string>;
          auto itr = dictionaryIn.Begin();
          auto end = dictionaryIn.End();
          //            itk::MetaDataDictionary outMetaData;

          std::string imagePositionPatient; // we might be able to get them this way, but can we set them?
          itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0032", imagePositionPatient);
          // perhaps we have to use the parsed values to write them again further down?
          double origin3D[3];
          sscanf(imagePositionPatient.c_str(), "%lf\\%lf\\%lf", &(origin3D[0]), &(origin3D[1]), &(origin3D[2]));
          // fprintf(stdout, "image position patient field: %lf, %lf, %lf\n", origin3D[0], origin3D[1], origin3D[2]);

          std::string imageOrientation;
          itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0037", imageOrientation);
          double imageOrientationField[6];
          sscanf(imageOrientation.c_str(), "%lf\\%lf\\%lf\\%lf\\%lf\\%lf", &(imageOrientationField[0]), &(imageOrientationField[1]),
                 &(imageOrientationField[2]), &(imageOrientationField[3]), &(imageOrientationField[4]), &(imageOrientationField[5]));
          // fprintf(stdout, "image orientation field: %lf, %lf, %lf, %lf, %lf, %lf\n", imageOrientationField[0], imageOrientationField[1],
          //        imageOrientationField[2], imageOrientationField[3], imageOrientationField[4], imageOrientationField[5]);

          std::string sliceThicknessString;
          double sliceThickness = 0.0;
          itk::ExposeMetaData<std::string>(dictionarySlice, "0018|0050", sliceThicknessString);
          sscanf(sliceThicknessString.c_str(), "%lf", &sliceThickness);

          std::string imageInstanceString;
          int imageInstance = 0;
          itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0013", imageInstanceString);
          sscanf(imageInstanceString.c_str(), "%d", &imageInstance);
          // fprintf(stdout, "FOUND INSTANCE: %d\n", imageInstance); // start counting with 0 when we use this value to pick the slice

          std::string sliceLocationString;
          float sliceLocation = 0.0f;
          itk::ExposeMetaData<std::string>(dictionarySlice, "0020|1041", sliceLocationString);
          sscanf(sliceLocationString.c_str(), "%f", &sliceLocation);

          std::string imageAcquisitionString;
          int acquisitionNumber = 0;
          itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0012", imageAcquisitionString);
          sscanf(imageAcquisitionString.c_str(), "%d", &acquisitionNumber);

          // go into slice by applying offset
          // here the problem is that slices can be in a different order from the number in sliceNames
          // we should therefore use an i that corresponds to the image instance number - not the slice name
          // sort order
          uint32_t len = inRegion.GetSize()[0] * inRegion.GetSize()[1] * 3 * sizeof(unsigned char);
          pixeldata.SetByteValue((char *)(((unsigned char *)buffer22) + i * len), len);
          image.SetDataElement(pixeldata);

          // create an image (see
          // http://gdcm.sourceforge.net/html/GenFakeImage_8cxx-example.html#_a1)
          gdcm::SmartPointer<gdcm::File> file = new gdcm::File; // empty file
          gdcm::FileDerivation fd;
          const char ReferencedSOPClassUID[] = "1.2.840.10008.5.1.4.1.1.7"; // Secondary Capture
          fd.AddReference(ReferencedSOPClassUID, frameOfReferenceUID.c_str());
          fd.SetPurposeOfReferenceCodeSequenceCodeValue(
              121324); // segmentation  (see
                       // https://github.com/malaterre/GDCM/blob/master/Source/MediaStorageAndFileFormat/gdcmFileDerivation.cxx)
          // CID 7203 Image Derivation
          // { "DCM",113072,"Multiplanar reformatting" },
          fd.SetDerivationCodeSequenceCodeValue(113076);
          fd.SetFile(*file);
          // If all Code Value are ok the filter will execute properly
          if (!fd.Derive()) {
            std::cerr << "Sorry could not derive using input info" << std::endl;
            return 1;
          }
          gdcm::DataSet &ds = fd.GetFile().GetDataSet();
          gdcm::Anonymizer ano;
          ano.SetFile(fd.GetFile());
          std::string seriesDescription;
          int seriesNumber;
          while (itr != end) {
            itk::MetaDataObjectBase::Pointer entry = itr->second;
            MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>(entry.GetPointer());
            if (entryvalue) {
              //                itk::EncapsulateMetaData< std::string >( outMetaData, "0020|000D",
              //                newSeriesUID );

              // itk::EncapsulateMetaData<std::string>( dictionaryOut, entryId, value );
              // gdcm::DataElement d( gdcm::Tag(0x7fe0,0x0010) );
              // gdcm::Attribute<0x0020, 0x4000> imagecomments;
              // imagecomments.SetValue("Hi");
              // parse tagkey for tags, get comma separated

              std::string tagkey = itr->first;
              std::string labelId;
              bool found = itk::GDCMImageIO::GetLabelFromTag(tagkey, labelId);
              std::string tagvalue = entryvalue->GetMetaDataObjectValue();
              if (strcmp(tagkey.c_str(), "0008|103e") == 0) {
                seriesDescription = tagvalue;
              }
              if (strcmp(tagkey.c_str(), "0020|0011") == 0) {
                seriesNumber = atoi(tagvalue.c_str());
              }
              //              if (strcmp(tagkey.c_str(), "0020|1041") == 0) {
              //                // don't overwrite the slice position
              //                ++itr;
              //                continue;
              //              }
              // change window level from -400..600 to 150..180
              if (strcmp(tagkey.c_str(), "0028|1050") == 0) {
                tagvalue = std::string("150");
              }
              if (strcmp(tagkey.c_str(), "0028|1051") == 0) {
                tagvalue = std::string("180");
              }

              unsigned int f1;
              unsigned int f2;
              sscanf(tagkey.c_str(), "%x|%x", &f1, &f2);
              // fprintf(stdout, "detected tags %4x %4x\n", f1, f2);

              // gdcm::Attribute<f1, f2> at1; // all tag information from old dataset
              // at1.SetValue( tagvalue.c_str() );
              // ds.Replace( at1.GetAsDataElement() );
              ano.Replace(gdcm::Tag(f1, f2), tagvalue.c_str());

              // d.SetValue( tagvalue.c_str() );
              // image.SetDataElement( d );
              /* if( found ) {
                  std::cout << "(" << tagkey << ") " << labelId;
                  std::cout << " = " << tagvalue.c_str() << std::endl;
                } else {
                  std::cout << "(" << tagkey <<  ") " << "Unknown";
                  std::cout << " = " << tagvalue.c_str() << std::endl;
                }*/
            }
            ++itr;
          }
          //            itk::EncapsulateMetaData< std::string >( outMetaData, "0020|000D",
          //            newSeriesUID ); image.SetMetaDataDictionary( outMetaData );

          // For the purpose of this execise we will pretend that this image is referencing
          // two source image (we need to generate fake UID for that).
          // move this inside the loop above to copy values over
          gdcm::Attribute<0x0008, 0x2111> at1; // Derivative Description
          at1.SetValue("Segmented Lung - Fused Segmentation");
          ds.Replace(at1.GetAsDataElement());

          gdcm::Attribute<0x0008, 0x0060> at2; // Derivative Description
          at2.SetValue("CT");
          ds.Replace(at2.GetAsDataElement());

          gdcm::Attribute<0x0020, 0x000E> at3;
          at3.SetValue(newSeriesUID2);
          ds.Replace(at3.GetAsDataElement());

          gdcm::Attribute<0x0008, 0x103E> at4;
          std::string extension = " (fused segmentation)";
          std::ostringstream value;
          value.str("");
          value << seriesDescription;
          // This is a long string and there is a 64 character limit in the
          // standard
          unsigned lengthDesc = value.str().length();
          std::string seriesDesc(value.str(), 0, lengthDesc + extension.length() > 64 ? 64 - extension.length() : lengthDesc + extension.length());
          // itk::EncapsulateMetaData<std::string>(dictionary, "0008|103e", seriesDesc + extension);
          at4.SetValue(seriesDesc + extension);
          ds.Replace(at4.GetAsDataElement());

          // seriesInstance
          gdcm::Attribute<0x0020, 0x0011> at5;
          at5.SetValue(1000 + seriesNumber + 3);
          ds.Replace(at5.GetAsDataElement());

          // use a unique SOPInstanceUID
          gdcm::Attribute<0x0008, 0x0018> at6;
          at6.SetValue(sopInstanceUID);
          ds.Replace(at6.GetAsDataElement());

          // image position patient from input
          // These values are actually not getting written to the files (RGB has no origin, values are 0\0\0, but see set origin further down)
          gdcm::Attribute<0x0020, 0x0032> at7;
          at7.SetValue(origin3D[0], 0);
          at7.SetValue(origin3D[1], 1);
          at7.SetValue(origin3D[2], 2);
          ds.Replace(at7.GetAsDataElement());
          std::ostringstream value2;
          value2.str("");
          at7.Print(value2);
          // fprintf(stdout, "origin is now supposed to be: %lf\\%lf\\%lf %s\n", origin3D[0], origin3D[1], origin3D[2], value2.str().c_str());
          // For RGB we can set this to make sure they show up at the right location in Horos/OsiriX
          image.SetOrigin(0, origin3D[0]);
          image.SetOrigin(1, origin3D[1]);
          image.SetOrigin(2, origin3D[2]);

          gdcm::Attribute<0x0018, 0x0050> at8;
          at8.SetValue(sliceThickness);
          ds.Replace(at8.GetAsDataElement());

          gdcm::Attribute<0x0020, 0x0037> at9;
          at9.SetValue(imageOrientationField[0], 0);
          at9.SetValue(imageOrientationField[1], 1);
          at9.SetValue(imageOrientationField[2], 2);
          at9.SetValue(imageOrientationField[3], 3);
          at9.SetValue(imageOrientationField[4], 4);
          at9.SetValue(imageOrientationField[5], 5);
          ds.Replace(at9.GetAsDataElement());

          // gdcm::Attribute<0x0020, 0x0013> at10;
          // at10.SetValue(imageInstance);
          // ds.Replace(at10.GetAsDataElement());

          gdcm::Attribute<0x0020, 0x1041> at11;
          at11.SetValue(sliceLocation);
          ds.Replace(at11.GetAsDataElement());

          gdcm::Attribute<0x0020, 0x0012> at12;
          at12.SetValue(1000 + acquisitionNumber + 3);
          ds.Replace(at12.GetAsDataElement());

          gdcm::Attribute<0x0020, 0x0013> at13;
          at13.SetValue(imageInstance); // count starts at 1 and increments for all slices
          ds.Replace(at13.GetAsDataElement());

          gdcm::Attribute<0x0020, 0x0052> at14;
          at14.SetValue(frameOfReferenceUID.c_str());
          ds.Replace(at14.GetAsDataElement());

          gdcm::ImageWriter writer;
          writer.SetImage(image);
          writer.SetFile(fd.GetFile());
          char pp[2048];
          sprintf(pp, "%s/dicom%05d.dcm", outputSeries.c_str(), i);
          std::string fname(pp);
          // std::ostringstream o;
          // o << outputSeries << "/dicom" << i << ".dcm";
          writer.SetFileName(fname.c_str());
          if (!writer.Write()) {
            return 1;
          }
          // here the only thing I can do right now is to write these files temporarily
          // and to read them back in as itk image, attach header and write again

          // we should take the image data from fused instead and make the image RGB
          /*             ImageType::Pointer nImage = final;
            InputImageType::PixelContainer* container2 = nImage->GetPixelContainer();
            InputImageType::PixelType* buffer3 = container2->GetBufferPointer();

            memcpy(buffer2, &(buffer3[i*size[0]*size[1]]), size[0]*size[1]*bla);

            typedef itk::MetaDataDictionary DictionaryType;
            DictionaryType & dictionary = inputImage->GetMetaDataDictionary();

            std::string studyUID;
            std::string sopClassUID;
            itk::ExposeMetaData<std::string>(dictionary, "0020|000d", studyUID);
            itk::ExposeMetaData<std::string>(dictionary, "0008|0016", sopClassUID);
            gdcmImageIO->KeepOriginalUIDOn();

            gdcm::UIDGenerator sopuid;
            std::string sopInstanceUID = sopuid.Generate();

            //std::string entryId( "0008|103e" );
            //std::string value( "Intensity Corrected" );
            //itk::EncapsulateMetaData<std::string>( dictionary, entryId, value );

            //DictionaryType *dict = new DictionaryType();

            // Copy the dictionary from the first slice
            //CopyDictionary (dictionary, *dict);

            // Set the UID's for the study, series, SOP  and frame of reference

            itk::EncapsulateMetaData<std::string>(dictionary,"0020|000d", studyUID);
            itk::EncapsulateMetaData<std::string>(dictionary,"0020|000e", newSeriesUID);
            itk::EncapsulateMetaData<std::string>(dictionary,"0020|0052", frameOfReferenceUID);

            // these keys don't exist - results in error
            //itk::EncapsulateMetaData<std::string>(dictionary,"0020|0052", "0"); // Intercept
            //itk::EncapsulateMetaData<std::string>(dictionary,"0020|0053", "1"); // Slope

            std::string oldSeriesDesc;
            itk::ExposeMetaData<std::string>(dictionary, "0008|103e", oldSeriesDesc);

            std::ostringstream value;
            value.str("");
            value << oldSeriesDesc  << " (lung fused)";
            // This is a long string and there is a 64 character limit in the
            // standard
            unsigned lengthDesc = value.str().length();

            std::string seriesDesc( value.str(), 0,
                                lengthDesc > 64 ? 64
                                : lengthDesc);
            itk::EncapsulateMetaData<std::string>(dictionary,"0008|103e", seriesDesc);

            // set a lung window -600 ... 1600
            itk::EncapsulateMetaData<std::string>(dictionary,"0028|1051", std::to_string(1400));
            itk::EncapsulateMetaData<std::string>(dictionary,"0028|1050", std::to_string(-500));

            // copy the values for this slice over
            //CopyDictionary (*dict, dictionary);

            // write out the result as a DICOM again
            typedef itk::ImageFileWriter< InputImageType >  Writer1Type;
            Writer1Type::Pointer writer1 = Writer1Type::New();

            writer1->SetInput( inputImage );
            std::ostringstream o;
            o << outputSeries << "/dicom" << i << ".dcm";
            writer1->SetFileName( o.str() );
            writer1->SetImageIO( gdcmImageIO );
            writer1->Update();
            //std::cout << "done with writing the image...";
            */
          // more conservative with information (only first and last in output)
          if (i == 0 || i == fileNames.size() - 1)
            resultJSON["output_label_images"].push_back(fname);
        }
      } // loop over series

      if (0) { /*
        // we should make the volume isotropic first, before we do the hessian
        typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
        ResampleFilterType::Pointer resampler = ResampleFilterType::New();
        typedef itk::IdentityTransform<double, Dimension> TransformType;

        TransformType::Pointer transform = TransformType::New();
        transform->SetIdentity();
        resampler->SetTransform(transform);

        typedef itk::WindowedSincInterpolateImageFunction<ImageType, 3>
        WindowedSincInterpolatorType; WindowedSincInterpolatorType::Pointer windowedSincInterpolator
        = WindowedSincInterpolatorType::New(); resampler->SetInterpolator(windowedSincInterpolator);

        resampler->SetDefaultPixelValue(-1024); // Hounsfield Units for Air

        const ImageType::SpacingType &inputSpacing = inputImage->GetSpacing();

        double minSpacing = itk::NumericTraits<double>::max();
        for (int i = 0; i < 3; i++) {
          minSpacing = (minSpacing > inputSpacing[i] ? inputSpacing[i] : minSpacing);
        }

        ImageType::SpacingType outputSpacing;
        outputSpacing[0] = minSpacing;
        outputSpacing[1] = minSpacing;
        outputSpacing[2] = minSpacing;

        resampler->SetOutputSpacing(outputSpacing);

        resampler->SetOutputOrigin(inputImage->GetOrigin());
        resampler->SetOutputDirection(inputImage->GetDirection());

        ImageType::SizeType inputSize = inputImage->GetLargestPossibleRegion().GetSize();

        typedef ImageType::SizeType::SizeValueType SizeValueType;

        const double dx = inputSize[0] * inputSpacing[0] / outputSpacing[0];
        const double dy = inputSize[1] * inputSpacing[1] / outputSpacing[1];
        const double dz = inputSize[2] * inputSpacing[2] / outputSpacing[2];

        ImageType::SizeType finalSize;

        finalSize[0] = static_cast<SizeValueType>(dx);
        finalSize[1] = static_cast<SizeValueType>(dy);
        finalSize[2] = static_cast<SizeValueType>(dz);

        fprintf(stdout, "finalSize of output is: %lu %lu %lu\n", finalSize[0], finalSize[1],
        finalSize[2]); resampler->SetSize(finalSize); resampler->SetInput(inputImage);

        // Bright plate like structures have low values of $\lambda_1$ and $\lambda_2$ and large
        negative values of $\lambda_3$ using HessianFilterType =
        itk::HessianRecursiveGaussianImageFilter<ImageType>; HessianFilterType::Pointer
        hessianFilter = HessianFilterType::New(); hessianFilter->SetInput(resampler->GetOutput());
        hessianFilter->SetSigma(static_cast<double>(1.0f));
        // now we need to get the three eigenvalues of the hessian to create our filter
        typedef itk::FixedArray<double, 3> EigenValueArrayType;
        typedef itk::Image<EigenValueArrayType, 3> EigenValueImageType;

        typedef itk::SymmetricEigenAnalysisImageFilter<HessianFilterType::OutputImageType,
        EigenValueImageType> EigenAnalysisFilterType; typename EigenAnalysisFilterType::Pointer
        m_SymmetricEigenValueFilter;

        m_SymmetricEigenValueFilter = EigenAnalysisFilterType::New();
        m_SymmetricEigenValueFilter->SetDimension(3);
        m_SymmetricEigenValueFilter->OrderEigenValuesBy(EigenAnalysisFilterType::FunctorType::OrderByValue);
        // sorting by value sorts from lowest (most negative to highest eigenvalue 01 < 02 < 03)

        // m_SymmetricEigenValueFilter->OrderEigenValuesBy(
        EigenAnalysisFilterType::FunctorType::OrderByMagnitude );

        m_SymmetricEigenValueFilter->SetInput(hessianFilter->GetOutput());
        m_SymmetricEigenValueFilter->Update();

        typedef typename EigenAnalysisFilterType::OutputImageType EigenValueOutputImageType;
        const typename EigenValueOutputImageType::ConstPointer eigenImage =
        m_SymmetricEigenValueFilter->GetOutput();

        EigenValueArrayType eigenValue;
        itk::ImageRegionConstIterator<EigenValueOutputImageType> it;
        it = itk::ImageRegionConstIterator<EigenValueOutputImageType>(eigenImage,
        eigenImage->GetLargestPossibleRegion());

        // create a new output image for the eigenvalue filter result
        ImageType::Pointer vessels = ImageType::New();
        ImageType::RegionType vesselRegion = eigenImage->GetLargestPossibleRegion();
        vessels->SetRegions(vesselRegion);
        vessels->Allocate();
        vessels->SetOrigin(inputImage->GetOrigin());
        vessels->SetSpacing(inputImage->GetSpacing());
        itk::ImageRegionIterator<ImageType> vesselIterator(vessels, vesselRegion);

        ImageType::Pointer walls = ImageType::New();
        ImageType::RegionType wallRegion = eigenImage->GetLargestPossibleRegion();
        walls->SetRegions(wallRegion);
        walls->Allocate();
        walls->SetOrigin(inputImage->GetOrigin());
        walls->SetSpacing(inputImage->GetSpacing());
        itk::ImageRegionIterator<ImageType> wallIterator(walls, wallRegion);

        // calculate the maximum eigenvalue (absolute first)
        double maxS = 0.0;
        double a = 0.0;
        it.GoToBegin();
        while (!it.IsAtEnd()) {
          // Get the eigen value
          eigenValue = it.Get();
          a = sqrt((eigenValue[0] * eigenValue[0]) + (eigenValue[1] * eigenValue[1]) +
        (eigenValue[2] * eigenValue[2])); if (a > maxS) maxS = a;
          ++it;
        }
        // fprintf(stdout, "maxS : %f\n", maxS);
        double ce = maxS / 2.0;

        it.GoToBegin();
        vesselIterator.GoToBegin();
        wallIterator.GoToBegin();
        double wallMax = 0.0;
        double wallMin = 0.0;
        bool firstTime = true;
        while (!it.IsAtEnd() && !vesselIterator.IsAtEnd()) {
          // Get the eigen value
          eigenValue = it.Get();

          // normalizeValue <= 0 for bright line structures
          double e0 = eigenValue[0];
          double e1 = eigenValue[1];
          double e2 = eigenValue[2];
          //
          //  Bright tubular structures will have low $\lambda_1$ and large negative values of
        $\lambda_2$ and $\lambda_3$.
          //  Conversely dark tubular structures will have a low value of $\lambda_1$ and large
        positive values of $\lambda_2$ and $\lambda_3$.
          //  Bright plate like structures have low values of $\lambda_1$ and $\lambda_2$ and large
        negative values of $\lambda_3$
          //  Dark plate like structures have low values of $\lambda_1$ and $\lambda_2$ and large
        positive values of $\lambda_3$
          //  Bright spherical (blob) like structures have all three eigen values as large negative
        numbers
          //  Dark spherical (blob) like structures have all three eigen values as large positive
        numbers
          //
          // fprintf(stdout, "%f %f %f", e0, e1, e2);
          double mag = sqrt((eigenValue[0] * eigenValue[0]) + (eigenValue[1] * eigenValue[1]) +
        (eigenValue[2] * eigenValue[2]));
          // double Rb = eigenValue[0] / ( .5 *(eigenValue[1] + eigenValue[2]) ); // if sort is by
        magnitude double Rb = .5 * std::exp(-(1 - 0.7)) * std::fabs(eigenValue[0]); double beta = 1;
          double scale = 1;
          vesselIterator.Set(scale * std::exp(-.5 * Rb * Rb / (beta * beta)) * (1 - std::exp(-.5 *
        mag * mag / (ce * ce))));

          // wall like structures
          double av = .5 * (eigenValue[1] + eigenValue[2]);
          if (fabs(av) < 1e-6)
            av = 1e-6;
          double Rw = (-1.0 * eigenValue[0]) / av;
          if (firstTime) {
            firstTime = false;
            wallMin = wallMax = Rw;
          }
          if (Rw < wallMin)
            wallMin = Rw;
          if (Rw > wallMax)
            wallMax = Rw;

          // Rw = 0.5 * std::exp(-(1 - 0.7)) * std::fabs(eigenValue[0]);
          // wallIterator.Set( scale * exp(-.5 * Rw * Rw / (beta * beta)) * (1 - exp(-.5 * mag * mag
        / (ce * ce))) + scale); scale = 1; wallIterator.Set(scale * Rw);

          ++it;
          ++vesselIterator;
          ++wallIterator;
        }
        fprintf(stdout, "Wall min: %f and wall max is: %f\n", wallMin, wallMax);

        typedef itk::Image<ImageType, 3> OutputImageType;
        typedef itk::CastImageFilter<ImageType, ImageType> CastFilterType;
        CastFilterType::Pointer caster = CastFilterType::New();
        caster->SetInput(final);
        // caster->SetInput( maskImage );
        // caster->SetInput( binaryErode4->GetOutput() );
        // caster->SetInput( sliceFilter->GetOutput() );
        // caster->SetInput( labelShapeKeepNObjectsImageFilter->GetOutput() ); // body mask
        // caster->SetInput( mask1 );
        // caster->SetInput( mask2 );
        caster->Update();

        std::cout << "Processing done" << std::endl;

        if (saveNifty) {
          typedef itk::ImageFileWriter<ImageType> WriterType;
          WriterType::Pointer writer = WriterType::New();

          std::string volFileName = output + "/" + seriesIdentifier + ".nii";
          path p(volFileName);
          create_directories(p.parent_path());

          writer->SetFileName(volFileName);
          writer->SetInput(final);

          std::cout << "Writing the intensity inside the lung as " << std::endl << std::endl;
          std::cout << volFileName << std::endl << std::endl;
          resultJSON["lung_intensity"] = std::string(volFileName);

          try {
            writer->Update();
          } catch (itk::ExceptionObject &ex) {
            std::cout << ex << std::endl;
            return EXIT_FAILURE;
          }

          // and again
          writer = WriterType::New();

          p = niftyfilename;
          create_directories(p.parent_path());

          writer->SetFileName(niftyfilename);
          writer->SetInput(vessels);

          std::cout << "Writing the vessel image as " << std::endl << std::endl;
          std::cout << niftyfilename << std::endl << std::endl;
          resultJSON["vessel_image"] = std::string(niftyfilename);

          try {
            writer->Update();
          } catch (itk::ExceptionObject &ex) {
            std::cout << ex << std::endl;
            return EXIT_FAILURE;
          }

          // typedef itk::ImageFileWriter< ImageType > WriterType;
          WriterType::Pointer writer2 = WriterType::New();

          p = niftyfilename2;
          create_directories(p.parent_path());

          writer2->SetFileName(niftyfilename2);
          writer2->SetInput(walls);

          std::cout << "Writing the wall image " << std::endl << std::endl;
          resultJSON["wall_image"] = std::string(niftyfilename2);

          try {
            writer2->Update();
          } catch (itk::ExceptionObject &ex) {
            std::cout << ex << std::endl;
            return EXIT_FAILURE;
          }
        } */
      }

      if (0 /* saveNifty */) { // don't save this because its not with lungs separated
        typedef itk::ImageFileWriter<ImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();

        std::string volFileName = output + "/" + seriesIdentifier + ".nii";
        path p(volFileName);
        create_directories(p.parent_path());

        writer->SetFileName(volFileName);
        writer->SetInput(final); // not the right field

        std::cout << "Writing the intensity inside the lung as " << std::endl << std::endl;
        std::cout << volFileName << std::endl << std::endl;
        resultJSON["lung_intensity"] = std::string(volFileName);

        try {
          writer->Update();
        } catch (itk::ExceptionObject &ex) {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }
      }

      // now save as DICOM (should be done later as well?)
      gdcm::UIDGenerator suid;
      std::string newSeriesUID = suid.Generate();
      int seriesNumber = 0;
      // frameOfReference as been defined below
      // std::string frameOfReferenceUID = fuid.Generate();

      // create the output directory for the DICOM data
      itksys::SystemTools::MakeDirectory(output);
      outputSeries = output + "/" + seriesIdentifier;
      itksys::SystemTools::MakeDirectory(outputSeries);
      // fprintf(stdout, "save data to %s\n", output.c_str());
      std::cout << "Writing output images to " << outputSeries << std::endl;
      resultJSON["output_series"] = std::string(outputSeries);
      resultJSON["output_images"] = json::array();
      // now read in each input file in a loop, copy the result data over and write out as DICOM
      for (int i = 0; i < fileNames.size(); i++) {
        // std::cout << "use slice: " << fileNames[i] << " as template for output" << std::endl;

        // this is 2D work
        typedef signed short InputPixelType;
        const unsigned int Dimension = 2;
        typedef itk::Image<InputPixelType, Dimension> InputImageType;

        typedef itk::ImageFileReader<InputImageType> ReaderType;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(fileNames[i]);

        typedef itk::GDCMImageIO ImageIOType;
        ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
        reader->SetImageIO(gdcmImageIO);

        try {
          reader->Update();
        } catch (itk::ExceptionObject &e) {
          std::cerr << "exception in file reader " << std::endl;
          std::cerr << e.GetDescription() << std::endl;
          std::cerr << e.GetLocation() << std::endl;
          return EXIT_FAILURE;
        }
        // ReaderType::DictionaryRawPointer inputDict =
        // (*(reader->GetMetaDataDictionaryArray()))[0];

        InputImageType::Pointer inputImage = reader->GetOutput();
        InputImageType::RegionType region;
        region = inputImage->GetBufferedRegion();
        InputImageType::SizeType size = region.GetSize();
        // std::cout << "size is: " << size[0] << " " << size[1] << std::endl;

        InputImageType::PixelContainer *container;
        container = inputImage->GetPixelContainer();
        container->SetContainerManageMemory(false);
        unsigned int bla = sizeof(InputImageType::PixelType);
        InputImageType::PixelType *buffer2 = container->GetBufferPointer();

        ImageType::Pointer nImage = final;
        InputImageType::PixelContainer *container2;
        container2 = nImage->GetPixelContainer();
        InputImageType::PixelType *buffer3 = container2->GetBufferPointer();

        memcpy(buffer2, &(buffer3[i * size[0] * size[1]]), size[0] * size[1] * bla);

        typedef itk::MetaDataDictionary DictionaryType;
        DictionaryType &dictionary = reader->GetOutput()->GetMetaDataDictionary();

        std::string studyUID;
        std::string sopClassUID;
        std::string seriesNumber;
        std::string oldSeriesDesc;
        std::string acquisitionNumber;
        std::string instanceNumber;
        itk::ExposeMetaData<std::string>(dictionary, "0008|103e", oldSeriesDesc);
        itk::ExposeMetaData<std::string>(dictionary, "0020|000d", studyUID);
        itk::ExposeMetaData<std::string>(dictionary, "0008|0016", sopClassUID);
        itk::ExposeMetaData<std::string>(dictionary, "0020|0011", seriesNumber);
        itk::ExposeMetaData<std::string>(dictionary, "0020|0012", acquisitionNumber);
        itk::ExposeMetaData<std::string>(dictionary, "0020|0013", instanceNumber);

        int newSeriesNumber = 1000 + atoi(seriesNumber.c_str()) + 1;
        int newAcquisitionNumber = 1000 + atoi(acquisitionNumber.c_str()) + 1;
        int newInstanceNumber = 1000 + atoi(instanceNumber.c_str()) + 1;

        gdcmImageIO->KeepOriginalUIDOn();

        gdcm::UIDGenerator sopuid;
        std::string sopInstanceUID = sopuid.Generate();

        // std::string entryId( "0008|103e" );
        // std::string value( "Intensity Corrected" );
        // itk::EncapsulateMetaData<std::string>( dictionary, entryId, value );

        // DictionaryType *dict = new DictionaryType();

        // Copy the dictionary from the first slice
        // CopyDictionary (dictionary, *dict);

        // Set the UID's for the study, series, SOP  and frame of reference

        itk::EncapsulateMetaData<std::string>(dictionary, "0020|000d", studyUID);
        itk::EncapsulateMetaData<std::string>(dictionary, "0020|000e", newSeriesUID);
        itk::EncapsulateMetaData<std::string>(dictionary, "0020|0011", std::to_string(newSeriesNumber));
        itk::EncapsulateMetaData<std::string>(dictionary, "0020|0012", std::to_string(newAcquisitionNumber));
        itk::EncapsulateMetaData<std::string>(dictionary, "0020|0013", std::to_string(newInstanceNumber));
        itk::EncapsulateMetaData<std::string>(dictionary, "0020|0052", frameOfReferenceUID);
        itk::EncapsulateMetaData<std::string>(dictionary, "0008|0018",
                                              sopInstanceUID); // make the images unique so not to confuse them with the existing images

        // these keys don't exist - results in error
        // itk::EncapsulateMetaData<std::string>(dictionary,"0020|0052", "0"); // Intercept
        // itk::EncapsulateMetaData<std::string>(dictionary,"0020|0053", "1"); // Slope

        std::ostringstream value;
        std::string extension = " (lung segmentation)";
        value.str("");
        value << oldSeriesDesc;
        // This is a long string and there is a 64 character limit in the
        // standard
        unsigned lengthDesc = value.str().length();

        std::string seriesDesc(value.str(), 0, lengthDesc + extension.length() > 64 ? 64 - extension.length() : lengthDesc + extension.length());
        itk::EncapsulateMetaData<std::string>(dictionary, "0008|103e", seriesDesc + extension);

        // set a lung window -600 ... 1600
        itk::EncapsulateMetaData<std::string>(dictionary, "0028|1051", std::to_string(1400));
        itk::EncapsulateMetaData<std::string>(dictionary, "0028|1050", std::to_string(-500));

        // copy the values for this slice over
        // CopyDictionary (*dict, dictionary);

        // write out the result as a DICOM again
        typedef itk::ImageFileWriter<InputImageType> Writer1Type;
        Writer1Type::Pointer writer1 = Writer1Type::New();

        writer1->SetInput(inputImage);
        std::ostringstream o;
        o << outputSeries << "/dicom" << i << ".dcm";
        writer1->SetFileName(o.str());
        writer1->SetImageIO(gdcmImageIO);
        writer1->Update();
        if (i == 0 || i == fileNames.size() - 1) // be conservative with output info
          resultJSON["output_images"].push_back(o.str());
        // std::cout << "done with writing the image...";
      }

      // now we want to re-slice the left and right lung and export them again as a new volume (also
      // as DICOM later) ImageType::Pointer finalLabelField ImageType::Pointer inputImage (this
      // comes directly from the individual DICOM files)
      // (showLeft true for label 1)
      if (command.GetOptionWasSet("SaveReslice")) {
        ImageType::Pointer reslicedLung1 = computeReslice(inputImage, finalLabelField, 1, true, verbose, nullptr); // 0 is for background, 3 for trachea
        ImageType::Pointer reslicedLung2 = computeReslice(inputImage, finalLabelField, 2, false, verbose, reslicedLung1);
        // and save the resulting image as another DICOM series
      }

    } // loop over series
  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  std::string res = resultJSON.dump(4) + "\n";
  std::ostringstream o;
  std::string si(resultJSON["series_identifier"]);
  si.erase(std::remove(si.begin(), si.end(), '\"'), si.end());
  o << output << "/" << si << ".json";
  std::ofstream out(o.str());
  out << res;
  out.close();

  fprintf(stdout, "%s", res.c_str());

  return EXIT_SUCCESS;
}

void CopyDictionary(itk::MetaDataDictionary &fromDict, itk::MetaDataDictionary &toDict) {
  typedef itk::MetaDataDictionary DictionaryType;

  DictionaryType::ConstIterator itr = fromDict.Begin();
  DictionaryType::ConstIterator end = fromDict.End();
  typedef itk::MetaDataObject<std::string> MetaDataStringType;

  while (itr != end) {
    itk::MetaDataObjectBase::Pointer entry = itr->second;

    MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>(entry.GetPointer());
    if (entryvalue) {
      std::string tagkey = itr->first;
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      itk::EncapsulateMetaData<std::string>(toDict, tagkey, tagvalue);
    }
    ++itr;
  }
}
