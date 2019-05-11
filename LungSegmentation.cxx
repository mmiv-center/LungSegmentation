#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkMetaDataObject.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
//#include "itkExtractImageFilter.h"
//#include "itkPasteImageFilter.h"
#include "itkSliceBySliceImageFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkSymmetricEigenAnalysisImageFilter.h"
#include "itkImageAdaptor.h"
//#include <itkPixelAccessor.h>
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "itkResampleImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkExtractImageFilter.h"



#include "itkShrinkImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "gdcmUIDGenerator.h"
#include "itkMetaDataDictionary.h"

#include "json.hpp"
#include "metaCommand.h"
#include <map>

using json = nlohmann::json;

// forward declaration
void CopyDictionary (itk::MetaDataDictionary &fromDict, itk::MetaDataDictionary &toDict);

template<typename TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate  Self;
  typedef itk::Command            Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {}

public:

  virtual void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
    {
    Execute( (const itk::Object *) caller, event);
    }

  virtual void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
    {
    const TFilter * filter =
      dynamic_cast< const TFilter * >( object );

    if ( typeid( event ) != typeid( itk::IterationEvent ) )
                                        { return; }
    if ( filter->GetElapsedIterations() == 1 ) {
      std::cout << "Current level = " << filter->GetCurrentLevel() + 1
                << std::endl;
    }
    std::cout << "  Iteration " << filter->GetElapsedIterations()
              << " (of "
              << filter->GetMaximumNumberOfIterations()[ filter->GetCurrentLevel() ]
              << ").  ";
    std::cout << " Current convergence value = "
              << filter->GetCurrentConvergenceMeasurement()
              << " (threshold = " << filter->GetConvergenceThreshold()
              << ")" << std::endl;
    }

};


template<typename TValue>
TValue Convert( std::string optionString )
{
  TValue             value;
  std::istringstream iss( optionString );

  iss >> value;
  return value;
}

template<typename TValue>
std::vector<TValue> ConvertVector( std::string optionString )
{
  std::vector<TValue>    values;
  std::string::size_type crosspos = optionString.find( 'x', 0 );

  if ( crosspos == std::string::npos )
    {
    values.push_back( Convert<TValue>( optionString ) );
    }
  else
    {
    std::string        element = optionString.substr( 0, crosspos );
    TValue             value;
    std::istringstream iss( element );
    iss >> value;
    values.push_back( value );
    while ( crosspos != std::string::npos )
      {
      std::string::size_type crossposfrom = crosspos;
      crosspos = optionString.find( 'x', crossposfrom + 1 );
      if ( crosspos == std::string::npos )
        {
        element = optionString.substr( crossposfrom + 1, optionString.length() );
        }
      else
        {
        element = optionString.substr( crossposfrom + 1, crosspos );
        }
      std::istringstream iss2( element );
      iss2 >> value;
      values.push_back( value );
      }
    }
  return values;
}

json resultJSON;

int main( int argc, char* argv[] ) {
  
  itk::MultiThreader::SetGlobalMaximumNumberOfThreads( 4 );

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  command.SetDescription("LungSegmentation on CT DICOM images. Read DICOM image series and perform a lung segmentation. Exports a new DICOM series.");
  command.AddField("indir", "Directory with input DICOM image series.", MetaCommand::STRING, true);
  command.AddField("outdir", "Directory for output DICOM image series.", MetaCommand::STRING, true);  

  command.SetOption("Mask", "m", false, "Provide a mask image as an additional input. If not supplied the mask will be calculated.");
  command.AddOptionField("Mask", "mask", MetaCommand::STRING, false);

  command.SetOption("Write", "s", false, "The shrink factor will make the problem easier to handle (sub-sample data). The larger the value the faster.");
  command.AddOptionField("Write", "shrinkFactor", MetaCommand::INT, false, "3");

  command.SetOption("SeriesName", "n", false, "Select series by series name (if more than one series is present).");
  command.AddOptionField("SeriesName", "seriesname", MetaCommand::STRING, false);
  
  command.SetOption("SaveLabelfield", "b", false, "Save the label field as a nifty file in the current directory");
  command.AddOptionField("SaveLabelfield", "labelfieldfilename", MetaCommand::STRING,true);

  command.SetOption("SaveNifty", "u", false, "Save the corrected dataset as a nifty image to the current directory");
  command.AddOptionField("SaveNifty", "niftyfilename", MetaCommand::STRING, true);

  command.SetOption("Verbose", "v", false, "Print more verbose output");


  if ( !command.Parse(argc,argv) ) {
      return 1;
  }

  std::string input  = command.GetValueAsString("indir");
  std::string output = command.GetValueAsString("outdir");

  bool verbose = false;
  if (command.GetOptionWasSet("Verbose"))
     verbose = true;

  bool saveLabelField = false;  
  bool saveNifty = false;
  bool seriesIdentifierFlag = false;

  if ( command.GetOptionWasSet("SaveLabelfield") )
    saveLabelField = true;
  if ( command.GetOptionWasSet("SaveNifty") ) {
    saveNifty = true;
    fprintf(stdout, "will save nifty\n");
  }
  if ( command.GetOptionWasSet("SeriesName") )
    seriesIdentifierFlag = true;
  std::string labelfieldfilename = command.GetValueAsString("SaveLabelfield", "labelfieldfilename");
  std::string niftyfilename = command.GetValueAsString("SaveNifty", "niftyfilename");
  std::string niftyfilename2 = niftyfilename + "_walls.nii";
  size_t lastdot = niftyfilename.find_last_of(".");
  if (lastdot == std::string::npos) 
    niftyfilename2 = niftyfilename + "_walls.nii";
  else
    niftyfilename2 = niftyfilename.substr(0, lastdot) + "_walls.nii"; 

  std::string seriesName    = command.GetValueAsString("SeriesName", "seriesname");

  // store information in the result json file
  resultJSON["command_line"] = json::array();
  for (int i = 0; i < argc; i++)
  {
    resultJSON["command_line"].push_back(std::string(argv[i]));
  }
  typedef signed short PixelType;
  typedef float           FloatPixelType;
  const unsigned int      Dimension = 3;
  
  typedef itk::Image< PixelType, Dimension >         ImageType;
  typedef itk::Image< FloatPixelType, Dimension >    FloatImageType;
  
  typedef itk::ImageSeriesReader< ImageType >        ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  typedef itk::Image<unsigned char, Dimension> MaskImageType;
  typedef itk::Image<unsigned char, 2> MaskSliceImageType;

  using StructuringElementType = itk::BinaryBallStructuringElement<
                      PixelType,
                      Dimension  >;
  using ErodeFilterType = itk::BinaryErodeImageFilter<
                            MaskImageType,
                            MaskImageType,
                            StructuringElementType >;
  using DilateFilterType = itk::BinaryDilateImageFilter<
                            MaskImageType,
                            MaskImageType,
                            StructuringElementType >;

  typedef itk::ConnectedComponentImageFilter <MaskImageType, ImageType >
    ConnectedComponentImageFilterType;

  typedef itk::GDCMImageIO       ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  
  reader->SetImageIO( dicomIO );
  
  
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  
  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->AddSeriesRestriction("0008|0021" );
  nameGenerator->SetRecursive( true );
  nameGenerator->SetDirectory( input );
    
  try {
      std::cout << "Found DICOM Series: ";
      std::cout << std::endl;
      
      typedef std::vector< std::string >    SeriesIdContainer;
      
      const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
      
      SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
      SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
      while( seriesItr != seriesEnd ) {
          std::cout << "  " << seriesItr->c_str() << std::endl;
          ++seriesItr;
      }     
      
      std::string seriesIdentifier;
      
      SeriesIdContainer runThese;
      if( seriesIdentifierFlag ) { // If no optional series identifier
          // seriesIdentifier = seriesName;
          runThese.push_back( seriesName );
      } else {
          // Todo: here we select only the first series. We should run
          // N3 on all series.
        
          seriesItr = seriesUID.begin();
          seriesEnd = seriesUID.end();
          // seriesIdentifier = seriesUID.begin()->c_str();
          while (seriesItr != seriesEnd) {
            runThese.push_back( seriesItr->c_str() );
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
      while( seriesItr != seriesEnd) {
        seriesIdentifier = seriesItr->c_str();
        ++seriesItr;
      
        std::cout << "Processing series: " << std::endl;
        std::cout << "  " << seriesIdentifier << std::endl;

        std::string outputSeries = output + "/" + seriesIdentifier;
        if (itksys::SystemTools::FileIsDirectory( outputSeries.c_str() ) ) {
            fprintf(stdout, "Skip this series %s, output directory exists already...\n", outputSeries.c_str());
            exit(0);
        }

        typedef std::vector< std::string >   FileNamesContainer;
        FileNamesContainer fileNames;
      
        fileNames = nameGenerator->GetFileNames( seriesIdentifier );
      
        if (fileNames.size() < 10) {
          std::cout << "skip processing, not enough images in this series..." << std::endl;
          continue;
        }
        fprintf(stdout, "sufficient number of images [%lu] in this series, start processing...\n", fileNames.size());
        resultJSON["series_identifier"] = seriesIdentifier;
        //for (int i = 0; i < fileNames.size(); i++) {
        //  resultJSON["file_names"].push_back(fileNames[i]);
        //}

        // here we read in all the slices as a single volume for processing
        // if we want to write them back out we have to read them slice by
        // slice and get a copy of the meta data for each slice
        reader->SetFileNames( fileNames );
      
        try {
          reader->Update();
        } catch (itk::ExceptionObject &ex) {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }

        // read the data dictionary      
        ImageType::Pointer inputImage = reader->GetOutput();
        typedef itk::MetaDataDictionary   DictionaryType;
        DictionaryType & dictionary = inputImage->GetMetaDataDictionary();
        fprintf(stdout, "pixel spacing of input is: %f %f %f\n", inputImage->GetSpacing()[0], inputImage->GetSpacing()[1], inputImage->GetSpacing()[2]);
        // overwrite a value in the dicom dictionary
        //std::string entryId( "0010|0010" );
        //std::string value( "MYNAME" );
        //itk::EncapsulateMetaData<std::string>( dictionary, entryId, value );

        //      
        // replace this with the N4 algorithm instead of the gaussian filter
        //       
        MaskImageType::Pointer maskImage = ITK_NULLPTR;

        // if we have a mask on the command line      
        if (command.GetOptionWasSet("Mask")) {
          std::string maskName    = command.GetValueAsString("Mask", "mask");
          typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
          MaskReaderType::Pointer maskreader = MaskReaderType::New();
          maskreader->SetFileName( maskName );
          try {
            maskreader->Update();
            maskImage = maskreader->GetOutput();
            maskImage->DisconnectPipeline();
          } catch( ... ) {
            maskImage = ITK_NULLPTR;
          }
        }

        if ( !maskImage ) {
          typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType >  gaussianFilterType;
          gaussianFilterType::Pointer gaussianFilter = gaussianFilterType::New();
          gaussianFilter->SetInput( inputImage );
          gaussianFilter->SetVariance(2.0f);

           typedef itk::BinaryThresholdImageFilter<ImageType, MaskImageType> ThresholderType;
           ThresholderType::Pointer thresh = ThresholderType::New();
           thresh->SetInput( gaussianFilter->GetOutput() );
           thresh->SetLowerThreshold( -20000 );
           thresh->SetUpperThreshold( -200 );
           thresh->SetOutsideValue( 1 );
           thresh->SetInsideValue( 0 );

           thresh->Update();
           maskImage = thresh->GetOutput();
           maskImage->DisconnectPipeline();
           fprintf(stdout, "pixel spacing of mask image is: %f %f %f\n", maskImage->GetSpacing()[0], maskImage->GetSpacing()[1], maskImage->GetSpacing()[2]);
        }

        //
        // improve the mask for the body by hole filling and smoothing (grow + shrink)
        //
        DilateFilterType::Pointer binaryDilate  = DilateFilterType::New();
        DilateFilterType::Pointer binaryDilate2 = DilateFilterType::New();
        ErodeFilterType::Pointer  binaryErode   = ErodeFilterType::New();
        ErodeFilterType::Pointer  binaryErode2  = ErodeFilterType::New();
        ErodeFilterType::Pointer  binaryErode3  = ErodeFilterType::New();
        ErodeFilterType::Pointer  binaryErode4  = ErodeFilterType::New();

        StructuringElementType  structuringElement;
        structuringElement.SetRadius( 1 );  // 3x3 structuring element
        structuringElement.CreateStructuringElement();
        binaryDilate->SetKernel( structuringElement );
        binaryDilate2->SetKernel( structuringElement );
        binaryErode->SetKernel(  structuringElement );
        binaryErode2->SetKernel(  structuringElement );
        binaryErode3->SetKernel(  structuringElement );
        binaryErode4->SetKernel(  structuringElement );

        binaryDilate->SetInput( maskImage );
        binaryDilate2->SetInput( binaryDilate->GetOutput() );
        binaryErode->SetInput( binaryDilate2->GetOutput() );
        //binaryErode2->SetInput( binaryErode->GetOutput() );
        //binaryErode3->SetInput( binaryErode->GetOutput() );
        binaryErode4->SetInput( binaryErode->GetOutput() );

        binaryDilate->SetDilateValue( 1 );
        binaryDilate2->SetDilateValue( 1 );
        binaryErode->SetErodeValue( 1 );
        binaryErode2->SetErodeValue( 1 );
        binaryErode3->SetErodeValue( 1 );
        binaryErode4->SetErodeValue( 1 );
        binaryDilate->Update();
        binaryDilate2->Update();
        binaryErode->Update();
        //binaryErode2->Update();
        //binaryErode3->Update();
        binaryErode4->Update();

        // 
        // fill holes in the binary image slice by slice
        //
        typedef itk::SliceBySliceImageFilter< MaskImageType, MaskImageType > SliceFilterType;
        SliceFilterType::Pointer sliceFilter = SliceFilterType::New();
        sliceFilter->SetInput( binaryErode4->GetOutput() );

        typedef itk::BinaryFillholeImageFilter< SliceFilterType::InternalInputImageType > HoleFillFilterType;    
        HoleFillFilterType::Pointer holefillfilter = HoleFillFilterType::New();
        holefillfilter->SetForegroundValue( 1 );

        sliceFilter->SetFilter( holefillfilter );
        sliceFilter->Update();

        ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
        connected->SetInput( sliceFilter->GetOutput() );
        connected->Update();

        // keep only the large connected component
        typedef itk::LabelShapeKeepNObjectsImageFilter< ImageType > LabelShapeKeepNObjectsImageFilterType;
        LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter = LabelShapeKeepNObjectsImageFilterType::New();
        labelShapeKeepNObjectsImageFilter->SetInput( connected->GetOutput() );
        labelShapeKeepNObjectsImageFilter->SetBackgroundValue( 0 );
        labelShapeKeepNObjectsImageFilter->SetNumberOfObjects( 1 );
        labelShapeKeepNObjectsImageFilter->SetAttribute( LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
        labelShapeKeepNObjectsImageFilter->Update();
        // the above image labelShapeKeepNObjectsImageFilter->GetOutput() is now a good mask for the body of the participant
        // next we need to segment air inside that mask (do we need to make the mask tighter?, do we need to smooth the input? do we need to calculate this on a smaller scale?)
        // use binaryErode2->GetOutput(), labelShapeKeepNObjectsImageFilter->GetOutput(), and inputImage to get lung mask

        MaskImageType::Pointer mask1 = binaryErode4->GetOutput();
        ImageType::Pointer mask2     = labelShapeKeepNObjectsImageFilter->GetOutput();  // body mask
        MaskImageType::Pointer mask = MaskImageType::New();
        MaskImageType::RegionType maskRegion  = inputImage->GetLargestPossibleRegion();
        MaskImageType::RegionType mask1Region = mask1->GetLargestPossibleRegion();
        ImageType::RegionType mask2Region     = mask2->GetLargestPossibleRegion();
  
        mask->SetRegions(maskRegion);
        mask->Allocate();
        mask->SetOrigin(inputImage->GetOrigin());
        mask->SetSpacing(inputImage->GetSpacing());
        MaskImageType::SizeType regionSize = maskRegion.GetSize();
        itk::ImageRegionIterator<MaskImageType> maskIterator(mask,maskRegion);
        itk::ImageRegionIterator<MaskImageType> inputMask1Iterator(mask1,mask1Region);
        itk::ImageRegionIterator<ImageType>     inputMask2Iterator(mask2,mask2Region);
        // everything that is 1 in mask2 and 0 in mask1
        // problem is that mask2 can also be a value larger than 1! example is 0263 where the value is == 2
        while (!maskIterator.IsAtEnd() && !inputMask1Iterator.IsAtEnd() && !inputMask2Iterator.IsAtEnd())  {
          if (inputMask2Iterator.Get() > 0 && inputMask1Iterator.Get() == 0) {
          // if (maskIterator.GetIndex()[0] > static_cast<ImageType::IndexValueType>(regionSize[0]) / 2) {
              maskIterator.Set(1);
          } else {
              maskIterator.Set(0);
          }
          ++maskIterator;
          ++inputMask1Iterator;
          ++inputMask2Iterator;
        }
        fprintf(stdout, "pixel spacing of grow/shrunk image is: %f %f %f\n", mask->GetSpacing()[0], mask->GetSpacing()[1], mask->GetSpacing()[2]);

        // now in mask we have every air filled cavety inside the body, look for the largest one, assume its the lungs
        ConnectedComponentImageFilterType::Pointer connected2 = ConnectedComponentImageFilterType::New ();
        connected2->SetInput( mask );
        connected2->Update();
        //std::cout << "Number of connected components: " << connected->GetObjectCount() << std::endl;

        //std::cout << "Number of connected components: " << connected->GetObjectCount() << std::endl;
        // we sometimes get the two lungs separated here, we should look for the size of connected components
        // we do have a volume range we are looking for and getting two lungs with different volumes would
        // be possible (see 0183)

        using LabelType = unsigned short;        
        using ShapeLabelObjectType = itk::ShapeLabelObject< LabelType, 3 >;
        using LabelMapType = itk::LabelMap< ShapeLabelObjectType >;
        using labelType = itk::LabelImageToShapeLabelMapFilter< ImageType, LabelMapType>;
        labelType::Pointer label = labelType::New();
        label->SetInput( connected2->GetOutput() );
        label->SetComputePerimeter(true);
        label->Update();

        int useNObjects = 0;
        LabelMapType *labelMap = label->GetOutput();
        if (labelMap->GetNumberOfLabelObjects() == 0) {
          // error case
          fprintf(stderr, "Could not find any object in the data using the current set of thresholds, we can try to start again by lowering the threshold?\n");
        }

        for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n) {
          ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
          fprintf(stdout, "object %d has %0.4f liters, %lu voxel\n", n, labelObject->GetPhysicalSize() / 1000000, labelObject->GetNumberOfPixels());
          
          if (labelObject->GetNumberOfPixels() > 150000) // magick number using sample of 1 - should be selected based on volume instead of number of pixel
            useNObjects++;
        }
        if (useNObjects < 0)
          useNObjects = 1;
        fprintf(stdout, "found %d object%s with number of voxel large enough to be of interest...\n", useNObjects, useNObjects==1?"":"s");

        // keep only the large connected component
        typedef itk::LabelShapeKeepNObjectsImageFilter< ImageType > LabelShapeKeepNObjectsImageMaskFilterType;
        LabelShapeKeepNObjectsImageMaskFilterType::Pointer labelShapeKeepNObjectsImageFilter2 = LabelShapeKeepNObjectsImageMaskFilterType::New();
        labelShapeKeepNObjectsImageFilter2->SetInput( connected2->GetOutput() );
        labelShapeKeepNObjectsImageFilter2->SetBackgroundValue( 0 );
        labelShapeKeepNObjectsImageFilter2->SetNumberOfObjects( useNObjects );
        labelShapeKeepNObjectsImageFilter2->SetAttribute( LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
        labelShapeKeepNObjectsImageFilter2->Update();

        // the label filter will keep the id of the label, its not 1, what is it?
        /*typedef itk::MinimumMaximumImageCalculator <ImageType> ImageCalculatorFilterType;
        ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
        imageCalculatorFilter->SetImage( labelShapeKeepNObjectsImageFilter2->GetOutput() );
        imageCalculatorFilter->ComputeMaximum();
        int maxValue = imageCalculatorFilter->GetMaximum();*/
        int maxValue = 1;

        // instead of keeping the labels we should get a single label here with the value 1
        ImageType::Pointer lungs = labelShapeKeepNObjectsImageFilter2->GetOutput();
        ImageType::RegionType lungRegion  = lungs->GetLargestPossibleRegion();
        itk::ImageRegionIterator<ImageType> lungIterator(lungs, lungRegion);
        while (!lungIterator.IsAtEnd())  {
          if (lungIterator.Get() > 0) {
              lungIterator.Set(1);
          }
          ++lungIterator;
        }

        // 
        // and fill holes in the final segmentation of the lung
        //
        typedef itk::SliceBySliceImageFilter< ImageType, ImageType > SliceFilterImageType;
        SliceFilterImageType::Pointer sliceFilter2 = SliceFilterImageType::New();
        sliceFilter2->SetInput( labelShapeKeepNObjectsImageFilter2->GetOutput() );

        typedef itk::BinaryFillholeImageFilter< SliceFilterImageType::InternalInputImageType > HoleFillFilterType2;
        HoleFillFilterType2::Pointer holefillfilter2 = HoleFillFilterType2::New();
        holefillfilter2->SetForegroundValue( maxValue );

        sliceFilter2->SetFilter( holefillfilter2 );
        sliceFilter2->Update();

        // now apply the lung mask to the input image and export that instead
        ImageType::Pointer final = sliceFilter2->GetOutput();
        ImageType::RegionType imRegion   = inputImage->GetLargestPossibleRegion();
        ImageType::RegionType maskRegion2 = final->GetLargestPossibleRegion();
        itk::ImageRegionIterator<ImageType>     imIterator(inputImage,imRegion);
        itk::ImageRegionIterator<ImageType>     maskIterator2(final,maskRegion2);

        MaskImageType::Pointer labelField = MaskImageType::New();
        MaskImageType::RegionType labelFieldRegion  = inputImage->GetLargestPossibleRegion();
        labelField->SetRegions(labelFieldRegion);
        labelField->Allocate();
        labelField->SetOrigin(inputImage->GetOrigin());
        labelField->SetSpacing(inputImage->GetSpacing());
        itk::ImageRegionIterator<MaskImageType> labelFieldIterator(labelField, labelFieldRegion);

        // as the background voxel we can use the first voxel overall
        int firstVoxel;
        bool gotFirstVoxel = false;
        size_t numberOfPixel = 0;
        while (!imIterator.IsAtEnd() && !maskIterator2.IsAtEnd())
        {
          if (!gotFirstVoxel) {
            gotFirstVoxel = true;
            firstVoxel = imIterator.Get();
          }

          if (maskIterator2.Get() != 0) {
          // if (maskIterator.GetIndex()[0] > static_cast<ImageType::IndexValueType>(regionSize[0]) / 2) {
              maskIterator2.Set( imIterator.Get() );
              labelFieldIterator.Set(1);
              numberOfPixel++;
          }
          else
          {
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
        double originx  = inputImage->GetOrigin()[0];
        double originy  = inputImage->GetOrigin()[1];
        double originz  = inputImage->GetOrigin()[2];
        resultJSON["voxel_size"] = json::array();
        resultJSON["voxel_size"].push_back(spacingx);
        resultJSON["voxel_size"].push_back(spacingy);
        resultJSON["voxel_size"].push_back(spacingz);

        // in liters
        resultJSON["lung_volume"] = numberOfPixel * spacingx * spacingy * spacingz * 1e-6;

        // save the label field (one lung)
        if (saveLabelField) {
           typedef itk::ImageFileWriter< MaskImageType > WriterType;
           WriterType::Pointer writer = WriterType::New();
           writer->SetFileName( labelfieldfilename );
           writer->SetInput( labelField );
      
           std::cout  << "Writing the label field as " << std::endl;
           std::cout  << labelfieldfilename << std::endl << std::endl;
           resultJSON["label_field_filename"] = std::string(labelfieldfilename);
           fprintf(stdout, "label field voxel size is: %f %f %f\n", labelField->GetSpacing()[0], labelField->GetSpacing()[1], labelField->GetSpacing()[2]);
           //std::string res2 = resultJSON.dump(4) + "\n";
           //fprintf(stdout, "%s", res2.c_str());

           try
           {
             writer->Update();
           } catch (itk::ExceptionObject &ex) {
             std::cout << ex << std::endl;
             return EXIT_FAILURE;
           }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////
        // we should now split the lungs into the airways and the lungs using region growing
        // for this we need to find out if we have a label in the top image... we should have 3
        // regions of interest if we start from the top, we should set a seed point into the middle
        // one
        using ExtractFilterType = itk::ExtractImageFilter< MaskImageType, MaskImageType >;
        ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
        extractFilter->SetDirectionCollapseToSubmatrix();
        const MaskImageType * inImage = labelField;
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
        for (; sliceNumber >= 0; sliceNumber--)
        {
          size[2] = 1; // we extract along z direction
          ImageType::IndexType start = inputRegion.GetIndex();
          start[2] = sliceNumber;
          ImageType::RegionType desiredRegion;
          desiredRegion.SetSize(  size  );
          desiredRegion.SetIndex( start );
          extractFilter->SetExtractionRegion( desiredRegion );
          extractFilter->SetInput( inImage );

          // Now we can run connected components on the 2D image of the first slice
          // we expect 3 distinct regions
          ConnectedComponentImageFilterType::Pointer connected3 = ConnectedComponentImageFilterType::New ();
          connected3->SetInput( extractFilter->GetOutput() );
          connected3->Update();

          //using LabelType = unsigned short;        
          //using ShapeLabelObjectType = itk::ShapeLabelObject< LabelType, 3 >;
          //using LabelMapType = itk::LabelMap< ShapeLabelObjectType >;
          //using labelType = itk::LabelImageToShapeLabelMapFilter< ImageType, LabelMapType>;
          labelType::Pointer labelTmp = labelType::New();
          labelTmp->SetInput( connected3->GetOutput() );
          labelTmp->SetComputePerimeter(true);
          labelTmp->Update();

          //int useNObjects = 0;
          LabelMapType *labelMapTmp = labelTmp->GetOutput();
          if (labelMapTmp->GetNumberOfLabelObjects() == 0) {
            // error case
            fprintf(stderr, "Look at slice %d, Could not find any object in the data using the current set of thresholds, we can try to start again by lowering the threshold?\n", sliceNumber);
          }

          for (unsigned int n = 0; n < labelMapTmp->GetNumberOfLabelObjects(); ++n) {
            ShapeLabelObjectType *labelObject = labelMapTmp->GetNthLabelObject(n);
            fprintf(stdout, "top slice: object %d has %0.4f liters, %lu voxel\n", n, labelObject->GetPhysicalSize() / 1000000, labelObject->GetNumberOfPixels());
            
            //if (labelObject->GetNumberOfPixels() > 150000) // magick number using sample of 1 - should be selected based on volume instead of number of pixel
            //  useNObjects++;
          }
          // it might be better to be able to stop early, there are some series with only a single lung 
          // for those we have to stop at 2. What we could do is stop as well if we are too far down in the image stack
          if (labelMapTmp->GetNumberOfLabelObjects() == 3 || sliceNumber < 2*lastSlice/3) {
            sliceSeedStart = sliceNumber;
            fprintf(stdout, "found %lu objects in slice %d of %d\n", labelMapTmp->GetNumberOfLabelObjects(), sliceNumber + 1, lastSlice);
            // seed point would be in the object that is closest to the center of mass in the image
            // we could go to the last slice that has 3 objects as well... 
            // what if we have only one lung? We should detect the diameter of the airways instead.. something like at least two objects
          
            // which one is closest to the center (center of mass)
            for (unsigned int n = 0; n < labelMapTmp->GetNumberOfLabelObjects(); n++)
            {
              ShapeLabelObjectType *labelObject = labelMapTmp->GetNthLabelObject(n);
              double x = labelObject->GetCentroid()[0];
              double y = labelObject->GetCentroid()[1];
              //fprintf(stdout, "location of object in pixel %d is %f x %f\n", n, x, y);
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

              double COMx = originx + (spacingx) * ((size[0]-1) / 2.0);
              double COMy = originy + (spacingy) * ((size[1]-1) / 2.0);
              int seedx = roundf((x - originx) / spacingx);
              int seedy = roundf((y - originy) / spacingy);

              fprintf(stdout, "%lu %lu %lu, %f %f %f, %f %f %f, center %f %f\n", size[0], size[1], size[2], 
                      originx, originy, originz, 
                      spacingx, spacingy, spacingz, 
                      COMx, COMy);

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
            }
            break;
          }
        }
        fprintf(stdout, "min object is: %d, seed is at %f %f slice %d\n", 
                minRegion, minSeedx, minSeedy, sliceSeedStart);
        resultJSON["trachea_slice_location_pixel"] = json::array();
        resultJSON["trachea_slice_location_pixel"].push_back(minSeedx);
        resultJSON["trachea_slice_location_pixel"].push_back(minSeedy);
        resultJSON["trachea_slice_location_pixel"].push_back(sliceSeedStart);
        resultJSON["trachea_slice_number_objects"] = minNumObjects;
        resultJSON["trachea_slice_min_region"] = minRegion;
        resultJSON["trachea_slice_min_distance"] = minDist;
        
        // before we try to separate the two lungs we should shrink our label again,
        // that should help in cases where the trachea touch the lung wall
        // if we do this we also have to grow the region again before we do an assignment to tissue type


        ///////////////////////////////////////////////////////////////////////////////
        // ok now do a region growing in labelField starting at minSeedx maxSeedy for all voxel that have the value 1
        // we should have a second label field here
        // how do we find neighbors?
        MaskImageType::Pointer regionGrowingField = MaskImageType::New();
        MaskImageType::RegionType regionGrowingFieldRegion  = inputImage->GetLargestPossibleRegion();
        regionGrowingField->SetRegions(regionGrowingFieldRegion);
        regionGrowingField->Allocate();
        regionGrowingField->SetOrigin(inputImage->GetOrigin());
        regionGrowingField->SetSpacing(inputImage->GetSpacing());
        // itk::ImageRegionIterator<MaskImageType> regionGrowingFieldIterator(regionGrowingField, regionGrowingFieldRegion);
        // how to we find neighboring voxel from current location minSeedx and minSeedy and sliceSeedStart?
        using PointType = itk::Point< int, MaskImageType::ImageDimension >;
        std::vector<PointType> front;
        PointType p1;
        p1[0] = minSeedx;
        p1[1] = minSeedy;
        p1[2] = sliceSeedStart;
        front.push_back(p1);
        // we should have a 1 at this place
        ImageType::IndexType pixelIndex = {{p1[0],p1[1],p1[2]}};
        regionGrowingField->SetPixel(pixelIndex, 1);
        unsigned int count = 0;
        unsigned int frontSize = 0;
        size = regionGrowingFieldRegion.GetSize();
        fprintf(stdout, "size of output: %lu %lu %lu\n", size[0], size[1], size[2]);
        while (1)
        {
          if (front.size() == 0)
            break;

          count++;
          //fprintf(stdout, "frontSize is %d\n", front.size());
          // how do we check the first value from front in labelField? We want to know if the value is 0
          PointType p;
          p[0] = front[0][0];
          p[1] = front[0][1];
          p[2] = front[0][2];
          front.erase(front.begin());
          //pixelIndex = {{p[0], p[1], p[2]}};
          //fprintf(stdout, "read pixel value at location %d %d %d as : %d\n", p[0], p[1], p[2], inImage->GetPixel(pixelIndex));

          // get the 6 neighbors
          PointType pp1;
          pp1[0] = p[0]+1;
          pp1[1] = p[1];
          pp1[2] = p[2];
          // inside boundaries?
          if (pp1[0] < size[0]) {
            // are we still in the region of interest?
            pixelIndex = {{pp1[0], pp1[1], pp1[2]}};
            if (inImage->GetPixel(pixelIndex) == 1 && regionGrowingField->GetPixel(pixelIndex) == 0)
            {
              // add to front as a new pixel
              regionGrowingField->SetPixel(pixelIndex, 1);
              front.push_back(pp1);
            }
          }

          PointType pp2;
          pp2[0] = p[0]-1;
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
            }
          }

          PointType pp3;
          pp3[0] = p[0];
          pp3[1] = p[1]+1;
          pp3[2] = p[2];
          // inside boundaries?
          if (pp3[1] < size[1]) {
            // are we still in the region of interest?
            pixelIndex = {{pp3[0], pp3[1], pp3[2]}};
            if (inImage->GetPixel(pixelIndex) == 1 && regionGrowingField->GetPixel(pixelIndex) == 0) {
              // add to front as a new pixel
              regionGrowingField->SetPixel(pixelIndex, 1);
              front.push_back(pp3);
            }
          }

          PointType pp4;
          pp4[0] = p[0];
          pp4[1] = p[1]-1;
          pp4[2] = p[2];
          // inside boundaries?
          if (pp4[1] >= 0) {
            // are we still in the region of interest?
            pixelIndex = {{pp4[0], pp4[1], pp4[2]}};
            if (inImage->GetPixel(pixelIndex) == 1 && regionGrowingField->GetPixel(pixelIndex) == 0) {
              // add to front as a new pixel
              regionGrowingField->SetPixel(pixelIndex, 1);
              front.push_back(pp4);
            }
          }

          PointType pp5;
          pp5[0] = p[0];
          pp5[1] = p[1];
          pp5[2] = p[2]+1;
          // inside boundaries?
          if (pp5[2] < size[2]) {
            // are we still in the region of interest?
            pixelIndex = {{pp5[0], pp5[1], pp5[2]}};
            if (inImage->GetPixel(pixelIndex) == 1 && regionGrowingField->GetPixel(pixelIndex) == 0) {
              // add to front as a new pixel
              regionGrowingField->SetPixel(pixelIndex, 1);
              front.push_back(pp5);
            }
          }

          PointType pp6;
          pp6[0] = p[0];
          pp6[1] = p[1];
          pp6[2] = p[2]-1;
          // inside boundaries?
          if (pp6[2] >= 0) {
            // are we still in the region of interest?
            pixelIndex = {{pp6[0], pp6[1], pp6[2]}};
            if (inImage->GetPixel(pixelIndex) == 1 && regionGrowingField->GetPixel(pixelIndex) == 0) {
              // add to front as a new pixel
              regionGrowingField->SetPixel(pixelIndex, 1);
              front.push_back(pp6);
            }
          }
          // how long is the list in front? This burn in phase might not be required,
          // if we can add all pixel of the initial trachea into the region of interest
          if (count == 1000) {
            frontSize = front.size();
            fprintf(stdout, "frontSize after 1000 is %d\n", frontSize);
          }
          if (count > 1000) {
            fprintf(stdout, "iteration on frontSize: %lu\n", front.size());
            if (front.size() > frontSize * 2)
            {
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
          typedef itk::ImageFileWriter< MaskImageType > WriterType;
          WriterType::Pointer writer = WriterType::New();
          std::string a;
          size_t lastdot = labelfieldfilename.find_last_of(".");
          if (lastdot == std::string::npos) 
            a = labelfieldfilename + "_trachea.nii";
          else
            a = labelfieldfilename.substr(0, lastdot) + "_trachea.nii";

          // std::string a(labelfieldfilename + "trachea.nii");
          writer->SetFileName(a);
          writer->SetInput( regionGrowingField );
      
          std::cout  << "Writing the trachea field as " << std::endl;
          std::cout  << a << std::endl << std::endl;
          resultJSON["trachea_filename"] = a;

          try  {
            writer->Update();
          } catch (itk::ExceptionObject &ex) {
            std::cout << ex << std::endl;
            return EXIT_FAILURE;
          }
        }

        // after separating the airways out we should check again if we have two large regions of interest left
        // those would be left and right lung respectively. Each one should get its own label.
        // new Mask for the lung (not in trachea)
        MaskImageType::Pointer no_trachea = MaskImageType::New();
        MaskImageType::RegionType no_tracheaRegion  = inputImage->GetLargestPossibleRegion();
        no_trachea->SetRegions(no_tracheaRegion);
        no_trachea->Allocate();
        no_trachea->SetOrigin(inputImage->GetOrigin());
        no_trachea->SetSpacing(inputImage->GetSpacing());
        itk::ImageRegionIterator<MaskImageType> no_tracheaIterator(no_trachea, no_tracheaRegion);
        itk::ImageRegionIterator<MaskImageType> tracheaIterator(regionGrowingField,no_tracheaRegion);
        itk::ImageRegionIterator<MaskImageType> lungMaskIterator(labelField,no_tracheaRegion);
        while (!lungMaskIterator.IsAtEnd() && !tracheaIterator.IsAtEnd() && !no_tracheaIterator.IsAtEnd())  {
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
          typedef itk::ImageFileWriter< MaskImageType > WriterType;
          WriterType::Pointer writer = WriterType::New();
          std::string a;
          size_t lastdot = labelfieldfilename.find_last_of(".");
          if (lastdot == std::string::npos) 
            a = labelfieldfilename + "_lung_no_trachea.nii";
          else
            a = labelfieldfilename.substr(0, lastdot) + "_lung_no_trachea.nii";

          // std::string a(labelfieldfilename + "trachea.nii");
          writer->SetFileName(a);
          writer->SetInput( no_trachea );
      
          std::cout  << "Writing the no trachea lung mask as " << std::endl;
          std::cout  << a << std::endl << std::endl;
          resultJSON["no_trachea_lung_mask_filename"] = a;

          try  {
            writer->Update();
          } catch (itk::ExceptionObject &ex) {
            std::cout << ex << std::endl;
            return EXIT_FAILURE;
          }
        }



        if (0) {
          // we should make the volume isotropic first, before we do the hessian
          typedef itk::ResampleImageFilter< ImageType, ImageType >  ResampleFilterType;
          ResampleFilterType::Pointer resampler = ResampleFilterType::New();
          typedef itk::IdentityTransform< double, Dimension >  TransformType;

          TransformType::Pointer transform = TransformType::New();
          transform->SetIdentity();
          resampler->SetTransform( transform );

          typedef itk::WindowedSincInterpolateImageFunction< ImageType, 3 >  WindowedSincInterpolatorType;
          WindowedSincInterpolatorType::Pointer windowedSincInterpolator = WindowedSincInterpolatorType::New();
          resampler->SetInterpolator( windowedSincInterpolator );

          resampler->SetDefaultPixelValue( -1024 ); // Hounsfield Units for Air

          const ImageType::SpacingType & inputSpacing = inputImage->GetSpacing();

          double minSpacing = itk::NumericTraits< double >::max();
          for (int i = 0; i < 3; i++) {
            minSpacing = (minSpacing > inputSpacing[i] ? inputSpacing[i] : minSpacing);
          }
          
          ImageType::SpacingType outputSpacing;
          outputSpacing[0] = minSpacing;
          outputSpacing[1] = minSpacing;
          outputSpacing[2] = minSpacing;

          resampler->SetOutputSpacing( outputSpacing );

          resampler->SetOutputOrigin( inputImage->GetOrigin() );
          resampler->SetOutputDirection( inputImage->GetDirection() );

          ImageType::SizeType   inputSize = inputImage->GetLargestPossibleRegion().GetSize();
          
          typedef ImageType::SizeType::SizeValueType SizeValueType;

          const double dx = inputSize[0] * inputSpacing[0] / outputSpacing[0];
          const double dy = inputSize[1] * inputSpacing[1] / outputSpacing[1];
          const double dz = inputSize[2] * inputSpacing[2] / outputSpacing[2];

          ImageType::SizeType   finalSize;

          finalSize[0] = static_cast<SizeValueType>( dx );
          finalSize[1] = static_cast<SizeValueType>( dy );
          finalSize[2] = static_cast<SizeValueType>( dz );

          fprintf(stdout, "finalSize of output is: %lu %lu %lu\n", finalSize[0], finalSize[1], finalSize[2]);
          resampler->SetSize( finalSize );
          resampler->SetInput( inputImage );


          // Bright plate like structures have low values of $\lambda_1$ and $\lambda_2$ and large negative values of $\lambda_3$
          using HessianFilterType = itk::HessianRecursiveGaussianImageFilter<ImageType>;
          HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
          hessianFilter->SetInput( resampler->GetOutput() );
          hessianFilter->SetSigma( static_cast< double >( 1.0f ) );
          // now we need to get the three eigenvalues of the hessian to create our filter
          typedef  itk::FixedArray< double, 3 > EigenValueArrayType;
          typedef  itk::Image< EigenValueArrayType, 3 > EigenValueImageType;

          typedef  itk::SymmetricEigenAnalysisImageFilter< HessianFilterType::OutputImageType, EigenValueImageType > EigenAnalysisFilterType;
          typename EigenAnalysisFilterType::Pointer m_SymmetricEigenValueFilter;

          m_SymmetricEigenValueFilter = EigenAnalysisFilterType::New();
          m_SymmetricEigenValueFilter->SetDimension( 3 );
          m_SymmetricEigenValueFilter->OrderEigenValuesBy( EigenAnalysisFilterType::FunctorType::OrderByValue );
          // sorting by value sorts from lowest (most negative to highest eigenvalue 01 < 02 < 03)

          //m_SymmetricEigenValueFilter->OrderEigenValuesBy( EigenAnalysisFilterType::FunctorType::OrderByMagnitude );

          m_SymmetricEigenValueFilter->SetInput( hessianFilter->GetOutput() );
          m_SymmetricEigenValueFilter->Update();

          typedef typename EigenAnalysisFilterType::OutputImageType EigenValueOutputImageType;
          const typename EigenValueOutputImageType::ConstPointer eigenImage = m_SymmetricEigenValueFilter->GetOutput();

          EigenValueArrayType eigenValue;
          itk::ImageRegionConstIterator<EigenValueOutputImageType> it;
          it = itk::ImageRegionConstIterator<EigenValueOutputImageType>( eigenImage, eigenImage->GetLargestPossibleRegion());
          
          // create a new output image for the eigenvalue filter result
          ImageType::Pointer vessels = ImageType::New();
          ImageType::RegionType vesselRegion  = eigenImage->GetLargestPossibleRegion();
          vessels->SetRegions(vesselRegion);
          vessels->Allocate();
          vessels->SetOrigin(inputImage->GetOrigin());
          vessels->SetSpacing(inputImage->GetSpacing());
          itk::ImageRegionIterator<ImageType> vesselIterator(vessels, vesselRegion);

          ImageType::Pointer walls = ImageType::New();
          ImageType::RegionType wallRegion  = eigenImage->GetLargestPossibleRegion();
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
            a = sqrt( (eigenValue[0] * eigenValue[0]) + (eigenValue[1] * eigenValue[1]) + (eigenValue[2] * eigenValue[2]));
            if (a > maxS)
              maxS = a;
            ++it;
          }
          fprintf(stdout, "maxS : %f\n", maxS);
          double ce = maxS/2.0;

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
            /*
            Bright tubular structures will have low $\lambda_1$ and large negative values of $\lambda_2$ and $\lambda_3$.
            Conversely dark tubular structures will have a low value of $\lambda_1$ and large positive values of $\lambda_2$ and $\lambda_3$.
            Bright plate like structures have low values of $\lambda_1$ and $\lambda_2$ and large negative values of $\lambda_3$
            Dark plate like structures have low values of $\lambda_1$ and $\lambda_2$ and large positive values of $\lambda_3$
            Bright spherical (blob) like structures have all three eigen values as large negative numbers
            Dark spherical (blob) like structures have all three eigen values as large positive numbers
            */
            //fprintf(stdout, "%f %f %f", e0, e1, e2);
            double mag = sqrt( (eigenValue[0] * eigenValue[0]) + (eigenValue[1] * eigenValue[1]) + (eigenValue[2] * eigenValue[2]) );
            // double Rb = eigenValue[0] / ( .5 *(eigenValue[1] + eigenValue[2]) ); // if sort is by magnitude
            double Rb = .5 * std::exp(-(1-0.7)) * std::fabs(eigenValue[0]);
            double beta = 1;
            double scale = 1;
            vesselIterator.Set( scale * std::exp(-.5 * Rb * Rb / (beta * beta)) * (1 - std::exp(-.5 * mag * mag / (ce * ce)))  );

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

            //Rw = 0.5 * std::exp(-(1 - 0.7)) * std::fabs(eigenValue[0]);
            //wallIterator.Set( scale * exp(-.5 * Rw * Rw / (beta * beta)) * (1 - exp(-.5 * mag * mag / (ce * ce))) + scale);
            scale = 1;
            wallIterator.Set( scale * Rw );

            ++it;
            ++vesselIterator;
            ++wallIterator;
          }
          fprintf(stdout, "Wall min: %f and wall max is: %f\n", wallMin, wallMax);

          typedef itk::Image< ImageType, 3 > OutputImageType;
          typedef itk::CastImageFilter< 
                          ImageType,
                          ImageType > CastFilterType;
          CastFilterType::Pointer  caster =  CastFilterType::New();
          caster->SetInput( final );
          //caster->SetInput( maskImage );
          //caster->SetInput( binaryErode4->GetOutput() );
          //caster->SetInput( sliceFilter->GetOutput() );
          // caster->SetInput( labelShapeKeepNObjectsImageFilter->GetOutput() ); // body mask
          //caster->SetInput( mask1 );
          //caster->SetInput( mask2 );
          caster->Update();

          std::cout  << "Processing done" << std::endl;

          if (saveNifty) {
            typedef itk::ImageFileWriter< ImageType > WriterType;
            WriterType::Pointer writer = WriterType::New();      
            writer->SetFileName( niftyfilename );
            writer->SetInput( vessels );
        
            std::cout  << "Writing the vessel image as " << std::endl << std::endl;
            std::cout  << argv[2] << std::endl << std::endl;
        
            try  {
              writer->Update();
            } catch (itk::ExceptionObject &ex) {
              std::cout << ex << std::endl;
              return EXIT_FAILURE;
            }

            //typedef itk::ImageFileWriter< ImageType > WriterType;
            WriterType::Pointer writer2 = WriterType::New();      
            writer2->SetFileName( niftyfilename2 );
            writer2->SetInput( walls );
        
            std::cout  << "Writing the wall image " << std::endl << std::endl;
        
            try  {
              writer2->Update();
            } catch (itk::ExceptionObject &ex) {
              std::cout << ex << std::endl;
              return EXIT_FAILURE;
            }

          }
        }      
        // now save as DICOM
        gdcm::UIDGenerator suid;
        std::string newSeriesUID = suid.Generate();
        gdcm::UIDGenerator fuid;
        std::string frameOfReferenceUID = fuid.Generate();
 
        // create the output directory for the DICOM data
        itksys::SystemTools::MakeDirectory( output );
        outputSeries = output + "/" + seriesIdentifier;
        itksys::SystemTools::MakeDirectory( outputSeries );
        //fprintf(stdout, "save data to %s\n", output.c_str());
        std::cout  << "Writing output images to " << outputSeries << std::endl;
        resultJSON["output_series"] = std::string(outputSeries);
        resultJSON["output_images"] = json::array();
        // now read in each input file in a loop, copy the result data over and write out as DICOM
        for (int i = 0; i < fileNames.size(); i++) {
           //std::cout << "use slice: " << fileNames[i] << " as template for output" << std::endl;
        
           // this is 2D work
           typedef signed short InputPixelType;
           const unsigned int   Dimension = 2;
           typedef itk::Image< InputPixelType, Dimension > InputImageType;

           typedef itk::ImageFileReader< InputImageType > ReaderType;
           ReaderType::Pointer reader = ReaderType::New();
           reader->SetFileName( fileNames[i] );

           typedef itk::GDCMImageIO           ImageIOType;
           ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
           reader->SetImageIO( gdcmImageIO );
       
           try {
              reader->Update();
           } catch (itk::ExceptionObject & e) {
              std::cerr << "exception in file reader " << std::endl;
              std::cerr << e.GetDescription() << std::endl;
              std::cerr << e.GetLocation() << std::endl;
              return EXIT_FAILURE;
           }
           // ReaderType::DictionaryRawPointer inputDict = (*(reader->GetMetaDataDictionaryArray()))[0];        
        
           InputImageType::Pointer inputImage = reader->GetOutput();
           InputImageType::RegionType region;
           region = inputImage->GetBufferedRegion();
           InputImageType::SizeType size  = region.GetSize();
           // std::cout << "size is: " << size[0] << " " << size[1] << std::endl;

           InputImageType::PixelContainer* container;
           container = inputImage->GetPixelContainer();
           container->SetContainerManageMemory( false );
           unsigned int bla = sizeof(InputImageType::PixelType);
           InputImageType::PixelType* buffer2 = container->GetBufferPointer();
         
           ImageType::Pointer nImage = final;
           InputImageType::PixelContainer* container2;
           container2 = nImage->GetPixelContainer();
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
           value << oldSeriesDesc  << " (lung segmentation)";
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
           resultJSON["output_images"].push_back(o.str());
           //std::cout << "done with writing the image...";
        }
      } // loop over series 
      
  } catch (itk::ExceptionObject &ex) {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
  }

  std::string res = resultJSON.dump(4) + "\n";
  std::ostringstream o;
  o << output << "/" << resultJSON["series_identifier"] << ".json";
  std::ofstream out(o.str());
  out << res;
  out.close();

  fprintf(stdout, "%s", res.c_str());

  return EXIT_SUCCESS;
}


 
void CopyDictionary (itk::MetaDataDictionary &fromDict, itk::MetaDataDictionary &toDict)
{
  typedef itk::MetaDataDictionary DictionaryType;
 
  DictionaryType::ConstIterator itr = fromDict.Begin();
  DictionaryType::ConstIterator end = fromDict.End();
  typedef itk::MetaDataObject< std::string > MetaDataStringType;
 
  while( itr != end )
    {
    itk::MetaDataObjectBase::Pointer  entry = itr->second;
 
    MetaDataStringType::Pointer entryvalue =
      dynamic_cast<MetaDataStringType *>( entry.GetPointer() ) ;
    if( entryvalue )
      {
      std::string tagkey   = itr->first;
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      itk::EncapsulateMetaData<std::string>(toDict, tagkey, tagvalue);
      }
    ++itr;
    }
}