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



#include "itkShrinkImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "gdcmUIDGenerator.h"
#include "itkMetaDataDictionary.h"

#include "metaCommand.h"
#include <map>

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


int main( int argc, char* argv[] ) {
  
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
  
  command.SetOption("SaveBiasfield", "b", false, "Save the biasfield as a nifty file in the current directory");
  command.AddOptionField("SaveBiasfield", "biasfieldfilename", MetaCommand::STRING,true);

  command.SetOption("SaveNifty", "n", false, "Save the corrected dataset as a nifty image to the current directory");
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

  bool saveBiasField = false;  
  bool saveNifty = false;
  bool seriesIdentifierFlag = false;

  if ( command.GetOptionWasSet("SaveBiasfield") )
    saveBiasField = true;
  if ( command.GetOptionWasSet("SaveNifty") )
    saveNifty = true;
  if ( command.GetOptionWasSet("SeriesName") )
    seriesIdentifierFlag = true;
  std::string biasfieldfilename = command.GetValueAsString("SaveBiasfield", "biasfieldfilename");
  std::string niftyfilename = command.GetValueAsString("SaveNifty", "niftyfilename");
  std::string seriesName    = command.GetValueAsString("SeriesName", "seriesname");
  
  
  typedef signed short    PixelType;
  const unsigned int      Dimension = 3;
  
  typedef itk::Image< PixelType, Dimension >         ImageType;
  
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
          fprintf(stdout, "object %d has %0.4f liters, %d voxel\n", n, labelObject->GetPhysicalSize() / 1000000, labelObject->GetNumberOfPixels());
          
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

        // as the background voxel we can use the first voxel overall
        int firstVoxel;
        bool gotFirstVoxel = false;
        while (!imIterator.IsAtEnd() && !maskIterator2.IsAtEnd())  {
          if (!gotFirstVoxel) {
            gotFirstVoxel = true;
            firstVoxel = imIterator.Get();
          }

          if (maskIterator2.Get() != 0) {
          // if (maskIterator.GetIndex()[0] > static_cast<ImageType::IndexValueType>(regionSize[0]) / 2) {
              maskIterator2.Set( imIterator.Get() );
          } else {
              maskIterator2.Set( firstVoxel ); // this depends on the pixel representation its ok
          }
          ++maskIterator2;
          ++imIterator;
        }

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
           writer->SetInput( caster->GetOutput() );
     
           //////////////////////////////////////////////  
 
           /* typedef itk::SmoothingRecursiveGaussianImageFilter< ImageType, ImageType >  FilterType;
      
           FilterType::Pointer filter = FilterType::New();
      
           filter->SetInput( reader->GetOutput() );
      
           const double sigma = atof( argv[3] );
           filter->SetSigma( sigma );

      
           typedef itk::ImageFileWriter< ImageType > WriterType;
           WriterType::Pointer writer = WriterType::New();
      
           writer->SetFileName( argv[2] );
           writer->SetInput( filter->GetOutput() );
           //writer->SetImageIO( dicomIO ); */
      
           std::cout  << "Writing the image as " << std::endl << std::endl;
           std::cout  << argv[2] << std::endl << std::endl;      
      
           try  {
             writer->Update();
           } catch (itk::ExceptionObject &ex) {
             std::cout << ex << std::endl;
             return EXIT_FAILURE;
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
         
           ImageType::Pointer nImage = caster->GetOutput();
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
          //std::cout << "done with writing the image...";
        
        }
      } // loop over series 
      
  } catch (itk::ExceptionObject &ex) {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
  }
    
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