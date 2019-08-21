# Basic Lung Segmentation using itk

This program is using a series of labelling and morphological operations to extract the Lung volume intensity image from chest CT scans. It was tested with the data from the LIDC-IDRI (Lung Image Database Consortium) project and depends on ITK/cmake.

![screenshot](img/screenshot.png)

After the initial step of extracting the intensities of the lungs and airways the algorithm attempts to separate the two lungs and the airways.

![DICOM output files](img/DICOMOutput.png)

## Build

This tool depends on itk and cmake. As a simple way to document (and build) the project file a Dockerfile is provided. It should be sufficient to build the docker container once
```
docker build -t lungsegmentation -f Dockerfile .
```
and to run it using
```
> docker run --rm -it lungsegmentation
Usage : ./LungSegmentation
 System tags: 
   [ -v ] or [ -h ]
      = List options in short format
   [ -V ] or [ -H ]
      = List options in long format
   [ -vxml ] or [ -hxml ] or [ -exportXML ]
      = List options in xml format for BatchMake
   [ --xml ]
      = List options in xml format for Slicer
   [ -vgad ] or [ -hgad ] or [ -exportGAD ]
      = List options in Grid Application Description format
   [ -version ]
      = return the version number
   [ -date ]
      = return the cvs checkout date
 Command tags: 
   [ -n [ seriesname ] ]
      = Select series by series name (if more than one series is present).
   [ -b < labelfieldfilename > ]
      = Save the label field as a nifty file in the current directory
   [ -u < niftyfilename > ]
      = Save the corrected dataset as a nifty image to the current directory
   [ -V ]
      = Print more verbose output
 Command fields: 
   < indir > 
      = Directory with input DICOM image series.
   < outdir > 
      = Directory for output DICOM image series.
```
In order to provide the data for processing the docker call should also include a '-v' to mount your data directory inside the docker container. Here an example call that assumes you have a 'data/folder_with_dicom' directory in your current directory. The following call will save the output in the same folder.
```
docker run --rm -it -v `pwd`/data:/data lungsegmentation /data/folder_with_dicom /data/folder_with_dicom_segmented
```


### Debug build

Adjust the CMakeLists.txt file with the path for your itk version.
```
cmake -DCMAKE_BUILD_TYPE=Debug .
make
```

