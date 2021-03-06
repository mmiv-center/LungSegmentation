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

## Generating training data for vessel tracing algorithms

The program FakeLungVolumes generates artificial vessel volumes. As those are generated by a noise smoothing process the data is artificial and only resembles vessels. Generation of such volumes is fast and can be tailored to resemble tissues of a given scale from micro-vascular trees to images resembling larger pulmonary vessel imaged with contrast CT.

The algorithm calculates the discrete intersection of two iso-surfaces from band-pass filtered white noise volumes. Based on the amount of band-pass filtering and the threshold for the iso-surface intersection detection blood-vessel like pattern appear.

To generate a 64 by 64 by 64 volume:
```
./FakeLungVolumes /tmp/output.nii 
```
To generate a higher resolution volume
```
./FakeLungVolumes -k 7 -t 0.0001 -f 1 -r 128x128x128 /tmp/output.nrrd
```
The generated output is artificially restricted to 12bit simulating common detector resolution. The intensity range is 0 to 4096. This does not correspond to HU but may be sufficient to generate test data for machine learning algorithms. 

To visualize the vessel generation process additionally to the vessels a four-class void space segmentation can be enabled. These void spaces are defined by the sign of the two band-pass filtered white noise volumes. Here an example:
```
./FakeLungVolumes -k 7 -t 0.0001 -w 0.0001 -f 1 -r 128x128x128 /tmp/output.nii
```

Here are all the options:
```
 Command tags: 
   [ -r [ resolution ] ]
      = Specify the resolution of the volume to be generated (in pixel as in 64x64x64).
   [ -k [ kernelSize ] ]
      = Specify the kernel size for the Gaussian in pixel (7).
   [ -i [ iterations ] ]
      = Specify the number of times the Gaussian kernels are applied (2).
   [ -t [ threshold ] ]
      = Specify the threshold for zero-crossing (0.0001).
   [ -f [ finalsmooth ] ]
      = Specify the kernel size of a smoothing with a Gaussian at the end of the process (0).
   [ -w [ voidspaces ] ]
      = Create void spaces with a given distance away from the lines. Default is that this option is not used.
   [ -s [ randomseed ] ]
      = Specify the value used for initialization of the random numbers (time based). The same value should produce the same fields.
   [ -z [ zero ] ]
      = Specify at what level the intersection should be performed (0).
   [ -f ]
      = Ignore existing files and force overwrite.
   [ -V ]
      = Print more verbose output
 Command fields: 
   < outfile > 
      = Exported file name.
```

![Fake vessel (gray) volume with 4 colored voids generated with ./FakeLungVolumes -k 7 -t 0.0001 -w 0.001 -r 128x128x128 /tmp/output.nii](https://github.com/mmiv-center/LungSegmentation/blob/master/img/FakeLungVoids.gif)


![Fake vessel volume generated with ./FakeLungVolumes -k 7 -t 0.0001 -f 1 -r 512x512x512 /tmp/output.nii](https://github.com/mmiv-center/LungSegmentation/blob/master/img/FakeVesselVolume.gif)

Computation time for the higher resolution volume is about 0.5seconds.

In order to explain the process geometrically the 'explain' sub-directory contains a website that performs this computation in JavaScript. By changing the number of iterations of the smoothing in 1-D, 2-D, and 3-D the different features can be visualized.

There are a number of extensions to this framework. One is that multiple vessel like structures can be generated with guaranteed properties. They will never intersect and keep a given distance from each other. Those hugging lines can be generated using the -z option to specify the crossing in conjunction with the seed option -s to make sure multiple runs of the program generate the same noise pattern. Given the same seed but different crossing levels vessel structures emerge. Here an example of the resulting 3 vessel like structures using Iso-Surface displays. This data was generated using three calls to the program:
```
./FakeLungVolumes -s 42 -z 0 -k 7 -t 0.0001 -f 0.4 -r 192x192x192 /tmp/output.nii
./FakeLungVolumes -s 42 -z 0.002 -k 7 -t 0.0001 -f 0.4 -r 192x192x192 /tmp/output2.nii
./FakeLungVolumes -s 42 -z -0.002 -k 7 -t 0.0001 -f 0.4 -r 192x192x192 /tmp/output3.nii
```

![Fake vessel volume with 3 sets of vessels generated with ./FakeLungVolumes -s 42 -z 0 -k 7 -t 0.0001 -f 0.4 -r 192x192x192 /tmp/output.nii; ./FakeLungVolumes -s 42 -z 0.002 -k 7 -t 0.0001 -f 0.4 -r 192x192x192 /tmp/output2.nii; ./FakeLungVolumes -s 42 -z -0.002 -k 7 -t 0.0001 -f 0.4 -r 192x192x192 /tmp/output3.nii](https://github.com/mmiv-center/LungSegmentation/blob/master/img/3setsNeverIntersectingIsoSurf.gif)

At lower resolution and in 2-D cross-section:

![3 sets of fake vessels in a cross-section.](https://github.com/mmiv-center/LungSegmentation/blob/master/img/3setsNeverIntersecting.gif)

As a final example here a closeup of a de novo in silico complex tissue.

![3 sets of fake vessels in a cross-section.](https://github.com/mmiv-center/LungSegmentation/blob/master/img/3setWithVoids.png)