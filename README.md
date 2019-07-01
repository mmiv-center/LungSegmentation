# Basic Lung Segmentation using itk

This program is using a series of labelling and morphological operations to extract the Lung volume intensity image from chest CT scans. It was tested with the data from the LIDC-IDRI (Lung Image Database Consortium) project and depends on ITK/cmake.

![screenshot](img/screenshot.png)

After the initial step of extracting the intensities of the lungs and airways the algorithm attempts to separate the two lungs and the airways.

![DICOM output files](img/DICOMOutput.png)

### Debug build

Adjust the CMakeLists.txt file with the path for your itk version.
```
cmake -DCMAKE_BUILD_TYPE=Debug .
make
```

