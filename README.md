# Basic Lung Segmentation using itk

This program is using a series of labelling and morphological operations to extract the Lung volume from chest CT scans. It was tested with the data from the LIDC-IDRI (Lung Image Database Consortium) data.

![screenshot](img/screenshot.png)

### Debug build

```
cmake -DCMAKE_BUILD_TYPE=Debug .
make
```

