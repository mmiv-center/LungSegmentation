#ifndef INCLUDE_MYTYPES_H_
#define INCLUDE_MYTYPES_H_

#include "itkImageRegionIterator.h"


// InputImageType::Pointer
// InputImageType::Pointer 
// ImageType::Pointer labelField

typedef signed short PixelType;
typedef itk::Image<PixelType, 3> ImageType;


#endif // INCLUDE_MYTYPES_H_