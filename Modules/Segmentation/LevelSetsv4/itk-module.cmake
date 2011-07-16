itk_module(ITKLevelSetsv4
  DEPENDS
    ITKCommon
    ITKStatistics
    ITKThresholding # itkBinaryThresholdImageFilter.h
    ITKDistanceMap
  TEST_DEPENDS
    ITKTestKernel
    ITKFastMarching
)
