itk_module(ITKLevelSetsv4
  DEPENDS
    ITKCommon
    ITKStatistics
    ITKLabelMap
    ITKThresholding # itkBinaryThresholdImageFilter.h
    ITKDistanceMap
  TEST_DEPENDS
    ITKTestKernel
    ITKFastMarching
)
