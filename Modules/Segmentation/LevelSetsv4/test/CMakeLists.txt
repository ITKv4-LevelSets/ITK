itk_module_test()
set(ITK-LevelSetsv4Tests
itkLevelSetsv4HeaderTest.cxx
itkLevelSetDomainPartitionBaseTest.cxx
itkLevelSetDomainMapImageFilterTest.cxx
itkLevelSetImageBaseTest.cxx
itkMultiLevelSetImageTest.cxx
itkLevelSetContainerBaseTest.cxx
itkMultiLevelSetChanAndVeseInternalTermTest.cxx
itkMultiLevelSetEvolutionTest.cxx
itkSingleLevelSetDense2DTest.cxx
itkTwoLevelSetDense2DTest.cxx
)

CreateTestDriver(ITK-LevelSetsv4  "${ITK-LevelSetsv4-Test_LIBRARIES}" "${ITK-LevelSetsv4Tests}")

add_test(NAME itkLevelSetsv4HeaderTest
      COMMAND ITK-LevelSetsv4TestDriver itkLevelSetsv4HeaderTest)
add_test(NAME itkLevelSetsv4DomainPartitionBaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkLevelSetDomainPartitionBaseTest)
add_test(NAME itkLevelSetsv4DomainMapImageFilterTest
      COMMAND ITK-LevelSetsv4TestDriver itkLevelSetDomainMapImageFilterTest)
add_test(NAME itkLevelSetImageBaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkLevelSetImageBaseTest)
add_test(NAME itkMultiLevelSetImageTest
      COMMAND ITK-LevelSetsv4TestDriver itkMultiLevelSetImageTest)
add_test(NAME itkLevelSetContainerBaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkLevelSetContainerBaseTest)
add_test(NAME itkMultiLevelSetChanAndVeseInternalTermTest
      COMMAND ITK-LevelSetsv4TestDriver itkMultiLevelSetChanAndVeseInternalTermTest)
add_test(NAME itkMultiLevelSetEvolutionTest
      COMMAND ITK-LevelSetsv4TestDriver itkMultiLevelSetEvolutionTest)
add_test(NAME itkSingleLevelSetDense2DTest
      COMMAND ITK-LevelSetsv4TestDriver itkSingleLevelSetDense2DTest ${ITK_DATA_ROOT}/Input/whiteSpot.png 25 25 10)
add_test(NAME itkTwoLevelSetDense2DTest
      COMMAND ITK-LevelSetsv4TestDriver itkTwoLevelSetDense2DTest ${ITK_DATA_ROOT}/Input/whiteSpot.png 25 25 10)
