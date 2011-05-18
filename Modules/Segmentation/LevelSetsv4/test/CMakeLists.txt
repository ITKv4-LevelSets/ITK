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
# Sparse
itkSparseLevelSetBaseTest.cxx
itkWhitakerSparseLevelSetBaseTest.cxx
itkShiSparseLevelSetBaseTest.cxx
itkMalcolmSparseLevelSetBaseTest.cxx
itkBinaryImageToWhitakerSparseLevelSetAdaptorTest.cxx
itkBinaryImageToMalcolmSparseLevelSetAdaptorTest.cxx
itkBinaryImageToShiSparseLevelSetAdaptorTest.cxx
itkUpdateWhitakerSparseLevelSetTest.cxx
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
      COMMAND ITK-LevelSetsv4TestDriver itkSingleLevelSetDense2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png 25 25 10)

# Sparse
add_test(NAME itkSparseLevelSetBaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkSparseLevelSetBaseTest)
add_test(NAME itkWhitakerSparseLevelSetBaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkWhitakerSparseLevelSetBaseTest)
add_test(NAME itkShiSparseLevelSetBaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkShiSparseLevelSetBaseTest)
add_test(NAME itkMalcolmSparseLevelSetBaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkMalcolmSparseLevelSetBaseTest)
add_test(NAME itkBinaryImageToWhitakerSparseLevelSetAdaptorTest
      COMMAND ITK-LevelSetsv4TestDriver
      itkBinaryImageToWhitakerSparseLevelSetAdaptorTest
      ${ITK_DATA_ROOT}/Input/circle.png
      output_binary_whitaker_adaptor.mha
      )
add_test(NAME itkBinaryImageToMalcolmSparseLevelSetAdaptorTest
      COMMAND ITK-LevelSetsv4TestDriver
      itkBinaryImageToMalcolmSparseLevelSetAdaptorTest
      ${ITK_DATA_ROOT}/Input/circle.png
      output_binary_malcolm_adaptor.mha
      )
add_test(NAME itkBinaryImageToShiSparseLevelSetAdaptorTest
      COMMAND ITK-LevelSetsv4TestDriver
      itkBinaryImageToShiSparseLevelSetAdaptorTest
      ${ITK_DATA_ROOT}/Input/circle.png
      output_binary_shi_adaptor.mha
      )
add_test(NAME itkUpdateWhitakerSparseLevelSetTest0
      COMMAND ITK-LevelSetsv4TestDriver
      itkUpdateWhitakerSparseLevelSetTest
      0
      output_update_whitaker_sparse0.mha
      status_update_whitaker_sparse0.mha
      )
add_test(NAME itkUpdateWhitakerSparseLevelSetTest1
      COMMAND ITK-LevelSetsv4TestDriver
      itkUpdateWhitakerSparseLevelSetTest
      1
      output_update_whitaker_sparse1.mha
      status_update_whitaker_sparse1.mha
     )
add_test(NAME itkUpdateWhitakerSparseLevelSetTest2
      COMMAND ITK-LevelSetsv4TestDriver
      itkUpdateWhitakerSparseLevelSetTest
      2
      output_update_whitaker_sparse2.mha
      status_update_whitaker_sparse2.mha
    )
