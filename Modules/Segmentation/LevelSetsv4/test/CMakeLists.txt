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
# Sparse
itkWhitakerSparseLevelSetBaseTest.cxx
itkShiSparseLevelSetBaseTest.cxx
itkMalcolmSparseLevelSetBaseTest.cxx
itkBinaryImageToWhitakerSparseLevelSetAdaptorTest.cxx
itkBinaryImageToMalcolmSparseLevelSetAdaptorTest.cxx
itkBinaryImageToShiSparseLevelSetAdaptorTest.cxx
itkUpdateWhitakerSparseLevelSetTest.cxx
itkUpdateShiSparseLevelSetTest.cxx
itkUpdateMalcolmSparseLevelSetTest.cxx
itkSingleLevelSetSparse2DTest.cxx
itkSingleLevelSetMalcolm2DTest.cxx
itkSingleLevelSetShi2DTest.cxx
itkTwoLevelSetSparse2DTest.cxx
itkTwoLevelSetMalcolm2DTest.cxx
itkTwoLevelSetShiSparse2DTest.cxx
itkSingleLevelSetSparseWithCurvature2DTest.cxx
itkSingleLevelSetSparseWithPropagation2DTest.cxx
itkSingleLevelSetSparseWithLaplacian2DTest.cxx
)

CreateTestDriver(ITK-LevelSetsv4  "${ITK-LevelSetsv4-Test_LIBRARIES}" "${ITK-LevelSetsv4Tests}")

add_test(NAME itkLevelSetsv4HeaderTest
      COMMAND ITK-LevelSetsv4TestDriver itkLevelSetsv4HeaderTest)
add_test(NAME itkLevelSetsv4DomainPartitionBaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkLevelSetDomainPartitionBaseTest)
add_test(NAME itkLevelSetsv4DomainMapImageFilterTest
      COMMAND ITK-LevelSetsv4TestDriver itkLevelSetDomainMapImageFilterTest)
add_test(NAME itkLevelSetsv4ImageBaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkLevelSetImageBaseTest)
add_test(NAME itkMultiLevelSetImageTest
      COMMAND ITK-LevelSetsv4TestDriver itkMultiLevelSetImageTest)
add_test(NAME itkLevelSetsv4ContainerBaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkLevelSetContainerBaseTest)
add_test(NAME itkMultiLevelSetsv4ChanAndVeseInternalTermTest
      COMMAND ITK-LevelSetsv4TestDriver itkMultiLevelSetChanAndVeseInternalTermTest)
add_test(NAME itkMultiLevelsv4SetEvolutionTest
      COMMAND ITK-LevelSetsv4TestDriver itkMultiLevelSetEvolutionTest)
add_test(NAME itkSingleLevelSetsv4Dense2DTest
      COMMAND ITK-LevelSetsv4TestDriver itkSingleLevelSetDense2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png 25 25 25)
add_test(NAME itkTwoLevelSetsv4Dense2DTest
      COMMAND ITK-LevelSetsv4TestDriver itkTwoLevelSetDense2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png 25 25 10)

# Sparse
add_test(NAME itkWhitakerSparseLevelSetsv4BaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkWhitakerSparseLevelSetBaseTest)
add_test(NAME itkShiSparseLevelSetsv4BaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkShiSparseLevelSetBaseTest)
add_test(NAME itkMalcolmSparseLevelSetsv4BaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkMalcolmSparseLevelSetBaseTest)
add_test(NAME itkBinaryImageToWhitakerSparseLevelSetsv4AdaptorTest
      COMMAND ITK-LevelSetsv4TestDriver
      itkBinaryImageToWhitakerSparseLevelSetAdaptorTest
      ${ITK_DATA_ROOT}/Input/circle.png
      output_binary_whitaker_adaptor.mha
      )
add_test(NAME itkBinaryImageToMalcolmSparseLevelSetsv4AdaptorTest
      COMMAND ITK-LevelSetsv4TestDriver
      itkBinaryImageToMalcolmSparseLevelSetAdaptorTest
      ${ITK_DATA_ROOT}/Input/circle.png
      output_binary_malcolm_adaptor.mha
      )
add_test(NAME itkBinaryImageToShiSparseLevelSetsv4AdaptorTest
      COMMAND ITK-LevelSetsv4TestDriver
      itkBinaryImageToShiSparseLevelSetAdaptorTest
      ${ITK_DATA_ROOT}/Input/circle.png
      output_binary_shi_adaptor.mha
      )
add_test(NAME itkUpdateWhitakerSparseLevelSetsv4Test0
      COMMAND ITK-LevelSetsv4TestDriver
      itkUpdateWhitakerSparseLevelSetTest
      0
      output_update_whitaker_sparse0.mha
      status_update_whitaker_sparse0.mha
      )
add_test(NAME itkUpdateWhitakerSparseLevelSetsv4Test1
      COMMAND ITK-LevelSetsv4TestDriver
      itkUpdateWhitakerSparseLevelSetTest
      1
      output_update_whitaker_sparse1.mha
      status_update_whitaker_sparse1.mha
     )
add_test(NAME itkUpdateWhitakerSparseLevelSetsv4Test2
      COMMAND ITK-LevelSetsv4TestDriver
      itkUpdateWhitakerSparseLevelSetTest
      2
      output_update_whitaker_sparse2.mha
      status_update_whitaker_sparse2.mha
    )

add_test(NAME itkUpdateShiSparseLevelSetsv4Test0
      COMMAND ITK-LevelSetsv4TestDriver
      itkUpdateShiSparseLevelSetTest
      0
      output_update_shi_sparse0.mha
      )
add_test(NAME itkUpdateShiSparseLevelSetsv4Test1
      COMMAND ITK-LevelSetsv4TestDriver
      itkUpdateShiSparseLevelSetTest
      1
      output_update_shi_sparse1.mha
     )
add_test(NAME itkUpdateShiSparseLevelSetsv4Test2
      COMMAND ITK-LevelSetsv4TestDriver
      itkUpdateShiSparseLevelSetTest
      2
      output_update_shi_sparse2.mha
    )

add_test(NAME itkUpdateMalcolmSparseLevelSetsv4Test0
      COMMAND ITK-LevelSetsv4TestDriver
      itkUpdateMalcolmSparseLevelSetTest
      0
      output_update_malcolm_sparse0.mha
      )
add_test(NAME itkUpdateMalcolmSparseLevelSetsv4Test1
      COMMAND ITK-LevelSetsv4TestDriver
      itkUpdateMalcolmSparseLevelSetTest
      1
      output_update_malcolm_sparse1.mha
     )
add_test(NAME itkUpdateMalcolmSparseLevelSetsv4Test2
      COMMAND ITK-LevelSetsv4TestDriver
      itkUpdateMalcolmSparseLevelSetTest
      2
      output_update_malcolm_sparse2.mha
    )
add_test(NAME itkSingleLevelSetsv4Sparse2DTest
      COMMAND ITK-LevelSetsv4TestDriver itkSingleLevelSetSparse2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png )

add_test(NAME itkSingleLevelSetsv4Malcolm2DTest
      COMMAND ITK-LevelSetsv4TestDriver itkSingleLevelSetMalcolm2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png )

add_test(NAME itkSingleLevelSetsv4Shi2DTest
      COMMAND ITK-LevelSetsv4TestDriver itkSingleLevelSetShi2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png )

add_test(NAME itkTwoLevelSetsv4Sparse2DTest
      COMMAND ITK-LevelSetsv4TestDriver itkTwoLevelSetSparse2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png 25 25 10)

add_test(NAME itkTwoLevelSetsv4Malcolm2DTest
      COMMAND ITK-LevelSetsv4TestDriver itkTwoLevelSetMalcolm2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png )

add_test(NAME itkTwoLevelSetsv4ShiSparse2DTest
      COMMAND ITK-LevelSetsv4TestDriver itkTwoLevelSetShiSparse2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png 25 25 10)

add_test(NAME itkSingleLevelSetsv4SparseWithCurvature2DTest
      COMMAND ITK-LevelSetsv4TestDriver itkSingleLevelSetSparseWithCurvature2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png )

add_test(NAME itkSingleLevelSetsv4SparseWithPropagation2DTest
      COMMAND ITK-LevelSetsv4TestDriver itkSingleLevelSetSparseWithPropagation2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png )

add_test(NAME itkSingleLevelSetsv4SparseWithLaplacian2DTest
      COMMAND ITK-LevelSetsv4TestDriver itkSingleLevelSetSparseWithLaplacian2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png )
