itk_module_test()
set(ITKLevelSetsv4Tests
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
#itkTwoLevelSetMalcolm2DTest.cxx
itkTwoLevelSetShiSparse2DTest.cxx
itkSingleLevelSetSparseWithCurvature2DTest.cxx
itkSingleLevelSetSparseWithPropagation2DTest.cxx
itkSingleLevelSetSparseWithLaplacian2DTest.cxx
)

CreateTestDriver(ITKLevelSetsv4  "${ITKLevelSetsv4-Test_LIBRARIES}" "${ITKLevelSetsv4Tests}")

add_test(NAME itkLevelSetsv4HeaderTest
      COMMAND ITKLevelSetsv4TestDriver itkLevelSetsv4HeaderTest)
add_test(NAME itkLevelSetsv4DomainPartitionBaseTest
      COMMAND ITKLevelSetsv4TestDriver itkLevelSetDomainPartitionBaseTest)
add_test(NAME itkLevelSetsv4DomainMapImageFilterTest
      COMMAND ITKLevelSetsv4TestDriver itkLevelSetDomainMapImageFilterTest)
add_test(NAME itkLevelSetsv4ImageBaseTest
      COMMAND ITKLevelSetsv4TestDriver itkLevelSetImageBaseTest)
add_test(NAME itkMultiLevelSetImageTest
      COMMAND ITKLevelSetsv4TestDriver itkMultiLevelSetImageTest)
add_test(NAME itkLevelSetsv4ContainerBaseTest
      COMMAND ITKLevelSetsv4TestDriver itkLevelSetContainerBaseTest)
add_test(NAME itkMultiLevelSetsv4ChanAndVeseInternalTermTest
      COMMAND ITKLevelSetsv4TestDriver itkMultiLevelSetChanAndVeseInternalTermTest)
add_test(NAME itkMultiLevelsv4SetEvolutionTest
      COMMAND ITKLevelSetsv4TestDriver itkMultiLevelSetEvolutionTest)
add_test(NAME itkSingleLevelSetsv4Dense2DTest
      COMMAND ITKLevelSetsv4TestDriver itkSingleLevelSetDense2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png 25 25 25)
add_test(NAME itkTwoLevelSetsv4Dense2DTest
      COMMAND ITKLevelSetsv4TestDriver itkTwoLevelSetDense2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png 25 25 10)

# Sparse
add_test(NAME itkWhitakerSparseLevelSetsv4BaseTest
      COMMAND ITKLevelSetsv4TestDriver itkWhitakerSparseLevelSetBaseTest)
add_test(NAME itkShiSparseLevelSetsv4BaseTest
      COMMAND ITKLevelSetsv4TestDriver itkShiSparseLevelSetBaseTest)
add_test(NAME itkMalcolmSparseLevelSetsv4BaseTest
      COMMAND ITKLevelSetsv4TestDriver itkMalcolmSparseLevelSetBaseTest)
add_test(NAME itkBinaryImageToWhitakerSparseLevelSetsv4AdaptorTest
      COMMAND ITKLevelSetsv4TestDriver
      itkBinaryImageToWhitakerSparseLevelSetAdaptorTest
      ${ITK_DATA_ROOT}/Input/circle.png
      ${ITK_TEST_OUTPUT_DIR}/output_binary_whitaker_adaptor.mha
      )

add_test(NAME itkBinaryImageToMalcolmSparseLevelSetsv4AdaptorTest
      COMMAND ITKLevelSetsv4TestDriver
      itkBinaryImageToMalcolmSparseLevelSetAdaptorTest
      ${ITK_DATA_ROOT}/Input/circle.png
      ${ITK_TEST_OUTPUT_DIR}/output_binary_malcolm_adaptor.mha
      )
add_test(NAME itkBinaryImageToShiSparseLevelSetsv4AdaptorTest
      COMMAND ITKLevelSetsv4TestDriver
      itkBinaryImageToShiSparseLevelSetAdaptorTest
      ${ITK_DATA_ROOT}/Input/circle.png
      ${ITK_TEST_OUTPUT_DIR}/output_binary_shi_adaptor.mha
      )
add_test(NAME itkUpdateWhitakerSparseLevelSetsv4Test0
      COMMAND ITKLevelSetsv4TestDriver
      itkUpdateWhitakerSparseLevelSetTest
      0
      ${ITK_TEST_OUTPUT_DIR}/output_update_whitaker_sparse0.mha
      ${ITK_TEST_OUTPUT_DIR}/status_update_whitaker_sparse0.mha
      )
add_test(NAME itkUpdateWhitakerSparseLevelSetsv4Test1
      COMMAND ITKLevelSetsv4TestDriver
      itkUpdateWhitakerSparseLevelSetTest
      1
      ${ITK_TEST_OUTPUT_DIR}/output_update_whitaker_sparse1.mha
      ${ITK_TEST_OUTPUT_DIR}/status_update_whitaker_sparse1.mha
     )
add_test(NAME itkUpdateWhitakerSparseLevelSetsv4Test2
      COMMAND ITKLevelSetsv4TestDriver
      itkUpdateWhitakerSparseLevelSetTest
      2
      ${ITK_TEST_OUTPUT_DIR}/output_update_whitaker_sparse2.mha
      ${ITK_TEST_OUTPUT_DIR}/status_update_whitaker_sparse2.mha
    )

add_test(NAME itkUpdateShiSparseLevelSetsv4Test0
      COMMAND ITKLevelSetsv4TestDriver
      itkUpdateShiSparseLevelSetTest
      0
      ${ITK_TEST_OUTPUT_DIR}/output_update_shi_sparse0.mha
      )
add_test(NAME itkUpdateShiSparseLevelSetsv4Test1
      COMMAND ITKLevelSetsv4TestDriver
      itkUpdateShiSparseLevelSetTest
      1
      ${ITK_TEST_OUTPUT_DIR}/output_update_shi_sparse1.mha
     )
add_test(NAME itkUpdateShiSparseLevelSetsv4Test2
      COMMAND ITKLevelSetsv4TestDriver
      itkUpdateShiSparseLevelSetTest
      2
      ${ITK_TEST_OUTPUT_DIR}/output_update_shi_sparse2.mha
    )

add_test(NAME itkUpdateMalcolmSparseLevelSetsv4Test0
      COMMAND ITKLevelSetsv4TestDriver
      itkUpdateMalcolmSparseLevelSetTest
      0
      ${ITK_TEST_OUTPUT_DIR}/output_update_malcolm_sparse0.mha
      )
add_test(NAME itkUpdateMalcolmSparseLevelSetsv4Test1
      COMMAND ITKLevelSetsv4TestDriver
      itkUpdateMalcolmSparseLevelSetTest
      1
      ${ITK_TEST_OUTPUT_DIR}/output_update_malcolm_sparse1.mha
     )
add_test(NAME itkUpdateMalcolmSparseLevelSetsv4Test2
      COMMAND ITKLevelSetsv4TestDriver
      itkUpdateMalcolmSparseLevelSetTest
      2
      ${ITK_TEST_OUTPUT_DIR}/output_update_malcolm_sparse2.mha
    )
add_test(NAME itkSingleLevelSetsv4Sparse2DTest
      COMMAND ITKLevelSetsv4TestDriver itkSingleLevelSetSparse2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png
      1
      ${ITK_TEST_OUTPUT_DIR}/whiteSpot_output_sparse_single.mha
)

add_test(NAME itkSingleLevelSetsv4Malcolm2DTest
      COMMAND ITKLevelSetsv4TestDriver itkSingleLevelSetMalcolm2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png
      1
      ${ITK_TEST_OUTPUT_DIR}/whiteSpot_output_malcolm_single.mha
 )

add_test(NAME itkSingleLevelSetsv4Shi2DTest
      COMMAND ITKLevelSetsv4TestDriver itkSingleLevelSetShi2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png
      1
      ${ITK_TEST_OUTPUT_DIR}/whiteSpot_output_shi_single.mha)

add_test(NAME itkTwoLevelSetsv4Sparse2DTest
      COMMAND ITKLevelSetsv4TestDriver itkTwoLevelSetSparse2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png )

#add_test(NAME itkTwoLevelSetsv4Malcolm2DTest
#      COMMAND ITKLevelSetsv4TestDriver itkTwoLevelSetMalcolm2DTest
#      ${ITK_DATA_ROOT}/Input/whiteSpot.png )

add_test(NAME itkTwoLevelSetsv4ShiSparse2DTest
     COMMAND ITKLevelSetsv4TestDriver itkTwoLevelSetShiSparse2DTest
     ${ITK_DATA_ROOT}/Input/whiteSpot.png
     1
     ${ITK_TEST_OUTPUT_DIR}/whiteSpot_output_shi_two.mha
)

add_test(NAME itkSingleLevelSetsv4SparseWithCurvature2DTest
      COMMAND ITKLevelSetsv4TestDriver itkSingleLevelSetSparseWithCurvature2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png )

add_test(NAME itkSingleLevelSetsv4SparseWithPropagation2DTest
      COMMAND ITKLevelSetsv4TestDriver itkSingleLevelSetSparseWithPropagation2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png )

add_test(NAME itkSingleLevelSetsv4SparseWithLaplacian2DTest
      COMMAND ITKLevelSetsv4TestDriver itkSingleLevelSetSparseWithLaplacian2DTest
      ${ITK_DATA_ROOT}/Input/whiteSpot.png )
