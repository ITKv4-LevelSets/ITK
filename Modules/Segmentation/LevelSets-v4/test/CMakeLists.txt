itk_module_test()
set(ITK-LevelSetsTests-v4
itkLevelSetDomainPartitionBaseTest.cxx
)

CreateTestDriver(ITK-LevelSets-v4  "${ITK-LevelSets-Test_LIBRARIES}" "${ITK-LevelSetsTests-v4}")

add_test(NAME itkLevelSetDomainPartitionBaseTest
      COMMAND ITK-LevelSets-v4TestDriver itkLevelSetDomainPartitionBaseTest)
