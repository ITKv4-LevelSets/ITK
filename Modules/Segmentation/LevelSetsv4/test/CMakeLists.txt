itk_module_test()
set(ITK-LevelSetsv4Tests
itkLevelSetDomainPartitionBaseTest.cxx
)

CreateTestDriver(ITK-LevelSetsv4  "${ITK-LevelSetsv4-Test_LIBRARIES}" "${ITK-LevelSetsv4Tests}")

add_test(NAME itkLevelSetDomainPartitionBaseTest
      COMMAND ITK-LevelSets-v4TestDriver itkLevelSetDomainPartitionBaseTest)
