itk_module_test()
set(ITK-LevelSetsv4Tests
itkLevelSetDomainPartitionBaseTest.cxx
itkLevelSetsv4HeaderTest.cxx
)

CreateTestDriver(ITK-LevelSetsv4  "${ITK-LevelSetsv4-Test_LIBRARIES}" "${ITK-LevelSetsv4Tests}")

add_test(NAME itkLevelSetsv4HeaderTest
      COMMAND ITK-LevelSetsv4TestDriver itkLevelSetsv4HeaderTest)
add_test(NAME itkLevelSetDomainPartitionBaseTest
      COMMAND ITK-LevelSetsv4TestDriver itkLevelSetDomainPartitionBaseTest)
