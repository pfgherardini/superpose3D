

IF(WIN32)
	ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(WIN32)


SET(superpose_SRCS
	Amino.cpp
	Amino.h
	Chain.cpp
	Chain.h
	ChainWrapper.h
	ChainWrapper.cpp
	Comparison.cpp
	Comparison.h
	DefinitionCenter.cpp
	DefinitionCenter.h
	Job.cpp
	Job.h
	Match.cpp
	Match.h
	PdbParser.cpp
	PdbParser.h
	PdbResidueDictionary.cpp
	PdbResidueDictionary.h
	Quaternion.cpp
	Quaternion.h
	ResidueDefinition.cpp
	ResidueDefinition.h
	SimilarityMatrix.cpp
	SimilarityMatrix.h
	superpose.cpp
	superpose.h
	Superposition.cpp
	Superposition.h
	Tree.cpp
	Tree.h
	common.h
	string_ops.cpp
	main.cpp
)

SET(superpose_pdb_files_SRCS
	string_ops.cpp
	common.h
	superpose_pdb_files.cpp)

ADD_EXECUTABLE(superpose ${superpose_SRCS})
ADD_EXECUTABLE(superpose_pdb_files ${superpose_pdb_files_SRCS})





