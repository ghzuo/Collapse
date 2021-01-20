###
# Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai, China.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@fudan.edu.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2017-09-01 15:54:30
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2021-01-20 18:06:53
###

INCLUDE_DIRECTORIES("../kit")

SET(LIBTAXSYS_SRC 
    taxsys.cpp taxsys.h
    reviseList.cpp  reviseList.h
    lineage.cpp lineage.h
    taxtree.cpp taxtree.h
    taxadb.cpp taxadb.h
    taxarank.cpp taxarank.h
)

SET(TAXHEADS taxsys.h reviseList.h lineage.h 
    taxtree.h taxadb.h taxarank.h)

SET(COLLAPSE_SRC ${KITHEADS} ${TAXHEADS}
    main.cpp
    collapse.cpp     collapse.h
    getLineage.cpp   getLineage.h
    getTaxaDB.cpp    getTaxaDB.h
    getLeafName.cpp  getLeafName.h
)

ADD_LIBRARY(taxsys STATIC ${LIBTAXSYS_SRC})

ADD_EXECUTABLE(cltree ${COLLAPSE_SRC})
TARGET_LINK_LIBRARIES(cltree kit taxsys z)

INSTALL(TARGETS cltree RUNTIME DESTINATION bin)
