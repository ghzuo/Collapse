###
# Copyright (c) 2018  Wenzhou Institute, University of Chinese Academy of Sciences.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2017-09-01 15:54:30
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2024-11-19 17:35:33
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

SET(LIBCLTREE_SRC 
  collapse.cpp       collapse.h
  queryLineage.cpp   queryLineage.h
  getTaxaDB.cpp      getTaxaDB.h
  getLeafName.cpp    getLeafName.h
  getRank.cpp        getRank.h
  searchLineage.cpp  searchLineage.h
  rooting.cpp        rooting.h
)

SET(CLTREEHEADS collapse.h queryLineage.h
  getTaxaDB.h getLeafName.h getRank.h 
  searchLineage.h)

SET(CLTREE_SRC ${KITHEADS} ${TAXHEADS} 
  ${CLTREEHEADS} main.cpp
)

ADD_LIBRARY(taxsys STATIC ${LIBTAXSYS_SRC})
ADD_LIBRARY(cltr   STATIC ${LIBCLTREE_SRC})

ADD_EXECUTABLE(cltree ${CLTREE_SRC})
TARGET_LINK_LIBRARIES(cltree cltr taxsys kit)

INSTALL(TARGETS cltree RUNTIME DESTINATION bin)
