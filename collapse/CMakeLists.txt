###
# Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai, China.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@fudan.edu.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2017-09-01 15:54:30
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2020-12-09 16:46:30
###

INCLUDE_DIRECTORIES("../kit")
LINK_LIBRARIES(libkit.a)
LINK_LIBRARIES(libtaxsys.a)

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
    collapse.cpp  collapse.h
)

SET(GETLNG_SRC ${KITHEADS} ${TAXHEADS} 
  getLineage.cpp   getLineage.h
  )

SET(GETNAME_SRC ${KITHEADS} ${TAXHEADS} 
  getLeafName.cpp
)

SET(GETDB_SRC ${KITHEADS} ${TAXHEADS} 
  getTaxaDB.cpp
)


ADD_LIBRARY(taxsys STATIC ${LIBTAXSYS_SRC})

ADD_EXECUTABLE(collapse ${COLLAPSE_SRC})
TARGET_LINK_LIBRARIES(collapse kit taxsys)

ADD_EXECUTABLE(getlng ${GETLNG_SRC})
TARGET_LINK_LIBRARIES(getlng kit taxsys)

ADD_EXECUTABLE(getname ${GETNAME_SRC})
TARGET_LINK_LIBRARIES(getname kit taxsys)

ADD_EXECUTABLE(getdb ${GETDB_SRC})
TARGET_LINK_LIBRARIES(getdb kit taxsys)

INSTALL(TARGETS collapse getlng getname getdb
  RUNTIME DESTINATION bin)
