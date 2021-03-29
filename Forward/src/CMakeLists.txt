# make libs of forward calculation

aux_source_directory(. DIR_LIB_SRCS_FORWARD)
add_library(forward ${DIR_LIB_SRCS_FORWARD})

install (TARGETS forward DESTINATION ${PROJECT_SOURCE_DIR}/lib)