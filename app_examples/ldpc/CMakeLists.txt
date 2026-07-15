add_eff_app(ldpc
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c ldpc.c ldpc_utils.c ldpc_utils.h ${CMAKE_CURRENT_BINARY_DIR}/generator_80211.c ${CMAKE_CURRENT_BINARY_DIR}/parity_80211.c
)

eff_subtarget_compile_options(ldpc_fabric PRIVATE -DUSE_PRAGMAS)

add_custom_command(
    OUTPUT generator_80211.c
           parity_80211.c
    COMMAND python3 generate_80211_matrices.py -d ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/generate_80211_matrices.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_source_files_properties(
    ${CMAKE_CURRENT_BINARY_DIR}/generator_80211.c
    ${CMAKE_CURRENT_BINARY_DIR}/parity_80211.c
    PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}
)
