# Regression test suite.
#
# Every case runs the solver binary on an instance and checks the exit code,
# that a complete result block was printed, the reported status and objective,
# and optionally a pattern on stderr (see run_case.cmake).
#
#   ctest --test-dir build -C Release                 # everything (about 40 s in Release)
#   ctest --test-dir build -C Release -L core         # the 30 OR-Library instances only
#   ctest --test-dir build -C Release -L error        # malformed-input handling
#
# Labels: core, demo, fractional, error, features, regression.

set(_data ${CMAKE_SOURCE_DIR}/data)
set(_tdata ${CMAKE_CURRENT_SOURCE_DIR}/data)
set(_driver ${CMAKE_CURRENT_SOURCE_DIR}/run_case.cmake)

# scpsol_add_case(<name> <instance> <expected primal>
#                 [ARGS "<solver flags>"] [STATUS <Optimal|TimeLimit|NodeLimit>]
#                 [STDERR_REGEX <regex>] [TIMEOUT <seconds>] [LABEL <label>])
function(scpsol_add_case name instance expect_primal)
    cmake_parse_arguments(C "" "ARGS;STATUS;STDERR_REGEX;TIMEOUT;LABEL" "" ${ARGN})
    if(NOT C_STATUS)
        set(C_STATUS Optimal)
    endif()
    if(NOT C_TIMEOUT)
        set(C_TIMEOUT 300)
    endif()
    add_test(NAME ${name}
        COMMAND ${CMAKE_COMMAND}
            -DSOLVER=$<TARGET_FILE:scpsol>
            -DINSTANCE=${instance}
            "-DARGS=${C_ARGS}"
            -DEXPECT_STATUS=${C_STATUS}
            "-DEXPECT_PRIMAL=${expect_primal}"
            "-DEXPECT_STDERR_REGEX=${C_STDERR_REGEX}"
            -DTIMEOUT=${C_TIMEOUT}
            -P ${_driver})
    set_tests_properties(${name} PROPERTIES TIMEOUT ${C_TIMEOUT} LABELS "${C_LABEL}")
endfunction()

# A malformed instance must be rejected with exit code 1 and an "Error:" line.
function(scpsol_add_error_case name instance)
    add_test(NAME ${name}
        COMMAND ${CMAKE_COMMAND}
            -DSOLVER=$<TARGET_FILE:scpsol>
            -DINSTANCE=${instance}
            -DEXPECT_FAIL=ON
            -DTIMEOUT=60
            -P ${_driver})
    set_tests_properties(${name} PROPERTIES TIMEOUT 60 LABELS "error")
endfunction()

# Registers "<label>.<name>" for every "<name> <optimum>" pair in the list.
function(scpsol_add_cases label dir pairs)
    list(LENGTH pairs n)
    math(EXPR last "${n} - 1")
    foreach(i RANGE 0 ${last} 2)
        math(EXPR j "${i} + 1")
        list(GET pairs ${i} name)
        list(GET pairs ${j} opt)
        scpsol_add_case(${label}.${name} ${dir}/${name}.txt ${opt} LABEL ${label})
    endforeach()
endfunction()

# ---- OR-Library sets 4, 5, A, B with the optima from Beasley (1990) ----
scpsol_add_cases(core ${_data}
    "scp41;429;scp42;512;scp43;516;scp44;494;scp45;512;scp46;560;scp47;430;scp48;492;scp49;641;scp410;514;scp51;253;scp52;302;scp53;226;scp54;242;scp55;211;scp56;213;scp57;293;scp58;288;scp59;279;scp510;265;scpa1;253;scpa2;252;scpa3;232;scpa4;234;scpa5;236;scpb1;69;scpb2;76;scpb3;80;scpb4;79;scpb5;72")

# ---- Small demo instances; optima verified with an independent exact enumeration ----
# The scp_demo_tiny* and demo00 instances are solved entirely by preprocessing
# (the greedy solution is optimal) and were reported as NoSolution before #1 was fixed.
scpsol_add_cases(demo ${_data}
    "scp_demo00;6;scp_demo01;5;scp_demo02;6;scp_demo03;6;scp_demo04;4;scp_demo05;6;scp_demo06;5;scp_demo07;5;scp_demo08;6;scp_demo09;6;scp_demo_tiny00;1;scp_demo_tiny01;2;scp_demo_tiny02;1;scp_demo_tiny03;2;scp_demo_tiny04;2;scp_demo_tiny05;1;scp_demo_tiny06;1;scp_demo_tiny07;1;scp_demo_tiny08;2;scp_demo_tiny09;1;demo00;2;ex_balas1;13")

# ---- Non-integer costs (issue #1) ----
scpsol_add_cases(fractional ${_tdata}
    "frac_4x4;4.4;int_4x4;44;scp41_half;214.5;scp_demo_tiny00_half;0.5")

# ---- Malformed input must fail cleanly (issue #8) ----
foreach(bad bad_index_range bad_index_zero bad_negative_count bad_duplicate_index bad_truncated bad_dimensions)
    scpsol_add_error_case(error.${bad} ${_tdata}/${bad}.txt)
endforeach()
add_test(NAME error.missing_file
    COMMAND ${CMAKE_COMMAND} -DSOLVER=$<TARGET_FILE:scpsol> -DINSTANCE=${_tdata}/does_not_exist.txt
            -DEXPECT_FAIL=ON -DTIMEOUT=60 -P ${_driver})
set_tests_properties(error.missing_file PROPERTIES TIMEOUT 60 LABELS "error")

# ---- Optional features: each must still reach the optimum of scp52 ----
set(_f52 ${_data}/scp52.txt)
scpsol_add_case(features.cuts            ${_f52} 302 ARGS "--cut-frequency 0.5"      LABEL features)
scpsol_add_case(features.diving          ${_f52} 302 ARGS "--diving-frequency 1"     LABEL features)
scpsol_add_case(features.lagrangian      ${_f52} 302 ARGS "--lagrangian-frequency 1" LABEL features)
scpsol_add_case(features.rins            ${_f52} 302 ARGS "--rins-frequency 1"       LABEL features)
scpsol_add_case(features.decomposition   ${_f52} 302 ARGS "--decomposition force"    LABEL features)
scpsol_add_case(features.most_fractional ${_f52} 302 ARGS "--branch most_fractional" LABEL features)
scpsol_add_case(features.no_balas        ${_f52} 302 ARGS "--balas-frequency 0"      LABEL features)
scpsol_add_case(features.no_preprocess   ${_f52} 302 ARGS "--preprocess none"        LABEL features)
scpsol_add_case(features.ipm_root        ${_f52} 302 ARGS "--root-lp-solver ipm"     LABEL features)
scpsol_add_case(features.show_solution   ${_f52} 302 ARGS "--show-solution"          LABEL features)

# ---- Specific code paths ----
# Node limit: the run must stop with a complete, non-optimal result block.
scpsol_add_case(regression.node_limit ${_data}/scpb2.txt "" ARGS "--max-nodes 3" STATUS NodeLimit LABEL regression)
# Mid-search column removal (issue #5). In these random instances an incumbent
# found during branch-and-bound makes a column removable while the tree is
# still open, so the model is renumbered mid-search. All optima were verified
# with an independent exact enumeration.
#   - midbnb_removal_frac: removal after an integral node LP (fractional costs)
#   - midbnb_budget_pruning: integer budget rule fires after a node heuristic
#   - midbnb_requeue_a/b: a node heuristic improves the incumbent while the
#     node itself stays open, so it has to be remapped and re-enqueued
scpsol_add_case(regression.midbnb_column_removal ${_tdata}/midbnb_removal_frac.txt 15.7
    ARGS "--verbosity 3" STDERR_REGEX "Mid-BnB reduction" LABEL regression)
scpsol_add_case(regression.midbnb_budget_pruning ${_tdata}/midbnb_budget_pruning.txt 8
    ARGS "--verbosity 3" STDERR_REGEX "Mid-BnB budget pruning" LABEL regression)
scpsol_add_case(regression.midbnb_requeue_a ${_tdata}/midbnb_requeue_a.txt 20.4
    ARGS "--verbosity 3" STDERR_REGEX "Current node re-enqueued" LABEL regression)
scpsol_add_case(regression.midbnb_requeue_b ${_tdata}/midbnb_requeue_b.txt 16.6
    ARGS "--verbosity 3" STDERR_REGEX "Current node re-enqueued" LABEL regression)
