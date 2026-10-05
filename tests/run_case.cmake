# Runs one solver invocation and checks its result. Invoked by CTest as
#   cmake -DSOLVER=<exe> -DINSTANCE=<file> [-DARGS="<flags>"] [-DEXPECT_STATUS=Optimal]
#         [-DEXPECT_PRIMAL=<value>] [-DEXPECT_STDERR_REGEX=<re>] [-DEXPECT_FAIL=ON]
#         [-DTIMEOUT=<sec>] -P run_case.cmake
#
# A passing case must exit with code 0, print a complete result block on
# stdout (a run that silently loses its output is a failure), report the
# expected status and primal value, and optionally match a regex on stderr.
# With EXPECT_FAIL=ON the solver must instead exit non-zero with an "Error:"
# line on stderr.

if(NOT SOLVER OR NOT INSTANCE)
    message(FATAL_ERROR "SOLVER and INSTANCE must be given")
endif()
if(NOT TIMEOUT)
    set(TIMEOUT 120)
endif()

separate_arguments(arg_list NATIVE_COMMAND "${ARGS}")
execute_process(
    COMMAND "${SOLVER}" "${INSTANCE}" ${arg_list}
    RESULT_VARIABLE rc
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err
    TIMEOUT ${TIMEOUT})

function(fail msg)
    message(STATUS "---- stdout ----\n${out}\n---- stderr ----\n${err}\n----------------")
    message(FATAL_ERROR "${msg}")
endfunction()

if(EXPECT_FAIL)
    if(rc EQUAL 0)
        fail("expected the solver to reject '${INSTANCE}' but it exited with 0")
    endif()
    if(NOT err MATCHES "Error:")
        fail("expected an 'Error:' message on stderr for '${INSTANCE}' (exit code ${rc})")
    endif()
    message(STATUS "rejected as expected (exit ${rc})")
    return()
endif()

if(NOT rc EQUAL 0)
    fail("solver exited with code ${rc}")
endif()
if(NOT out MATCHES "Status:[ \t]+([A-Za-z]+)")
    fail("no 'Status:' line on stdout (output lost or truncated)")
endif()
set(status "${CMAKE_MATCH_1}")
if(NOT out MATCHES "Primal:[ \t]+([^ \t\r\n]+)")
    fail("no 'Primal:' line on stdout")
endif()
set(primal "${CMAKE_MATCH_1}")

if(NOT EXPECT_STATUS)
    set(EXPECT_STATUS Optimal)
endif()
if(NOT status STREQUAL EXPECT_STATUS)
    fail("status is '${status}', expected '${EXPECT_STATUS}'")
endif()

if(DEFINED EXPECT_PRIMAL AND NOT EXPECT_PRIMAL STREQUAL "")
    set(ok FALSE)
    if(primal STREQUAL EXPECT_PRIMAL)
        set(ok TRUE)
    elseif(primal EQUAL EXPECT_PRIMAL)
        set(ok TRUE)
    endif()
    if(NOT ok)
        fail("primal is ${primal}, expected ${EXPECT_PRIMAL}")
    endif()
endif()

if(DEFINED EXPECT_STDERR_REGEX AND NOT EXPECT_STDERR_REGEX STREQUAL "")
    if(NOT err MATCHES "${EXPECT_STDERR_REGEX}")
        fail("stderr does not match '${EXPECT_STDERR_REGEX}'")
    endif()
endif()

message(STATUS "ok: status=${status} primal=${primal}")
