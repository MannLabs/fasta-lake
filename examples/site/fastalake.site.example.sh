#!/usr/bin/env bash
# Optional environment template for a local installation.
# Copy this file under a new local name and use paths on your own computer.
# The manifest runner itself takes data paths as explicit arguments.
# FASTALAKE_BIN_DIR is useful when all compiled executables have been collected
# in one directory. Otherwise add the Rust target/release directories to PATH.
# FASTALAKE_CONVERSION_SCRIPTS refers to optional independently maintained raw
# conversion and de novo scripts. It is unnecessary for the supported workflow.
# Optional study tests accept FL_JANKO_DIR or FASTALAKE_TEST_DATA. The synthetic
# demo and unit fixtures are shipped and do not require these study data.
# export FASTALAKE_BIN_DIR="/path/to/fastalake/bin"
# export FASTALAKE_CONVERSION_SCRIPTS="/path/to/optional/conversion/scripts"
# export FASTALAKE_TEST_DATA="/path/to/optional/study/fixtures"
