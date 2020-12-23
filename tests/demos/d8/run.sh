#!/bin/bash
set -e
cd $(dirname "$0")
. ./setenv.sh
sobjects prepare
sobjects run -j 
