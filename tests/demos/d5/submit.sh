#!/bin/bash
set -e
cd $(dirname "$0")
sobjects prepare
sobjects submit -j
