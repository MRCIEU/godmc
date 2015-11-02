#!/bin/bash

set -e
source config

echo "Collecting log files"
tar czf uploads_logfiles.tgz log_files config resources/parameters

echo "Collecting results"
tar czvf uploads_results.tgz results config resources/parameters

echo "Successfully created results archives"