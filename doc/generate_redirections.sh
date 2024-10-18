#!/bin/bash

# Check if the version number is provided as a command line argument
if [ $# -eq 0 ]; then
  echo "Please provide the version number as a command line argument."
  exit 1
fi

# Get the version number from the command line argument
version=$1

# Check if the 'latest' folder already exists and remove it if it does
if [ -d "latest" ]; then
  rm -rf latest
fi

# Recursively copy the directories and .html files from the specified version folder to 'latest'
rsync -av --include="*/" --include="*.html" --exclude="*" "$version/" latest/

# Save the path to all html files in 'latest' in the 'files' variable
files=$(find latest -type f -name "*.html")

# Loop over each HTML file and overwrite its contents
for file in $files; do
  redirect_path=${file//latest/$version}
  levels=$(grep -o '/' <<< "$redirect_path" | wc -l)
  prefix=""
  for ((i=0; i<$levels; i++)); do
    prefix="../$prefix"
  done
  echo '<meta http-equiv="refresh" content="0; URL='$prefix$redirect_path'"/>' > $file
done
