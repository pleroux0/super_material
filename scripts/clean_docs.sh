#!/bin/sh

DIR="./dist/docs/"

if [[ -d "${DIR}" ]]
then
  rm -r "${DIR}"
fi
