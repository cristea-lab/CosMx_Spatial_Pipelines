#!/bin/bash

DIR=$PWD

if [ $# -lt 2 ]; then
    echo "ERROR: missing param- $0 requires a nanostring username and an atomix output folder name (which it will try to download into the current directory). "
    echo "Example: download_atomx.bash lentaing@mail.dfci.harvard.edu Alex_PCA84_PCA397_resegmented_11_04_2024_23_02_16_253 <optional: output path>"
    exit
fi

#sftp -r \"$user\"@na.export.atomx.nanostring.com:$2 .
if [ $# -eq 3 ]; then
    sftp -r -o user=$1 na.export.atomx.nanostring.com:$2 $3
else
    sftp -r -o user=$1 na.export.atomx.nanostring.com:$2 .
fi
