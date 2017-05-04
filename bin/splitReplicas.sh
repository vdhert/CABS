#! /bin/bash

awk -v fn=$1 '/\./{
 fname=sprintf("%s_%i",fn,$NF)
}{
 print $0 >fname;
}' $1
