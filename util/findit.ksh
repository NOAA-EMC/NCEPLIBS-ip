#!/bin/ksh

#--------------------------------------------------------------------
# Determine whether OPS executable use certain IPLIB routines.
#--------------------------------------------------------------------

#set -x

iplib_files="gausslat_ gcdist_ gdsawt_ ijkgds_ \
             ipmerge2_ ipsector_ ipspaste_ \
             ipxetas_ ipxwafs2_ ipxwafs3_ ipxwafs_ \
             makgds_ movect_ polateg0_ polateg1_ polateg4_"

#--------------------------------------------------------------------
# Look for all directories with executables.  Assume they are named 'exec'.
# This can take a while, so save the directories to a log file the
# first time this script runs.  Subsequent runs of this script
# will then skip this step.  Also, the user can manualy edit
# this file to add/delete directories as desired.
#--------------------------------------------------------------------

if [[ ! -s ./exec.dir.log ]]; then
  for directory in /nwprod/*
  do
    find $directory -type d -name "exec" -print >> ./exec.dir.log
  done
fi

while read exec_dir
do
  echo $exec_dir
  for files in $exec_dir/*
  do
    if [[ -f $files ]]
    then
      echo CHECK $files
      nm -g $files > ./temp
      for ipfile in $iplib_files
      do
        grep -q "${ipfile}" ./temp
        status=$?
        if [[ $status == 0 ]];then
          echo FOUND ${ipfile} IN $files
        fi
      done  # loop over ip routines
      rm -f ./temp
    fi
  done
done < exec.dir.log

exit 0
