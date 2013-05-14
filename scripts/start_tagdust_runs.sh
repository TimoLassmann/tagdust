#!/bin/sh

#  start_tagdust_runs.sh
#  tagdust2
#
#  Created by lassmann on 5/14/13.
#  Copyright (c) 2013 lassmann. All rights reserved.

#!/bin/sh

#
#	Set parameters below to fit your system
#

if [ $# -ne 1 ]
then
echo Usage is : $0 directory_name
fi

for file in $1/*
do
if [ -f $file ]; then
if [[ $file =~ tagdust_command.sh$ ]]; then
sh $file;
fi
fi
done
