#! /usr/bin/bash

## 1) input directories
### check all file permissions are 766
declare -i PERMISSION
PERMISSION=8#766

# loop over directories
for direct in "S$@";
do
    for filename in $direct/*;
    do
	declare -i file_stat # declare file permission then set as octal
	file_stat=8#$(stat -c "%a" $filename)

	if [ "$file_stat" -gt "$PERMISSION" ]
	then
	    echo "Files with 766 permission: $filename"
	else
	    echo "Setting correct permissions: $filename"
	    chmod 766 $filename
	fi
    done
done
