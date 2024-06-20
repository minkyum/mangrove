#!/bin/bash -l

site=/projectnb/modislc/users/mkmoon/Planet/rawImage/ord_test/$1

# for site in $PARENT/*
# do
	
	echo Processing site: $site

	download_paths=($(ls $site/data))

	for dir in $site/data/*
	do
		echo Processing: $dir
		for file_path in $dir/PSScene4Band/*
		do
			filename=$(basename $file_path)
			echo Processing: $filename;
			year=${filename:0:4}
			
			echo Move to $site/$year
			mkdir -p $site/$year
			mv $file_path $site/$year/.
		done
	done
#done
