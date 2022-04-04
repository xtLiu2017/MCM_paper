#!/bin/bash
a="mystring"
echo $a
for entry in "~/Documents/MCM_paper/data/dap_data_v4/motifs/WRKY_tnt"/*
do
	if grep colamp entry
		echo $entry
	done