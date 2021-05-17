#!/bin/bash

PID=Pat100
STUDY=phs000452

if [ ! -e /home/sbowler/biocore_lts/scott/phs000452/HRD/Pat11.out.seqz.gz ] ; then
	echo "not hardcoded"
else
	echo "hardcoded"
fi

if [ ! -e /home/sbowler/biocore_lts/scott/phs000452/HRD/"${PID}".out.seqz.gz ] ; then
	echo "not PID"
else
	echo "PID"
fi

if [ ! -e /home/sbowler/biocore_lts/scott/"${STUDY}"/HRD/"${PID}".out.seqz.gz ] ; then
	echo "not STUDY"
else
	echo "STUDY"
fi
