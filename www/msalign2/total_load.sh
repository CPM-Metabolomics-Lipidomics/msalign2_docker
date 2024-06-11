#!/bin/bash

ps aux | awk -vNPROCS=$(cat /proc/cpuinfo  | grep CPU | wc -l|sed -e 's/ //g' ) 'BEGIN{sum=0;}{if (NR > 1) sum=sum+$3;}END{print sum/NPROCS " %"}'
