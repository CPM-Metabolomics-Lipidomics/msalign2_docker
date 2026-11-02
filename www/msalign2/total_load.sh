#!/bin/bash

ps aux | awk -vNPROCS=$(nproc) 'BEGIN{sum=0;}{if (NR > 1) sum=sum+$3;}END{print sum/NPROCS " %"}'
