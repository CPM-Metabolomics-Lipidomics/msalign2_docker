#!/bin/bash

# start script with:
# 1: temp file to read which contains all the filenames for aligning
# 2: Master file which you use to align everything
# 3: mass error
# 4: timewindow [s]
# 5: number of features
# 6: costs per breakpoint
# 7: fitness
# 8: X option for msalign [s], 0 = end of run
# 9: Y option for msalign [s], 0 = end of run
# 10: start retention time [s]
# 11: end retention time [s], 0 = end of run

LINES="$( cat "$1" | awk '{ print $0 }' )"

# the commands below are run by a new bash, so quote everything that goes into them
# (file names may contain spaces)
QMASTER=$( printf '%q' "$2" )
QPARAMS=$( printf '%q ' "-e$3" "-f$5" "-c$6" "-d$7" "-X$8" "-Y$9" "-R${10},${11}" )

if awk -v t="$4" 'BEGIN { exit !(t + 0 == 0) }'; then   # also works for decimal seconds
  TIMEWINDOW=""
else
  TIMEWINDOW=$( printf '%q' "-l${4}" )
fi

MAXJOBS=6

function spawnjob()
{
  echo "$1" | bash
}

function clearToSpawn
{
  local JOBCOUNT="$( jobs -r | grep -c . )"
  if [ $JOBCOUNT -lt $MAXJOBS ]; then
    echo 1;
    return 1;
  fi
  echo 0;
  return 0;
}

JOBLIST=""

IFS="
"
for LINE in $LINES; do
  QLINE=$( printf '%q' "$LINE" )
  OUTPUTFILE=$( printf '%q' "${LINE%.*}_out.txt" )
  ALIGNFILE=$( printf '%q' "$LINE.alignment" )
  PLFFILE=$( printf '%q' "$LINE.alignment.plf" )
  # R splits its arguments on spaces, so the file name goes to wrap_mzXML.R in an environment variable
  COMMANDLIST+="./msalign2 -1$QLINE -2$QMASTER $QPARAMS $TIMEWINDOW > $OUTPUTFILE && ( echo 'Start creating aligned file' >> $OUTPUTFILE; MSALIGN_FILE=$QLINE R --no-restore --no-save --no-readline < ./wrap_mzXML.R; echo 'Done creating aligned file.' >> $OUTPUTFILE; ./FeaturesPlot.sh $ALIGNFILE $PLFFILE; echo 'Finished aligning!' >> $OUTPUTFILE )
  "
done

IFS="
"
for COMMAND in $COMMANDLIST; do
  while [ `clearToSpawn` -ne 1 ]; do
    sleep 1
  done
  spawnjob "$COMMAND" &
  LASTJOB=$!
  JOBLIST="$JOBLIST $LASTJOB"
done

IFS=" "
for JOB in $JOBLIST; do
  wait $JOB
  echo "Job $JOB exited with status $?"
done

echo "Done."


