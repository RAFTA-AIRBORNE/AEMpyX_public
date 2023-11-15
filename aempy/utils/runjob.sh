#!/bin/bash
myenv="AEMpyX"
pwd
if [ -z "$CONDA_DEFAULT_ENV" ]
then
echo "$myenv environment is not set!"
   . ~/bin/conda_init
   conda activate $myenv
else
   echo "conda active environment is:"
   echo $CONDA_DEFAULT_ENV
fi


if [ "$CONDA_DEFAULT_ENV" = $myenv ]; then
    echo "$myenv environment is set!"
else
    echo "$myenv environment is not set!"
    exit
fi


for f in $@
do
   echo "job file is $f"
   nohup python3 -u $f.py > $f.log &
   sleep 1
done

myjobs | grep python
