## generate step1 shell scripts
## for each sample and cellline, generate a shell script, eg. 
## you need to change the sample and cellline for different samples
bash step1.sh P1339-1 HEK > HEK_P1339-1.sh && bash HEK_P1339-1.sh

## generate step2 shell scripts and run them
## you need to provide a sample file contain all sample names, one per line without header, only one column
python3 step2.py
## then run the shell scripts, eg.
bash tk.1.sh

## get N17-1bp ins data

julia scripts/find_1bpINS.jl

## if you want to use get_indel_profile.py, you need to build it because it use cython
cd get_indel_profile
python3 setup.py build
