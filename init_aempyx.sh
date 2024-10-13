#! /bin/bash
export MY_ROOT=$HOME

export CONDA_ROOT=${HOME}/.Miniconda2024/

export AEMPYX_ENVI=AEM
export AEMPYX_ROOT=$MY_ROOT/AEMpyX_public/
export AEMPYX_DATA=$MY_ROOT/AEM_Data/

# # >>> conda initialize >>>
# # !! Contents within this block are managed by 'conda init' !!
# # __conda_setup="$('${HOME}.Miniconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
# __conda_setup="$('${CONDA_ROOT}/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
# if [ $? -eq 0 ]; then
#     eval "$__conda_setup"
# else
#     if [ -f ${CONDA_ROOT}/etc/profile.d/conda.sh ]; then
#         . ${CONDA_ROOT}/etc/profile.d/conda.sh
#     else
#         export PATH='${CONDA_ROOT}/bin:$PATH'
#     fi
# fi
# unset __conda_setup
# # <<< conda initialize <<<

cd ${AEMPYX_ROOT}
conda activate ${AEMPYX_ENVI}
