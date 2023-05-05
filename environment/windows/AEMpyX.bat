@REM Store AEMpyX env vars and set to this conda env
@REM so other AEMpyX installs don't pollute the environment

@set "AEMPYX_ROOT=C:\Users\33642\AEMpyX"
@set "AEMPYX_DATA=C:\Users\33642\AEMpyX\data"

conda activate AEMpyX
mamba update --all

