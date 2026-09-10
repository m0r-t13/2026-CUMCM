@echo off
setlocal
pushd "%~dp0"
python -X utf8 "%~dp0launch_viewer.py" %*
set "viewer_exit=%errorlevel%"
popd
if not "%viewer_exit%"=="0" pause
exit /b %viewer_exit%
