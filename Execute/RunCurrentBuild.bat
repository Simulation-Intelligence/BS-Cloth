@echo off
REM Usage: script_name.bat <folder>

rmdir /s /q %1
mkdir %1

copy B-Spline-IPC\config.yaml %1\

mkdir %1\models
mkdir %1\quad_data
mkdir %1\seam

xcopy /s /e /i B-Spline-IPC\models\* %1\models\
xcopy /s /e /i B-Spline-IPC\quad_data\* %1\quad_data\
xcopy /s /e /i B-Spline-IPC\seam\* %1\seam\

xcopy /s /e /i bin\Release-windows-x86_64\B-Spline-IPC\* %1\

cd %1
.\B-Spline-IPC.exe

pause
