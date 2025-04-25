@echo off
setlocal




REM Define versions and paths
set "R_VERSION=4.4.3"
set "R_URL=https://cran.r-project.org/bin/windows/base/R-4.4.3-win.exe"
set "INSTALL_PATH_R=C:\Program Files\R\R-%R_VERSION%"
set "R_EXEC=%INSTALL_PATH_R%\bin\R.exe"
set "RSTUDIO_EXEC=C:\Program Files\RStudio\rstudio.exe"
SET R_USERLIB=%USERPROFILE%\AppData\Local\R\win-library\4.4
SET R_LIBS_USER=%R_USERLIB%

:: Ensure the directory exists
if not exist "%R_USERLIB%" mkdir "%R_USERLIB%"



REM Check if R is installed
if not exist "%R_EXEC%" (
    echo R %R_VERSION% not found. Downloading and installing...
    powershell -Command "& {Invoke-WebRequest -Uri '%R_URL%' -OutFile 'R-%R_VERSION%-win.exe'}"
    start /wait R-%R_VERSION%-win.exe /SILENT /DIR="%INSTALL_PATH_R%"
    del /f /q "%~dp0R-%R_VERSION%-win.exe"
    if not exist "%R_EXEC%" (
        echo R installation failed. Exiting.
        exit /b 1
    )
    echo R %R_VERSION% installed successfully!
) else (
    echo R %R_VERSION% is already installed.
)

REM Run the R script automatically
echo Running the R script...

cd /d "%~dp0"  # Changes to the directory where the .bat file is located
start /wait "" "C:\Program Files\R\R-4.4.3\bin\x64\Rterm.exe" --no-save --quiet -f "%~ddlp.R"



:: Open RGui and preload the script
start "" "C:\Program Files\R\R-4.4.3\bin\x64\Rgui.exe"

pause
