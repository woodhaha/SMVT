@echo off
REM SMVT MD Environment Setup (Windows)
REM ===================================
echo ========================================
echo SMVT MD Environment Setup (Windows)
echo ========================================

REM Check conda
where conda >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: conda not found. Install Miniconda from:
    echo   https://docs.conda.io/en/latest/miniconda.html
    pause
    exit /b 1
)
echo [OK] conda found

REM Create environment
call conda env list | findstr "smvt-md" >nul
if %ERRORLEVEL% NEQ 0 (
    echo Creating conda environment 'smvt-md'...
    call conda create -n smvt-md python=3.11 -y
)
echo [OK] Environment ready

REM Install packages
echo Installing packages...
call conda activate smvt-md
call conda install -c conda-forge -y openmm pdbfixer mdtraj numpy pandas matplotlib seaborn scipy jupyter
call pip install openmmforcefields openff-toolkit rdkit

REM Verify
echo.
echo Verifying...
call python -c "import openmm; print(f'  openmm {openmm.__version__} OK')"
call python -c "import mdtraj; print(f'  mdtraj {mdtraj.__version__} OK')"
call python -c "import rdkit; print(f'  rdkit {rdkit.__version__} OK')"

echo.
echo ========================================
echo Setup complete!
echo Activate: conda activate smvt-md
echo Run MD:   python scripts\run_smvt_md.py --prepare ^&^& python scripts\run_smvt_md.py --compound HYDROMORPHONE --ns 100
echo ========================================
pause
