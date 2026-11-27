<#
.SYNOPSIS
    Build physicsforge wheels for one or more Python versions.

.DESCRIPTION
    For each requested Python version, creates an isolated virtual environment,
    builds an sdist and a wheel into ./dist, and optionally installs the wheel
    for a chosen version into ./.venv.

.PARAMETER PythonVersion
    Python versions to build for (resolved through the `py` launcher).

.PARAMETER Clean
    Remove the per-version virtual environments and intermediate build folders afterwards.

.PARAMETER PipInstall
    Python version whose freshly built wheel should be installed into ./.venv.

.PARAMETER Test
    Run the pytest suite against each freshly built wheel.

.EXAMPLE
    ./build_all.ps1 -PythonVersion 3.14 -Test -Clean
#>
param(
    [string[]]$PythonVersion = @("3.12", "3.13", "3.14"),
    [switch]$Clean,
    [switch]$Test,
    [string]$PipInstall
)

$ErrorActionPreference = "Stop"

$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$distDir = Join-Path $projectRoot "dist"
New-Item -ItemType Directory -Force -Path $distDir | Out-Null

foreach ($version in $PythonVersion)
{
    Write-Host "=== Building for Python $version ===" -ForegroundColor Cyan

    $pythonExe = & py -$version -c "import sys; print(sys.executable)" 2> $null
    if (-not $pythonExe)
    {
        Write-Warning "Python $version not found. Skipping."
        continue
    }

    $venvPath = Join-Path $projectRoot ".venv$version"
    if (Test-Path $venvPath)
    {
        Remove-Item -Recurse -Force $venvPath
    }
    & $pythonExe -m venv $venvPath
    $venvPython = Join-Path $venvPath "Scripts\python.exe"

    & $venvPython -m pip install --upgrade pip -q
    & $venvPython -m pip install build pytest -q

    & $venvPython -m build --sdist --wheel --outdir $distDir $projectRoot
    if ($LASTEXITCODE -ne 0)
    {
        Write-Error "Build failed for Python $version (exit code $LASTEXITCODE)"
    }

    if ($Test)
    {
        $tag = "cp" + ($version -replace '\.', '')
        $wheel = Get-ChildItem -Path $distDir -Filter "physicsforge-*-$tag-*.whl" |
                Sort-Object LastWriteTime -Descending | Select-Object -First 1
        & $venvPython -m pip install $wheel.FullName --force-reinstall -q
        # Run from outside the repo so the installed wheel is imported, not ./python.
        Push-Location $env:TEMP
        try
        {
            & $venvPython -m pytest (Join-Path $projectRoot "tests\python") -q -p no:cacheprovider -o addopts="" --rootdir $env:TEMP
            if ($LASTEXITCODE -ne 0)
            {
                Write-Error "Tests failed for Python $version"
            }
        }
        finally
        {
            Pop-Location
        }
    }

    if ($Clean)
    {
        Write-Host "Deleting build files..." -ForegroundColor Yellow
        Remove-Item -Recurse -Force -LiteralPath $venvPath
        foreach ($dir in @("build", "python/physicsforge.egg-info"))
        {
            $path = Join-Path $projectRoot $dir
            if (Test-Path $path)
            {
                Remove-Item -Recurse -Force -LiteralPath $path
            }
        }
    }

    Write-Host "Build for Python $version complete." -ForegroundColor Green
}

Write-Host "All builds complete. Wheels are in $distDir"

if ($PipInstall)
{
    if ($PythonVersion -notcontains $PipInstall)
    {
        Write-Warning "Python $PipInstall was not part of this build; nothing to install."
        return
    }

    $tag = "cp" + ($PipInstall -replace '\.', '')
    $wheel = Get-ChildItem -Path $distDir -Filter "physicsforge-*-$tag-*.whl" |
            Sort-Object LastWriteTime -Descending | Select-Object -First 1
    if (-not $wheel)
    {
        Write-Warning "No wheel found for Python $PipInstall in $distDir"
        return
    }

    $venvPath = Join-Path $projectRoot ".venv"
    if (-not (Test-Path $venvPath))
    {
        $pythonExe = & py -$PipInstall -c "import sys; print(sys.executable)"
        & $pythonExe -m venv $venvPath
    }
    Write-Host "Installing $( $wheel.Name ) into .venv..." -ForegroundColor Yellow
    & (Join-Path $venvPath "Scripts\python.exe") -m pip install $wheel.FullName --force-reinstall
    Write-Host "Installation complete." -ForegroundColor Cyan
}
