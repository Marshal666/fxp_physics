param(
    [ValidateSet('Debug','Release')][string]$Configuration='Debug',
    [switch]$Portable,
    [switch]$Clang,
    [switch]$BuildOnly
)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$vswhere='C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe'
$installation=& $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (!$installation) { throw 'Visual Studio C++ tools are required.' }
$devcmd=Join-Path $installation 'Common7\Tools\VsDevCmd.bat'
$environment=& cmd.exe /d /s /c """$devcmd"" -no_logo -arch=x64 -host_arch=x64 >nul && set"
if ($LASTEXITCODE -ne 0) { throw 'Could not initialize Visual Studio environment.' }
foreach ($line in $environment) {
    if ($line -match '^([^=]+)=(.*)$') { [Environment]::SetEnvironmentVariable($matches[1],$matches[2],'Process') }
}
$cmake=Join-Path $installation 'Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe'
$ctest=Join-Path (Split-Path $cmake) 'ctest.exe'
$suffix=if($Portable){'-portable'}else{''}
if($Clang){$suffix+='-clang'}
$build=Join-Path $projectRoot ('build/windows-'+$Configuration.ToLowerInvariant()+$suffix)
$configure=@('-S',$projectRoot,'-B',$build,'-G','Ninja',"-DCMAKE_BUILD_TYPE=$Configuration")
if($Portable){$configure+=@('-DFXP_PORTABLE_ONLY=ON','-DFXP_DISABLE_SSE2=ON')}
if($Clang){$configure+=('-DCMAKE_CXX_COMPILER='+ (Join-Path $installation 'VC\Tools\Llvm\x64\bin\clang-cl.exe'))}
& $cmake @configure
if ($LASTEXITCODE -ne 0) { throw 'Configure failed.' }
& $cmake --build $build --parallel
if ($LASTEXITCODE -ne 0) { throw 'Build failed.' }
if (!$BuildOnly) {
    & $ctest --test-dir $build --output-on-failure
    if ($LASTEXITCODE -ne 0) { throw 'Tests failed.' }
}
