param(
    [string]$tb = "verilog/tb_top.sv",
    [string]$top = "",
    [string]$srcs = "verilog",
    [string]$out = "icarus/tb_top.vvp",
    [string]$vcd = "tb_top.vcd",
    [switch]$OpenWave,
    [string]$gtkwavePath = "gtkwave"
)

$ErrorActionPreference = 'Stop'

Set-Location -LiteralPath $PSScriptRoot

# Create output directory if not exists
if (-not (Test-Path "icarus")) {
    New-Item -ItemType Directory -Path "icarus" | Out-Null
}

$simTemp = (Resolve-Path "icarus").Path
$env:TMP = $simTemp
$env:TEMP = $simTemp

$tbPath = (Resolve-Path $tb).Path
if ([string]::IsNullOrWhiteSpace($top)) {
    $top = [System.IO.Path]::GetFileNameWithoutExtension($tbPath)
}

# Resolve source files (recursively if directory), excluding other testbenches.
# The selected testbench is passed explicitly as $tbPath, and -s binds the top.
if (Test-Path $srcs -PathType Container) {
    $srcFiles = Get-ChildItem -Path $srcs -Recurse -Include "*.sv", "*.v" -File |
        Where-Object {
            $_.FullName -ne $tbPath -and
            -not ($_.Name -match '^tb_.*\.(sv|v)$')
        } |
        Sort-Object FullName |
        Select-Object -ExpandProperty FullName
} else {
    $srcFiles = Resolve-Path $srcs | Select-Object -ExpandProperty Path
}

Write-Host "Top: $top"
Write-Host "Compiling: iverilog -g2012 -s $top -o $out $tb $srcFiles"
& iverilog -g2012 -s $top -o $out $tbPath $srcFiles
if ($LASTEXITCODE -ne 0) { Write-Error "iverilog failed (exit $LASTEXITCODE)"; exit $LASTEXITCODE }

Write-Host "Running: vvp $out"
& vvp $out
if ($LASTEXITCODE -ne 0) { Write-Error "vvp failed (exit $LASTEXITCODE)"; exit $LASTEXITCODE }

if ((Test-Path $vcd) -and $OpenWave) {
    Write-Host "VCD created: $vcd"
    try {
        Write-Host "Opening GTKWave (path: $gtkwavePath)"
        Start-Process -FilePath $gtkwavePath -ArgumentList $vcd -ErrorAction Stop
    } catch {
        Write-Warning "Failed to start GTKWave. You can open $vcd manually with GTKWave."
        Write-Host "If GTKWave is installed, provide its full path: .\run_sim.ps1 -gtkwavePath 'C:\Program Files\gtkwave\gtkwave.exe'"
    }
} elseif (-not (Test-Path $vcd)) {
    Write-Warning "VCD not found: $vcd"
} else {
    Write-Host "VCD created: $vcd"
}
