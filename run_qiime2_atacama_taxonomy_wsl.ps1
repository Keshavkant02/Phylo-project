param(
    [string]$DistroName = "qiime2-jammy",
    [string]$StageRoot = "E:\qiime2_atacama_taxonomy",
    [string]$WslRoot = "E:\WSL\qiime2-jammy",
    [string]$RepoRoot = "C:\Users\DELL\OneDrive\Documents\New project"
)

$ErrorActionPreference = "Stop"
if (Get-Variable PSNativeCommandUseErrorActionPreference -ErrorAction SilentlyContinue) {
    $PSNativeCommandUseErrorActionPreference = $false
}

$UbuntuRootfsUrl = "https://cloud-images.ubuntu.com/wsl/releases/22.04/current/ubuntu-jammy-wsl-amd64-wsl.rootfs.tar.gz"
$ClassifierUrl = "https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"
$ClassifierSha256 = "c08a1aa4d56b449b511f7215543a43249ae9c54b57491428a7e5548a62613616"
$QiimeEnvUrl = "https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml"

function New-RequiredDirectory {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path)) {
        New-Item -ItemType Directory -Force -Path $Path | Out-Null
    }
}

function Assert-FileExists {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path)) {
        throw "Required file missing: $Path"
    }
}

New-RequiredDirectory "$StageRoot\downloads"
New-RequiredDirectory "$StageRoot\inputs"
New-RequiredDirectory "$StageRoot\outputs"
New-RequiredDirectory $WslRoot

$RepSeqs = Join-Path $RepoRoot "tmp\atacama_qiime2_source\atacama-rep-seqs.qza"
$Table = Join-Path $RepoRoot "tmp\atacama_qiime2_source\atacama-table.qza"
$Metadata = Join-Path $RepoRoot "tmp\atacama_qiime2_source\sample_metadata.tsv"
Assert-FileExists $RepSeqs
Assert-FileExists $Table
Assert-FileExists $Metadata

Copy-Item -LiteralPath $RepSeqs -Destination "$StageRoot\inputs\atacama-rep-seqs.qza" -Force
Copy-Item -LiteralPath $Table -Destination "$StageRoot\inputs\atacama-table.qza" -Force
Copy-Item -LiteralPath $Metadata -Destination "$StageRoot\inputs\sample_metadata.tsv" -Force

$Rootfs = "$StageRoot\downloads\ubuntu-jammy-wsl-amd64-wsl.rootfs.tar.gz"
if (-not (Test-Path -LiteralPath $Rootfs)) {
    curl.exe -L -o $Rootfs $UbuntuRootfsUrl
}
Assert-FileExists $Rootfs

$InstalledDistros = & wsl.exe -l -q 2>$null
if ($InstalledDistros -notcontains $DistroName) {
    & wsl.exe --import $DistroName $WslRoot $Rootfs --version 2
}

$LinuxScript = @"
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y wget curl ca-certificates bzip2 coreutils

if [ ! -s /opt/miniforge3/bin/conda ] || [ ! -s /opt/miniforge3/etc/profile.d/conda.sh ] || ! /opt/miniforge3/bin/conda --version >/dev/null 2>&1; then
  rm -rf /opt/miniforge3
  cd /root
  rm -f Miniforge3-Linux-x86_64.sh
  wget -O Miniforge3-Linux-x86_64.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash Miniforge3-Linux-x86_64.sh -b -p /opt/miniforge3
fi

source /opt/miniforge3/etc/profile.d/conda.sh
conda --version

if ! conda env list | awk '{print `$1}' | grep -qx qiime2-amplicon-2024.10; then
  conda clean --packages --tarballs -y || true
  conda env create -n qiime2-amplicon-2024.10 --file "$QiimeEnvUrl"
fi

conda activate qiime2-amplicon-2024.10

mkdir -p /root/atacama-taxonomy-work/downloads /root/atacama-taxonomy-work/inputs /root/atacama-taxonomy-work/outputs
cp /mnt/e/qiime2_atacama_taxonomy/inputs/atacama-rep-seqs.qza /root/atacama-taxonomy-work/inputs/

cd /root/atacama-taxonomy-work/downloads
if [ ! -f silva-138-99-nb-classifier.qza ]; then
  wget -O silva-138-99-nb-classifier.qza "$ClassifierUrl"
fi
echo "$ClassifierSha256  silva-138-99-nb-classifier.qza" | sha256sum -c -

cd /root/atacama-taxonomy-work
qiime --version | tee /mnt/e/qiime2_atacama_taxonomy/outputs/qiime-version.txt
qiime info > /mnt/e/qiime2_atacama_taxonomy/outputs/qiime-info.txt

qiime feature-classifier classify-sklearn \
  --i-classifier downloads/silva-138-99-nb-classifier.qza \
  --i-reads inputs/atacama-rep-seqs.qza \
  --o-classification outputs/atacama-taxonomy.qza \
  --p-n-jobs 1

rm -rf outputs/exported-taxonomy
qiime tools export \
  --input-path outputs/atacama-taxonomy.qza \
  --output-path outputs/exported-taxonomy

cp outputs/atacama-taxonomy.qza /mnt/e/qiime2_atacama_taxonomy/outputs/
rm -rf /mnt/e/qiime2_atacama_taxonomy/outputs/exported-taxonomy
cp -r outputs/exported-taxonomy /mnt/e/qiime2_atacama_taxonomy/outputs/
"@

$ScriptPath = "$StageRoot\run_qiime2_taxonomy_inside_wsl.sh"
$LinuxScript = $LinuxScript -replace "`r`n", "`n"
[System.IO.File]::WriteAllText($ScriptPath, $LinuxScript, [System.Text.UTF8Encoding]::new($false))

$RunLog = "$StageRoot\outputs\run-log-$(Get-Date -Format 'yyyyMMdd-HHmmss').txt"
$PreviousErrorActionPreference = $ErrorActionPreference
$ErrorActionPreference = "Continue"
& wsl.exe -d $DistroName -- bash -lc "bash /mnt/e/qiime2_atacama_taxonomy/run_qiime2_taxonomy_inside_wsl.sh" *> $RunLog
$WslExitCode = $LASTEXITCODE
$ErrorActionPreference = $PreviousErrorActionPreference
if ($WslExitCode -ne 0) {
    throw "WSL QIIME run failed with exit code $WslExitCode. See $RunLog"
}
Copy-Item -LiteralPath $RunLog -Destination "$StageRoot\outputs\run-log.latest.txt" -Force

$TaxonomyPath = "$StageRoot\outputs\exported-taxonomy\taxonomy.tsv"
Assert-FileExists $TaxonomyPath

Copy-Item -LiteralPath $TaxonomyPath -Destination (Join-Path $RepoRoot "tmp\atacama_qiime2_source\atacama-taxonomy.tsv") -Force
Copy-Item -LiteralPath $TaxonomyPath -Destination (Join-Path $RepoRoot "soil_16s_class_cache\goal2_atacama_qiime_taxonomy.tsv") -Force

Push-Location $RepoRoot
try {
    python build_atacama_soil_asv_phylogeny_colab.py
    python verify_atacama_soil_asv_phylogeny_colab.py
}
finally {
    Pop-Location
}

Write-Host "QIIME taxonomy complete: $TaxonomyPath"
