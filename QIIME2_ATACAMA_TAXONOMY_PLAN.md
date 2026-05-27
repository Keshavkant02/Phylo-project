# QIIME 2 Atacama Taxonomy Run Plan

## Goal

Produce a real QIIME 2 taxonomy artifact for the Atacama ASV representative sequences, export `taxonomy.tsv`, cache it in this repository, rebuild the Colab notebook, verify it, and push the updated notebook.

The student-facing Colab should remain browser-only. This QIIME 2 run is a one-time backend build step.

## Current Local State

- Workspace: `C:\Users\DELL\OneDrive\Documents\New project`
- E drive free space checked on 2026-05-26: about 53 GB.
- C drive free space checked on 2026-05-26: about 1 GB.
- Docker is not installed.
- WSL2 is available, but no Linux distribution is installed.
- Conda is installed on Windows, but QIIME 2 should not be installed natively on Windows.
- Real Atacama inputs already present:
  - `tmp\atacama_qiime2_source\atacama-rep-seqs.qza`
  - `tmp\atacama_qiime2_source\atacama-table.qza`
  - `tmp\atacama_qiime2_source\sample_metadata.tsv`

## Docs-Based Decisions

1. Use QIIME 2 `2024.10` Amplicon.
   - The Atacama artifacts we use are from QIIME 2 `2024.10`.
   - The `2024.10` Amplicon distribution includes `q2-feature-classifier`.
   - Official install command for Linux/WSL:

   ```bash
   conda env create -n qiime2-amplicon-2024.10 \
     --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
   ```

2. Use the pre-trained SILVA 138 Naive Bayes classifier for QIIME `2024.5-2026.1`.
   - Classifier URL:

   ```text
   https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza
   ```

   - Size confirmed by HTTP headers: `218245868` bytes.
   - SHA256 from QIIME Library page: `c08a1aa4d56b449b511f7215543a43249ae9c54b57491428a7e5548a62613616`.

3. Use `classify-sklearn` as the ideal method.
   - This gives a formal QIIME taxonomy artifact with confidence values.
   - The existing nearest-reference SILVA cache is only an interim fallback.

4. Stage the environment on `E:`.
   - `C:` is too full for a safe WSL/QIIME install.
   - Use WSL import to create a Linux distro directly under `E:\WSL\qiime2-jammy`.

## Proposed Disk Layout

```text
E:\qiime2_atacama_taxonomy\
  downloads\
    ubuntu-jammy-wsl-amd64-wsl.rootfs.tar.gz
    silva-138-99-nb-classifier.qza
  inputs\
    atacama-rep-seqs.qza
    atacama-table.qza
    sample_metadata.tsv
  outputs\
    atacama-taxonomy.qza
    exported-taxonomy\taxonomy.tsv
    qiime-info.txt
    run-log.txt

E:\WSL\qiime2-jammy\
  ext4.vhdx and WSL filesystem
```

## Execution Plan

### 1. Create staging directories on E

PowerShell:

```powershell
New-Item -ItemType Directory -Force E:\qiime2_atacama_taxonomy\downloads
New-Item -ItemType Directory -Force E:\qiime2_atacama_taxonomy\inputs
New-Item -ItemType Directory -Force E:\qiime2_atacama_taxonomy\outputs
New-Item -ItemType Directory -Force E:\WSL\qiime2-jammy
```

Copy the real Atacama inputs:

```powershell
Copy-Item -LiteralPath "C:\Users\DELL\OneDrive\Documents\New project\tmp\atacama_qiime2_source\atacama-rep-seqs.qza" -Destination E:\qiime2_atacama_taxonomy\inputs\atacama-rep-seqs.qza
Copy-Item -LiteralPath "C:\Users\DELL\OneDrive\Documents\New project\tmp\atacama_qiime2_source\atacama-table.qza" -Destination E:\qiime2_atacama_taxonomy\inputs\atacama-table.qza
Copy-Item -LiteralPath "C:\Users\DELL\OneDrive\Documents\New project\tmp\atacama_qiime2_source\sample_metadata.tsv" -Destination E:\qiime2_atacama_taxonomy\inputs\sample_metadata.tsv
```

### 2. Import Ubuntu WSL directly onto E

PowerShell:

```powershell
curl.exe -L -o E:\qiime2_atacama_taxonomy\downloads\ubuntu-jammy-wsl-amd64-wsl.rootfs.tar.gz https://cloud-images.ubuntu.com/wsl/releases/22.04/current/ubuntu-jammy-wsl-amd64-wsl.rootfs.tar.gz

wsl.exe --import qiime2-jammy E:\WSL\qiime2-jammy E:\qiime2_atacama_taxonomy\downloads\ubuntu-jammy-wsl-amd64-wsl.rootfs.tar.gz --version 2

wsl.exe -d qiime2-jammy -- bash -lc "cat /etc/os-release"
```

This avoids a normal Store/default WSL install landing on `C:`.

### 3. Install Miniforge and QIIME 2 inside WSL

PowerShell launching WSL:

```powershell
wsl.exe -d qiime2-jammy -- bash -lc "apt-get update && apt-get install -y wget curl ca-certificates bzip2 coreutils"
```

Inside WSL:

```bash
cd /root
wget -O Miniforge3-Linux-x86_64.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p /opt/miniforge3
source /opt/miniforge3/etc/profile.d/conda.sh
conda env create -n qiime2-amplicon-2024.10 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
conda activate qiime2-amplicon-2024.10
qiime --version
qiime info
```

Expected risk: this is the slowest step and may take 45-120 minutes depending on network and solver speed.

### 4. Download and verify the classifier

Inside WSL:

```bash
mkdir -p /root/atacama-taxonomy-work/downloads /root/atacama-taxonomy-work/inputs /root/atacama-taxonomy-work/outputs
cp /mnt/e/qiime2_atacama_taxonomy/inputs/atacama-rep-seqs.qza /root/atacama-taxonomy-work/inputs/

cd /root/atacama-taxonomy-work/downloads
wget -O silva-138-99-nb-classifier.qza https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza
echo "c08a1aa4d56b449b511f7215543a43249ae9c54b57491428a7e5548a62613616  silva-138-99-nb-classifier.qza" | sha256sum -c -
```

### 5. Classify Atacama ASVs

Inside WSL:

```bash
source /opt/miniforge3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10
cd /root/atacama-taxonomy-work

qiime feature-classifier classify-sklearn \
  --i-classifier downloads/silva-138-99-nb-classifier.qza \
  --i-reads inputs/atacama-rep-seqs.qza \
  --o-classification outputs/atacama-taxonomy.qza \
  --p-n-jobs 1

qiime tools export \
  --input-path outputs/atacama-taxonomy.qza \
  --output-path outputs/exported-taxonomy

cp outputs/atacama-taxonomy.qza /mnt/e/qiime2_atacama_taxonomy/outputs/
cp -r outputs/exported-taxonomy /mnt/e/qiime2_atacama_taxonomy/outputs/
qiime info > /mnt/e/qiime2_atacama_taxonomy/outputs/qiime-info.txt
```

### 6. Bring taxonomy back into the notebook cache

PowerShell from repo root:

```powershell
Copy-Item -LiteralPath E:\qiime2_atacama_taxonomy\outputs\exported-taxonomy\taxonomy.tsv -Destination .\tmp\atacama_qiime2_source\atacama-taxonomy.tsv
Copy-Item -LiteralPath E:\qiime2_atacama_taxonomy\outputs\exported-taxonomy\taxonomy.tsv -Destination .\soil_16s_class_cache\goal2_atacama_qiime_taxonomy.tsv
python build_atacama_soil_asv_phylogeny_colab.py
python verify_atacama_soil_asv_phylogeny_colab.py
```

Then inspect:

```powershell
Get-Content .\soil_16s_class_cache\goal2_atacama_feature_key.csv -TotalCount 15
```

### 7. Publish

Copy updated files into `tmp\phylo_project_publish`, verify there, commit, and push.

## Success Criteria

The QIIME step is successful only if all of these are true:

- `E:\qiime2_atacama_taxonomy\outputs\atacama-taxonomy.qza` exists and is a real QIIME artifact.
- `E:\qiime2_atacama_taxonomy\outputs\exported-taxonomy\taxonomy.tsv` exists.
- `soil_16s_class_cache\goal2_atacama_qiime_taxonomy.tsv` exists for published rebuilds.
- `taxonomy.tsv` contains columns such as `Feature ID`, `Taxon`, and `Confidence`.
- Feature IDs in `taxonomy.tsv` overlap the real Atacama representative sequence IDs.
- Rebuilt notebook verifier passes with 0 execution errors.
- Notebook manifest says taxonomy came from QIIME `classify-sklearn`, not the interim nearest-reference cache.
- Student-facing wording still says closest taxonomic match, not species proof.

## Execution Notes

- First run imported WSL onto `E:` and installed Ubuntu/Miniforge.
- A partial Miniforge install left zero-byte conda files; the wrapper now detects and replaces broken conda installs.
- Conda package cache corruption from the interrupted install was fixed by cleaning package/tarball caches before environment creation.
- QIIME 2 `2024.10.1` is now installed successfully inside WSL.
- The SILVA classifier downloaded and passed SHA256 verification.
- `classify-sklearn` was killed once under WSL's 3.8 GB memory limit with `--p-n-jobs 2`; the plan now uses `--p-n-jobs 1`.
- The successful run exported taxonomy and the builder now reports `Taxonomy loaded from goal2_atacama_qiime_taxonomy.tsv.`

## Fallbacks

If QIIME 2 install fails:

1. Try QIIME 2 `2025.10`/newer only if the classifier compatibility and command syntax are confirmed from docs.
2. Try QIIME `classify-consensus-vsearch` using the static SILVA 515F/806R sequence and taxonomy artifacts.
3. Keep the current nearest-reference SILVA cache as an honest fallback, explicitly labeled as not QIIME Naive Bayes.

## Source Links

- QIIME 2 2024.10 install overview: https://docs.qiime2.org/2024.10/install/
- QIIME 2 native/WSL install command: https://docs.qiime2.org/2024.10/install/native/
- QIIME 2 feature-classifier tutorial: https://docs.qiime2.org/2024.10/tutorials/feature-classifier/
- QIIME 2 taxonomy workflow overview: https://docs.qiime2.org/2024.10/tutorials/overview/
- QIIME 2 classifier resources: https://library.qiime2.org/data-resources
- Microsoft WSL install docs: https://learn.microsoft.com/en-us/windows/wsl/install
- Ubuntu WSL rootfs index: https://cloud-images.ubuntu.com/wsl/releases/22.04/current/
