#!/usr/bin/env bash
# install.sh — bootstrap the Plasmodium long-amplicon Speciation pipeline environment
#
# One-shot install script for collaborators. Installs vvg-box into the
# project directory, then populates it with every tool the pipeline needs.
# Idempotent: rerunning after a partial install picks up where it left off.
#
# Usage:
#   cd /path/to/Speciation
#   bash envs/install.sh
#
# Scope: all installs stay inside Speciation/. Nothing is written to
# /usr/local, /opt, ~/.zshrc, or any system directory.

set -euo pipefail

#---------------------------------------------------------------------------
# Locate the project root. This script assumes it's run from Speciation/
# or from envs/ inside the project.
#---------------------------------------------------------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"

# Sanity check
if [[ ! -f "$PROJECT_ROOT/INSTRUCTIONS.md" ]]; then
  echo "ERROR: could not locate project root. Expected INSTRUCTIONS.md at $PROJECT_ROOT"
  exit 1
fi

cd "$PROJECT_ROOT"
echo "==> Project root: $PROJECT_ROOT"
echo "==> Architecture: $(uname -s)-$(uname -m)"

# vvg-box and several of its helpers use `ln -sr` (GNU coreutils flag for
# relative symlinks), `readlink -f`, etc. macOS BSD ln/readlink reject these
# flags. If GNU coreutils is installed via Homebrew, prepend its gnubin so
# the install script and post-install helpers find GNU versions transparently.
if [[ "$(uname -s)" == "Darwin" ]]; then
  for gnubin in /opt/homebrew/opt/coreutils/libexec/gnubin /usr/local/opt/coreutils/libexec/gnubin; do
    if [[ -d "$gnubin" ]]; then
      export PATH="$gnubin:$PATH"
      echo "==> Prepending GNU coreutils to PATH: $gnubin"
      break
    fi
  done
  if ! ln --help 2>&1 | grep -q -- '--relative'; then
    echo "ERROR: macOS detected but GNU coreutils not found."
    echo "  Install with: brew install coreutils"
    echo "  Then rerun this script."
    exit 1
  fi
fi

#---------------------------------------------------------------------------
# Step 1: install vvg-box if not already present
#---------------------------------------------------------------------------
# Two pins matter here:
#   (a) MAMBA_ROOT_PREFIX. vvg-box's installer respects a pre-existing value.
#       If the user's shell has one set (e.g. from a prior `mamba init`), the
#       conda env lands there — outside this project — and collides with any
#       other project's vvg-box. Unset it so vvg-box defaults to project-local.
#   (b) The vvg-box commit. The default `bash <(curl .../main/install.sh)`
#       fetches HEAD, which currently sources `etc/functions` containing
#       `[[ -v VAR ]]` (bash 4.2+). macOS ships bash 3.2 and won't upgrade
#       (GPLv3 license). Pop-gen Indonesia bootstrapped its working vvg-box
#       at commit a0fb993 (2026-04-12), pre-dating that breaking change.
#       Pin to the same commit so this install reproduces what worked there.
unset MAMBA_ROOT_PREFIX
VVG_BOX_PIN="a0fb9939c97568eb8c3c7a37115463ff7c3607fe"

if [[ ! -f "$PROJECT_ROOT/vvg-box/bin/activate" && ! -f "$PROJECT_ROOT/activate" ]]; then
  echo "==> Installing vvg-box into $PROJECT_ROOT (project-local, pinned to ${VVG_BOX_PIN:0:7})"
  # Wipe any partially-bootstrapped vvg-box dir from a prior aborted run —
  # vvg-box's installer otherwise reuses an inconsistent state.
  if [[ -d "$PROJECT_ROOT/vvg-box" ]]; then
    echo "    (clearing partial $PROJECT_ROOT/vvg-box from a prior run)"
    rm -rf "$PROJECT_ROOT/vvg-box"
  fi
  # Fetch the installer at the pinned commit.
  TMP_INST="$(mktemp -t vvg_install.XXXXXX)"
  curl -fsSL \
    "https://raw.githubusercontent.com/vivaxgen/vvg-box/${VVG_BOX_PIN}/install.sh" \
    -o "$TMP_INST"
  # The pinned installer also clones HEAD of vvg-box by default. Patch the
  # `git clone --depth 1` line so the clone is at the same pinned commit;
  # otherwise we get the new (broken) repo even with the old installer.
  sed -i.bak \
    "s|git clone --depth 1 https://github.com/vivaxgen/vvg-box.git \${ENVS_DIR}/vvg-box|git clone https://github.com/vivaxgen/vvg-box.git \${ENVS_DIR}/vvg-box \&\& git -C \${ENVS_DIR}/vvg-box checkout ${VVG_BOX_PIN}|" \
    "$TMP_INST"
  if ! grep -q "checkout ${VVG_BOX_PIN}" "$TMP_INST"; then
    echo "ERROR: could not patch the vvg-box installer's clone line. Upstream"
    echo "       must have changed. Inspect $TMP_INST and update envs/install.sh."
    exit 1
  fi
  # Use bash explicitly (zsh refuses some vvg-box syntax). Redirect stdin
  # from /dev/null so any prompts auto-accept defaults.
  bash "$TMP_INST" < /dev/null
  rm -f "$TMP_INST" "$TMP_INST.bak"
else
  echo "==> vvg-box already installed — skipping bootstrap"
fi

#---------------------------------------------------------------------------
# Step 2: activate vvg-box for this install session
#---------------------------------------------------------------------------
# The activation path depends on what vvg-box's installer produced. Common
# shapes are ./activate or ./vvg-box/bin/activate. Probe both.
if [[ -f "$PROJECT_ROOT/activate" ]]; then
  ACTIVATE="$PROJECT_ROOT/activate"
elif [[ -f "$PROJECT_ROOT/vvg-box/bin/activate" ]]; then
  ACTIVATE="$PROJECT_ROOT/vvg-box/bin/activate"
elif [[ -f "$PROJECT_ROOT/bin/activate" ]]; then
  ACTIVATE="$PROJECT_ROOT/bin/activate"
else
  echo "ERROR: could not find vvg-box activation script. Expected one of:"
  echo "  $PROJECT_ROOT/activate"
  echo "  $PROJECT_ROOT/vvg-box/bin/activate"
  echo "  $PROJECT_ROOT/bin/activate"
  echo
  echo "If vvg-box put it somewhere else, update envs/install.sh with the"
  echo "right path, and also envs/activate.sh so collaborators can source it."
  exit 1
fi

# vvg-box's activate script (and its bashrc.d profiles) reference unset vars
# and call `history -r` on a non-existent file. Both trip `set -euo pipefail`.
# Relax both around the source, then re-enable.
# shellcheck disable=SC1090
set +eu
source "$ACTIVATE"
set -eu
echo "==> Activated vvg-box from $ACTIVATE"

# Record the activation path for future sessions
cat > "$PROJECT_ROOT/envs/activate.sh" <<EOF
#!/usr/bin/env bash
# Source this file to activate the pipeline environment.
# Generated by envs/install.sh on $(date).
source "$ACTIVATE"
EOF
chmod +x "$PROJECT_ROOT/envs/activate.sh"

#---------------------------------------------------------------------------
# Step 3: install Phase-1/2/3 tools into the vvg-box env
#---------------------------------------------------------------------------
# vvg-box is micromamba-based under the hood. Tools come from bioconda + conda-forge.
#
# Phase-1/2/3 tools (for replication of speciation_long):
#   - blast+:           blastn, makeblastdb
#   - seqkit:           FASTA subset, length tabulation
#   - mafft:            multiple-sequence alignment
#   - snakemake:        pipeline orchestration
#   - quarto:           report rendering (Quarto book)
#   - r-base + tidyverse + officer (for .docx reading in the Phase-5 report)
#   - python deps:      biopython, pandas, numpy
#
# Pinning strategy: latest stable on osx-arm64 / linux-64 unless a deviation
# from a prior HPC version is needed. Recorded in the software/versions notes.

# Prefer `mamba` if vvg-box exposed it; fall back to micromamba or conda.
if command -v mamba &>/dev/null; then
  INSTALLER="mamba install -y -c bioconda -c conda-forge"
elif command -v micromamba &>/dev/null; then
  INSTALLER="micromamba install -y -c bioconda -c conda-forge"
elif command -v conda &>/dev/null; then
  INSTALLER="conda install -y -c bioconda -c conda-forge"
else
  echo "ERROR: no mamba/micromamba/conda on PATH after activating vvg-box."
  echo "Inspect what vvg-box exposes and update this script accordingly."
  exit 1
fi

echo "==> Installing core CLI tools (blast+, seqkit, mafft, snakemake, quarto)"
$INSTALLER \
  blast \
  seqkit \
  mafft \
  snakemake-minimal \
  quarto \
  r-base=4.3

echo "==> Installing Python deps (biopython, pandas, numpy)"
$INSTALLER \
  biopython \
  pandas \
  numpy

echo "==> Installing R packages (tidyverse + .docx reading + plotting essentials)"
# r-officer reads .docx for the published-panel comparison chapter (officer
# is on conda-forge for osx-arm64). r-officedown / r-rvg are NOT on
# conda-forge osx-arm64 — installed via CRAN below. r-readxl is included
# in case any supplementary tables ship as .xlsx.
$INSTALLER \
  r-tidyverse \
  r-readxl \
  r-janitor \
  r-officer \
  r-data.table \
  r-viridis \
  r-svglite \
  r-biocmanager

echo "==> Installing Bioconductor Biostrings via BiocManager (if not on bioconda for this arch)"
# Biostrings is sometimes useful for FASTA / alignment manipulation in R; install
# only if it's not already pulled in by a tidy package.
Rscript -e '
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    BiocManager::install("Biostrings", ask = FALSE, update = FALSE)
  } else {
    message("Biostrings already installed — skipping")
  }
'

echo "==> Installing officedown + rvg via CRAN (no osx-arm64 build on conda-forge)"
# r-rvg is not packaged for conda-forge osx-arm64, which prevents
# r-officedown from solving via mamba. CRAN ships source builds that compile
# locally. Recorded as a deviation in the software/versions notes.
Rscript -e '
  pkgs <- c("rvg", "officedown")
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(to_install) > 0) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
  } else {
    message("officedown and rvg already installed — skipping")
  }
'

#---------------------------------------------------------------------------
# Step 4: Phase 5 primer-design tools (primer3, MFEprimer)
#---------------------------------------------------------------------------
# Added for Phase 5 Step 3 (primer-pair selection). primer3 is the CLI;
# primer3-py wraps the same thermodynamic core so we can compute per-primer
# hairpin ΔG and cross-pair heterodimer ΔG without a second CLI pass. Both
# resolve under vvg-box/opt/umamba/envs/vvg-box/bin/. MFEprimer is Step 5's
# in-silico PCR engine — installed here in the same pass so no re-solve is
# needed later. Off-target reference genomes stay Jacob-staged (embargo).
echo "==> Installing Phase-5 primer-design tools (primer3, MFEprimer)"
$INSTALLER \
  primer3 \
  mfeprimer

# primer3-py exposes primer3's thermodynamic functions from Python. Install
# via pip because bioconda's osx-arm64 build lags. Skip if already present.
echo "==> Installing primer3-py (pip)"
python -c 'import primer3' 2>/dev/null \
  && echo "    primer3-py already installed — skipping" \
  || pip install --quiet 'primer3-py>=2.0'

#---------------------------------------------------------------------------
# Step 5: write a lockfile so future installs are reproducible
#---------------------------------------------------------------------------
LOCKFILE="$PROJECT_ROOT/envs/environment.lock.yaml"
echo "==> Writing lockfile to $LOCKFILE"
if command -v micromamba &>/dev/null; then
  micromamba env export --explicit > "$LOCKFILE" 2>/dev/null || micromamba env export > "$LOCKFILE"
elif command -v mamba &>/dev/null; then
  mamba env export > "$LOCKFILE"
elif command -v conda &>/dev/null; then
  conda env export > "$LOCKFILE"
fi

#---------------------------------------------------------------------------
# Step 6: stage references from legacy into data/reference
#---------------------------------------------------------------------------
# The pipeline expects references in data/reference/. The legacy folder
# under scripts/legacy/speciation_long/data/ is read-only, so we copy
# rather than symlink (collaborators on systems with strict symlink
# semantics get the same files this way).
echo "==> Staging references from legacy into data/reference/"
LEGACY_DATA="$PROJECT_ROOT/scripts/legacy/speciation_long/data"
DATA_REF="$PROJECT_ROOT/data/reference"
if [[ -d "$LEGACY_DATA" ]]; then
  mkdir -p "$DATA_REF"
  for f in mit_all.fasta mit_all.target.fasta 18S_ref_db.fasta 18S_ref_db.target.fasta; do
    if [[ -f "$LEGACY_DATA/$f" && ! -f "$DATA_REF/$f" ]]; then
      cp -v "$LEGACY_DATA/$f" "$DATA_REF/$f"
    fi
  done
  # BLAST DBs — multiple files per DB. Copy the whole set if present.
  for dbprefix in homology_db homology_db_18S; do
    if compgen -G "$LEGACY_DATA/${dbprefix}.*" > /dev/null; then
      for f in "$LEGACY_DATA/${dbprefix}".*; do
        [[ -f "$DATA_REF/$(basename "$f")" ]] || cp -v "$f" "$DATA_REF/"
      done
    fi
  done
else
  echo "    (no legacy data dir found at $LEGACY_DATA — references will need to be staged manually)"
fi

#---------------------------------------------------------------------------
# Step 7: sanity check
#---------------------------------------------------------------------------
echo
echo "==> Sanity-check installed tools:"
for tool in blastn makeblastdb seqkit mafft snakemake quarto Rscript python; do
  if command -v "$tool" &>/dev/null; then
    ver="$("$tool" --version 2>&1 | head -1 || true)"
    printf "    %-14s  %s\n" "$tool" "$ver"
  else
    printf "    %-14s  NOT FOUND\n" "$tool"
  fi
done

echo
echo "==> Done. Activate in future sessions with:"
echo "     source envs/activate.sh"
echo
echo "==> Fill in the version table in the software/versions notes with the values printed above."
echo "==> Then proceed to the next pipeline step."
