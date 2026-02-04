#!/usr/bin/env bash
#SBATCH --job-name=mph_prep
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=csd795
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --export=ALL

set -eo pipefail

# Expect CONFIG_FILE passed in from submission script
: "${CONFIG_FILE:?CONFIG_FILE not set (pass via sbatch --export=ALL,CONFIG_FILE=/path/to/config.sh)}"

# shellcheck disable=SC1090
source "$CONFIG_FILE"

: "${MPH_DIR:?MPH_DIR must be set in config}"

# Repo directory
# NOTE: Slurm runs a spooled copy of this script on the compute node, so
# deriving paths from ${BASH_SOURCE[0]} would point into /cm/local/apps/slurm/var/spool/...
# The submission wrapper exports PIPELINE_REPO_DIR and also sets --chdir.
REPO_DIR="${PIPELINE_REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
if [[ ! -d "${REPO_DIR}" ]]; then
  echo "ERROR: REPO_DIR does not exist: ${REPO_DIR}" >&2
  echo "Hint: submit via bin/submit_mph_pipeline.sh (it exports PIPELINE_REPO_DIR)." >&2
  exit 2
fi
cd "${REPO_DIR}"

# Optional conda env activation (for R packages)
# NOTE: avoid sourcing ~/.bashrc in Slurm jobs (can break in non-interactive shells).
if [[ -n "${CONDA_ENV:-}" ]]; then
  CONDA_SH_PATH="${CONDA_SH:-}"

  if [[ -z "${CONDA_SH_PATH}" ]]; then
    # Common defaults
    if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
      CONDA_SH_PATH="$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
      CONDA_SH_PATH="$HOME/anaconda3/etc/profile.d/conda.sh"
    fi
  fi

  if [[ ! -f "${CONDA_SH_PATH}" ]]; then
    echo "ERROR: CONDA_ENV is set (${CONDA_ENV}) but conda.sh could not be found." >&2
    echo "Set CONDA_SH in your config (run 'conda info --base' on the login node)." >&2
    exit 2
  fi

  # shellcheck disable=SC1090
  source "${CONDA_SH_PATH}"
  conda activate "${CONDA_ENV}"

  echo "[env] Activated conda env: ${CONDA_ENV}"
  echo "[env] Rscript: $(command -v Rscript || echo NA)"
fi

# Hard requirements for this stage
command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found. Set CONDA_ENV in config or load an R module." >&2; exit 1; }
command -v "${PLINK_BIN:-plink}" >/dev/null 2>&1 || { echo "ERROR: plink not found (PLINK_BIN=${PLINK_BIN:-plink})." >&2; exit 1; }

echo "[prep] CONFIG_FILE=${CONFIG_FILE}"
echo "[prep] MPH_DIR=${MPH_DIR}"

mkdir -p "$MPH_DIR"/{pheno,genotypes,grm,reml,gsem,logs,tmp}

MAPPING_FILE="${MPH_DIR}/pheno/pheno_cohort_project_dict.csv"

# -------------------------------------------------------------------
# Step 0: genotype IDs are "king" (filter processed pheno early)
# -------------------------------------------------------------------
PHENO_SOURCE="${PHENO_SOURCE:-processed}"

# The genotype file used for MPH is whatever the user supplies in the config:
# - If SUBSET_GENO_BFILE is provided, that is the analysis genotype set.
# - Otherwise GENO_BFILE is the full set that will be subsetted.
KING_BFILE=""
if [[ -n "${SUBSET_GENO_BFILE:-}" ]]; then
  KING_BFILE="${SUBSET_GENO_BFILE}"
else
  : "${GENO_BFILE:?GENO_BFILE must be set in config (or set SUBSET_GENO_BFILE)}"
  KING_BFILE="${GENO_BFILE}"
fi

KING_FAM="${KING_BFILE}.fam"
if [[ ! -f "${KING_FAM}" ]]; then
  echo "ERROR: King genotype .fam not found: ${KING_FAM}" >&2
  exit 1
fi

KING_IDS="${MPH_DIR}/genotypes/king_genotyped_ids.txt"
awk '{print $1}' "${KING_FAM}" | sort -u > "${KING_IDS}"
echo "[prep] Genotyped ID universe: $(wc -l < "${KING_IDS}" | awk '{print $1}') (from ${KING_FAM})"

# If we're building phenotypes from the processed wide file, filter that file FIRST
# to genotyped IDs (rfid). This prevents phenotype-only IDs from entering cov/pheno outputs.
if [[ "${PHENO_SOURCE}" == "processed" ]]; then
  need_processed=0
  if [[ ! -f "${MAPPING_FILE}" ]]; then
    need_processed=1
  fi
  if [[ ! -f "${MPH_DIR}/pheno/combined_phenotype.csv" || ! -f "${MPH_DIR}/pheno/indicator_covariates.csv" ]]; then
    need_processed=1
  fi

  if [[ "${need_processed}" -eq 1 ]]; then
    : "${PROCESSED_FILE:?PROCESSED_FILE must be set in config for PHENO_SOURCE=processed}"
    if [[ ! -f "${PROCESSED_FILE}" ]]; then
      echo "ERROR: PROCESSED_FILE not found: ${PROCESSED_FILE}" >&2
      exit 1
    fi

    FILTERED_PROCESSED="${MPH_DIR}/tmp/processed_data_ready.genotyped.csv"

    if [[ ! -f "${FILTERED_PROCESSED}" || "${FILTERED_PROCESSED}" -ot "${PROCESSED_FILE}" || "${FILTERED_PROCESSED}" -ot "${KING_FAM}" ]]; then
      echo "[prep] Filtering processed file to genotyped IDs via column 'rfid' -> ${FILTERED_PROCESSED}"
      Rscript - <<RS
suppressPackageStartupMessages(library(data.table))
processed_file <- "${PROCESSED_FILE}"
ids_file       <- "${KING_IDS}"
out_file       <- "${FILTERED_PROCESSED}"

dt <- fread(processed_file)
if (!("rfid" %in% names(dt))) {
  stop("Processed file must contain column 'rfid': ", processed_file)
}
dt[, rfid := as.character(rfid)]

ids <- fread(ids_file, header = FALSE)[[1]]
ids <- as.character(ids)

dt2 <- dt[rfid %in% ids]
fwrite(dt2, out_file, sep = ",")
cat("[prep] Filtered processed rows: ", nrow(dt), " -> ", nrow(dt2), "\n", sep = "")
RS
    else
      echo "[prep] Using existing filtered processed file: ${FILTERED_PROCESSED}"
    fi

    # Use the filtered file for downstream mapping/stacking in this job
    PROCESSED_FILE="${FILTERED_PROCESSED}"
  fi
fi

# -------------------------------------------------------------------
# Step 1: build mapping file (if needed)
# -------------------------------------------------------------------
if [[ -f "$MAPPING_FILE" ]]; then
  echo "[prep] Found existing mapping file: $MAPPING_FILE"
else
  : "${PROCESSED_FILE:?PROCESSED_FILE must be set in config OR provide an existing $MAPPING_FILE}"
  if [[ ! -f "$PROCESSED_FILE" ]]; then
    echo "ERROR: PROCESSED_FILE not found: $PROCESSED_FILE" >&2
    exit 1
  fi

  include_arg=""
  exclude_arg=""

  if [[ -n "${INCLUDE_TRAITS_FILE:-}" ]]; then
    include_arg="--include_traits=@${INCLUDE_TRAITS_FILE}"
  elif [[ -n "${INCLUDE_TRAITS:-}" ]]; then
    include_arg="--include_traits=${INCLUDE_TRAITS}"
  fi

  if [[ -n "${EXCLUDE_TRAITS_FILE:-}" ]]; then
    exclude_arg="--exclude_traits=@${EXCLUDE_TRAITS_FILE}"
  elif [[ -n "${EXCLUDE_TRAITS:-}" ]]; then
    exclude_arg="--exclude_traits=${EXCLUDE_TRAITS}"
  fi

  echo "[prep] Creating mapping file -> $MAPPING_FILE"
  Rscript "${REPO_DIR}/scripts/make_pheno_cohort_project_dict.R"     --processed="${PROCESSED_FILE}"     --outdir="${MPH_DIR}/pheno"     --outname="pheno_cohort_project_dict.csv"     ${include_arg}     ${exclude_arg}
fi

# -------------------------------------------------------------------
# Step 2: stack phenotype + build cohort indicator covariates
# -------------------------------------------------------------------
PHENO_OUT="${MPH_DIR}/pheno/combined_phenotype.csv"
COV_OUT="${MPH_DIR}/pheno/indicator_covariates.csv"

if [[ -f "$PHENO_OUT" && -f "$COV_OUT" ]]; then
  echo "[prep] Found existing stacked phenotype + covariates. Skipping."
else
  echo "[prep] Stacking phenotypes and creating covariates (PHENO_SOURCE=${PHENO_SOURCE}) ..."

  if [[ "${PHENO_SOURCE}" == "processed" ]]; then
    : "${PROCESSED_FILE:?PROCESSED_FILE must be set for PHENO_SOURCE=processed}"
    if [[ ! -f "$PROCESSED_FILE" ]]; then
      echo "ERROR: PROCESSED_FILE not found: $PROCESSED_FILE" >&2
      exit 1
    fi

    Rscript "${REPO_DIR}/scripts/stack_pheno_create_cov.R"       --processed_file="${PROCESSED_FILE}"       --pheno_cohort_file="${MAPPING_FILE}"       --output_pheno_file="${PHENO_OUT}"       --output_cov_file="${COV_OUT}"

  elif [[ "${PHENO_SOURCE}" == "gcta" ]]; then
    : "${PHENO_DIR:?PHENO_DIR must be set for PHENO_SOURCE=gcta}"
    if [[ ! -d "$PHENO_DIR" ]]; then
      echo "ERROR: PHENO_DIR not found: $PHENO_DIR" >&2
      exit 1
    fi

    Rscript "${REPO_DIR}/scripts/stack_pheno_create_cov.R"       --pheno_dir="${PHENO_DIR}"       --pheno_cohort_file="${MAPPING_FILE}"       --output_pheno_file="${PHENO_OUT}"       --output_cov_file="${COV_OUT}"
  else
    echo "ERROR: Unknown PHENO_SOURCE='${PHENO_SOURCE}'. Use 'processed' or 'gcta'." >&2
    exit 1
  fi
fi

# -------------------------------------------------------------------
# Step 2b: enforce genotype-supported sample set (safeguard)
# -------------------------------------------------------------------
# Even if the user re-runs (and pheno/cov already exist), ensure we never carry
# phenotype-only IDs forward.
if [[ -f "${COV_OUT}" ]]; then
  echo "[prep] Filtering pheno/cov outputs to genotyped IDs (safeguard) ..."
  tmp_cov="${COV_OUT}.tmp"
  awk -F',' 'FNR==NR{ids[$1]=1; next} FNR==1{print; next} ($1 in ids){print}' "${KING_IDS}" "${COV_OUT}" > "${tmp_cov}"
  mv -f "${tmp_cov}" "${COV_OUT}"
fi
if [[ -f "${PHENO_OUT}" ]]; then
  tmp_pheno="${PHENO_OUT}.tmp"
  awk -F',' 'FNR==NR{ids[$1]=1; next} FNR==1{print; next} ($1 in ids){print}' "${KING_IDS}" "${PHENO_OUT}" > "${tmp_pheno}"
  mv -f "${tmp_pheno}" "${PHENO_OUT}"
fi

# -------------------------------------------------------------------
# Step 3: create ID list for genotype subsetting
# -------------------------------------------------------------------
IDS_1COL="${MPH_DIR}/genotypes/subset_ids.txt"
KEEP_FILE="${IDS_1COL}"

# Always regenerate: it is derived from the current covariate file.
echo "[prep] Creating subset_ids.txt from covariates..."
cut -d',' -f1 "${COV_OUT}" | tail -n +2 > "${IDS_1COL}"

if [[ "${PLINK_KEEP_FLAG:-"--keep-fam"}" == "--keep" ]]; then
  KEEP_FILE="${MPH_DIR}/genotypes/subset_keep.txt"
  echo "[prep] Creating subset_keep.txt (FID=IID=ID) for PLINK --keep ..."
  # naive default: FID=IID=ID
  awk '{print $1, $1}' "${IDS_1COL}" > "${KEEP_FILE}"
fi

# -------------------------------------------------------------------
# Step 4: subset genotypes to matching IDs (or use pre-subset prefix)
# -------------------------------------------------------------------
: "${PLINK_BIN:=plink}"
: "${CHR_RANGE:=1-20}"
: "${PLINK_KEEP_FLAG:=--keep-fam}"

SUBSET_PREFIX="${MPH_DIR}/genotypes/subset_genotypes"

# If subset already exists in the run directory, keep it.
if [[ -f "${SUBSET_PREFIX}.bed" && -f "${SUBSET_PREFIX}.bim" && -f "${SUBSET_PREFIX}.fam" ]]; then
  echo "[prep] Found existing subset genotypes in run dir: ${SUBSET_PREFIX}.*"

# Else, if user provided an already-subset genotype prefix elsewhere, link/copy it in.
elif [[ -n "${SUBSET_GENO_BFILE:-}" ]]; then
  if [[ ! -f "${SUBSET_GENO_BFILE}.bed" || ! -f "${SUBSET_GENO_BFILE}.bim" || ! -f "${SUBSET_GENO_BFILE}.fam" ]]; then
    echo "ERROR: SUBSET_GENO_BFILE is set but missing one of .bed/.bim/.fam: ${SUBSET_GENO_BFILE}.*" >&2
    exit 1
  fi

  echo "[prep] Using pre-subset genotypes: ${SUBSET_GENO_BFILE}.*"
  # Symlink if possible; fall back to copy if symlink fails.
  for ext in bed bim fam; do
    rm -f "${SUBSET_PREFIX}.${ext}" || true
    ln -s "${SUBSET_GENO_BFILE}.${ext}" "${SUBSET_PREFIX}.${ext}" 2>/dev/null || cp -f "${SUBSET_GENO_BFILE}.${ext}" "${SUBSET_PREFIX}.${ext}"
  done

else
  : "${GENO_BFILE:?GENO_BFILE must be set in config (or set SUBSET_GENO_BFILE)}"
  echo "[prep] Subsetting genotypes with PLINK..."
  "${PLINK_BIN}"     --bfile "${GENO_BFILE}"     ${PLINK_KEEP_FLAG} "${KEEP_FILE}"     --chr "${CHR_RANGE}"     --make-bed     --out "${SUBSET_PREFIX}"
fi

# -------------------------------------------------------------------
# Step 5: sanity check sample set match (MPH requirement)
# -------------------------------------------------------------------
if [[ ! -f "${SUBSET_PREFIX}.fam" ]]; then
  echo "ERROR: Missing subset .fam: ${SUBSET_PREFIX}.fam" >&2
  exit 1
fi

n_pheno=$(tail -n +2 "${COV_OUT}" | wc -l | awk '{print $1}')
n_fam=$(wc -l < "${SUBSET_PREFIX}.fam" | awk '{print $1}')

if [[ "${n_pheno}" -ne "${n_fam}" ]]; then
  echo "ERROR: Sample mismatch between phenotype IDs and genotype IDs." >&2
  echo "  phenotype IDs (indicator_covariates.csv): ${n_pheno}" >&2
  echo "  genotype IDs   (subset_genotypes.fam):    ${n_fam}" >&2
  echo "" >&2
  echo "MPH requires the GRM and phenotype/covariate files to have the SAME sample set." >&2
  echo "Fix the mismatch (or regenerate the subset) and re-run." >&2
  exit 1
fi

echo "[prep] Sample counts match (N=${n_fam})."
echo "[prep] Done."
