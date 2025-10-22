# MGall_Analysis

Collection of scripts for analysis of Mycoplasma gallisepticum data: manifest generation, MD5 matching, GFF parsing & GWAS Manhattan plotting, ODE SEIRS modelling, parameter sweeps and sensitivity analysis.

Author: Dr. Alberto Torcello‑Requena — Environment and Sustainability Institute, University of Exeter

---

## Repository layout (complete / updated)

- Scripts/
  - MG_genome/
    - Manifest_files.py
      - Builds sample manifest CSV from raw read folders.
      - Reads `.md5` text files and replaces the `forward_file_md5` / `reverse_file_md5` columns in the template with the actual MD5 sum (file content).
      - Normalises sample name matching by stripping `.fastq.gz.md5` and `.fastq.md5` suffixes so it matches current FASTQ filenames (e.g. `.fastq`).
      - Produces `MGall_Manifest.csv` and `MGall_Manifest_with_MD5.csv`.
  - Python_scripts/
    - Manhattan_annotated.py
      - Produces Manhattan plots and annotates Bonferroni‑significant SNPs with nearby genes.
      - Parses GFF with Biopython / BCBio and uses `feature.location.strand` (avoids BiopythonDeprecationWarning).
      - Handles chromosome/name mismatches and reports useful debug messages; configurable distance window for annotation.
    - Rename_dictionary.py
      - Utilities to build replacement dictionaries for sample / folder names used by manifest creation.
  - MG_modelling/
    - Models/
      - params.py
        - Centralised default parameters and initial compartment counts. All modelling scripts import from this module.
      - SEIRS_Models.py
        - The SEIRS model implementation used by all modelling scripts (functions: `SEIRS_first_model`, `SEIRS_second_model`, etc.).
    - First_draft.py
      - Original ODE model file (kept for reference); now loads parameters from `Models.params`.
    - Paramenter_sensitivity_analysis_m.py
      - One‑at‑a‑time sensitivity analysis using the model in `Models.SEIRS_Models`.
      - Beautified horizontal bar (tornado) plot implemented with seaborn; updated to use `hue=...` + `legend=False` to avoid seaborn deprecation warnings.
      - Saves plots in headless mode and falls back to interactive display when an interactive backend is available.
    - Sweep_phi_theta_heatmap.py
      - 2D parameter sweep (phi_transmission × theta), computes peak and equilibrium metrics for high‑virulence strain and produces heatmaps and CSV table.
      - Reads default grids/time from `Models.params` when present.
    - Model_deeper_dive.py
      - Extended analysis: time series, Jacobian local stability, phase plane, theta sweeps and summary tables. Uses `Models.params` for defaults.

---

## Key design changes (what's different now)

- Centralised parameters:
  - All modelling and sweep scripts import model defaults and initial compartment counts from `Scripts/MG_modelling/Models/params.py`. Edit parameters in a single place.
  - `params.py` also optionally includes `phi_vals`, `theta_vals`, `t_max`, `t_steps` to control sweeps and time grids.

- MD5 / manifest handling:
  - Scripts now read `.md5` files (text) and replace the manifest MD5 columns with the actual MD5 value read from those files.
  - Matching between `.md5` filenames and FASTQ filenames normalises suffixes to handle both `.fastq` and legacy `.fastq.gz` names.

- GFF parsing:
  - Uses Biopython/BCBio parsing; accesses strand via `feature.location.strand` to avoid deprecation warnings.
  - `find_nearby_genes()` logic improved to compute distances and return meaningful defaults when no genes found.

- Plot backends & Qt issues:
  - Scripts detect/set matplotlib backend at startup:
    - Prefer interactive Qt backend (`Qt5Agg`) if available.
    - Fall back to non‑interactive `Agg` for headless servers / CI.
  - Sensible behavior: save plots automatically when `Agg` is used; show interactive windows only when a Qt/Tk backend is available.
  - If you need interactive plotting on Linux desktops, install a Qt binding:
    - pip: `pip install PyQt5` or `pip install PySide6`
    - conda: `conda install pyqt -c conda-forge`
  - For Jupyter interactive plots, install `ipympl` and enable `%matplotlib widget`.

- Seaborn deprecation: when using `palette` provide `hue` (we use `hue=param` and set `legend=False`) or supply a single `color=`.

---

## Requirements

Python 3.8+ recommended. Core packages:

- numpy, pandas, scipy, matplotlib, seaborn, biopython, bcbio‑gff, openpyxl

Install with pip:
```
pip install numpy pandas scipy matplotlib seaborn biopython bcbio-gff openpyxl
```

Optional (interactive plotting):
```
pip install PyQt5        # or PySide6
pip install ipympl        # for Jupyter widget backend
```

Conda example:
```
conda create -n MGall_Analysis python=3.8 numpy pandas scipy matplotlib seaborn biopython openpyxl -c conda-forge
conda activate MGall_Analysis
conda install pyqt -c conda-forge    # optional interactive backend
```

---

## Typical workflows / commands

1. Build manifest and replace MD5 columns with actual MD5 sums
```
python Scripts/MG_genome/Manifest_files.py
# Outputs: MGall_Manifest.csv and MGall_Manifest_with_MD5.csv
```

2. Produce Manhattan plot and annotate significant SNPs
```
python Scripts/Python_scripts/Manhattan_annotated.py --manhattan path/to/gwas.tsv --gff path/to/annotations.gff --output out.png
```
Notes:
- Ensure GWAS table contains columns `CHR`, `BP`, `P` or precomputed `LOG10P`.
- If annotation shows "Unknown (0bp)", check GFF chromosome names, coordinate system, and distance parameter.

3. Run baseline ODE model
```
python Scripts/MG_modelling/First_draft.py
```
(uses `Models.params` defaults)

4. Sensitivity analysis (tornado plot)
```
python Scripts/MG_modelling/Paramenter_sensitivity_analysis_m.py
# Saves sensitivity_tornado_plot.png
```

5. 2D sweep (phi × theta heatmap)
```
python Scripts/MG_modelling/Sweep_phi_theta_heatmap.py
# Saves ./Figures/phi_theta_heatmap_firstmodel.png and ./Tables/phi_theta_heatmap_results_firstmodel.csv
```

6. In‑depth model exploration
```
python Scripts/MG_modelling/Model_deeper_dive.py
# Saves multiple diagnostic figures (dynamics, phase plane, bifurcation-like plots)
```

---

## Debugging & tips

- If you see Qt plugin errors (wayland / xcb), the fastest solution for script runs is to set `MPLBACKEND=Agg` or let scripts fall back to `Agg`. For interactive plotting, install PyQt/PySide as noted above.
- Biopython warning: replace `.strand` with `.location.strand` when using SeqFeature objects.
- MD5 mapping:
  - Provide an Excel (`md5_values.xlsx`) with columns mapping MD5 filenames → MD5 value, or store individual `.md5` files in the reads folder. The manifest script supports both approaches.
- If a parameter change is needed across all analyses, edit `Scripts/MG_modelling/Models/params.py` once.
- For reproducibility, code seeds numpy RNG where stochastic behaviour may exist.

---

## Outputs

- MGall_Manifest.csv, MGall_Manifest_with_MD5.csv
- Manhattan plot PNGs (Manhattan_annotated.py)
- Sensitivity plots (sensitivity_tornado_plot.png)
- Sweep heatmaps (phi_theta_heatmap_firstmodel.png) and CSV results
- Diagnostic figures produced by Model_deeper_dive.py

---

## Contact & citation

For questions about the code or data layout contact:
Dr. Alberto Torcello‑Requena — alberto.torcello-requena@exeter.ac.uk

License: code provided "as is" (no license specified).

