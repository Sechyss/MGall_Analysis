# MGall_Analysis

Collection of scripts used for analysis of Mycoplasma gallisepticum genomes and associated analyses (manifest generation, GFF parsing, GWAS Manhattan plotting, ODE modelling and sensitivity analysis).

Author: Dr. Alberto Torcello-Requena — Environment and Sustainability Institute, University of Exeter

---

## Repository layout (key scripts)

- Scripts/MG_genome/
  - Manifest_files.py
    - Build sample manifest CSV from raw folders.
    - Replace `.md5` filenames in the template with actual MD5 sums (reads `.md5` file contents).
    - Produces `MGall_Manifest.csv` and `MGall_Manifest_with_MD5.csv`.

- Scripts/Python_scripts/
  - Manhattan_annotated.py
    - Builds Manhattan plots from GWAS data.
    - Parses GFF files using `bcbio-gff` / Biopython (uses `BCBio.GFF`) to annotate significant SNPs.
    - Applies Bonferroni threshold and annotates nearby genes.
  - Rename_dictionary.py
    - Utilities to build replacement dictionaries for sample / folder names.

- Scripts/MG_modelling/
  - First_draft.py
    - ODE compartmental model implementation (integration via `scipy.integrate.odeint`).
  - Paramenter_sensitivity_analysis_m.py
    - One-at-a-time parameter sensitivity analysis and plotting.
    - Beautified sensitivity plot (uses seaborn, updated to avoid deprecated palette usage by assigning `hue` and `legend=False`).

---

## Requirements

Recommended Python >= 3.8. Example packages:

- numpy
- pandas
- scipy
- matplotlib
- seaborn
- biopython
- bcbio-gff
- openpyxl (for reading Excel MD5 list)

Install with pip:
```
pip install numpy pandas scipy matplotlib seaborn biopython bcbio-gff openpyxl
```

Or create a conda environment:
```
conda create -n MGall_Analysis python=3.8 numpy pandas scipy matplotlib seaborn biopython openpyxl -c conda-forge
```

---

## Usage notes and examples

1. Manifest generation and MD5 replacement
   - Edit `base` paths and file paths at the top of `Scripts/MG_genome/Manifest_files.py` to match your filesystem.
   - The script expects an Excel with MD5 mapping (`md5_values.xlsx`) and a template TSV (`fastq2_template_*.tsv`).
   - The script normalizes MD5 keys by removing `.fastq.gz.md5` and `.fastq.md5` suffixes so they match current fastq filenames (`.fastq`).
   - Run:
     ```
     python Scripts/MG_genome/Manifest_files.py
     ```
   - Output: `MGall_Manifest.csv` and `MGall_Manifest_with_MD5.csv`.

2. Manhattan plotting and GFF annotation
   - Ensure input GWAS file contains required columns (e.g., `CHR`, `BP`, `P`, `LOG10P`) or adjust `Manhattan_annotated.py` accordingly.
   - Provide a GFF file; the parser uses `feature.location` (Biopython API) and reads qualifiers for `gene` / `product`.
   - Run:
     ```
     python Scripts/Python_scripts/Manhattan_annotated.py --manhattan path/to/file.tsv --gff path/to/annotations.gff --output out.png
     ```
   - Notes:
     - Bonferroni threshold computed as `-log10(alpha / n_tests)`.
     - If no nearby genes are found the script returns a default annotation; increase `distance` if needed.
     - Uses `BCBio.GFF` to avoid fragile text parsing.

3. ODE model and sensitivity analysis
   - `First_draft.py` contains the model. `Paramenter_sensitivity_analysis_m.py` runs sensitivity experiments and saves plots.
   - Non-interactive runs should set matplotlib backend to `Agg` (already set in `Paramenter_sensitivity_analysis_m.py`) to avoid Qt plugin errors:
     - If you see `qt.qpa.plugin: Could not find the Qt platform plugin "wayland"`, either use `matplotlib.use('Agg')` or install Qt packages (for interactive sessions).
   - Run:
     ```
     python Scripts/MG_modelling/Paramenter_sensitivity_analysis_m.py
     ```

---

## Implementation & debugging tips

- GFF parsing: use `.location.strand` instead of `.strand` to avoid Biopython deprecation warnings.
- MD5 matching:
  - Strip both `.fastq.gz.md5` and `.fastq.md5` when deriving base sample names from `.md5` filenames to match `.fastq` files.
  - To replace the MD5 filename columns with the content of `.md5` files, build a dictionary mapping filename → md5_sum and then `old_df[col] = old_df[col].map(md5_dict).fillna(old_df[col])`.
- When iterating folders, list files inside the folder (use `os.listdir(folder_path)`) — do not iterate the folder name string.
- For seaborn barplots: to avoid the `palette` deprecation, assign the `y` variable to `hue` and set `legend=False`, or provide a color list via `color=`.

---

## Output files

- MGall_Manifest.csv — manifest built from raw reads folders.
- MGall_Manifest_with_MD5.csv — manifest with MD5 sums replaced.
- Manhattan plots — PNG files saved by `Manhattan_annotated.py`.
- Sensitivity plots — PNG files saved by `Paramenter_sensitivity_analysis_m.py`.

---

## Contact
For questions about the code or data layout, contact:
Dr. Alberto Torcello-Requena — alberto.torcello-requena@exeter.ac.uk

License: (none specified) — code provided "as is".

