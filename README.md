# MCR-ALS and MSPC monitoring of Septoria tritici blotch in wheat leaves: analysis code

MATLAB code for the analyses of:

> Latchoumane L., Ballini E., Ecarnot M., Roger J.-M., Bendoula R., Mas Garcia S. Process monitoring of Septoria tritici blotch in wheat leaves using hyperspectral imaging, Multivariate Curve Resolution-Alternating Least Squares and Multivariate Statistical Process Control. Submitted to *Computers and Electronics in Agriculture* (under revision).

Please cite the article once it is published; the reference will be added here.

The hyperspectral images, the leaf masks and the file `Data_description.csv` (genotype, treatment and confidence level of the infection status of each leaf) are available from Data INRAE (Latchoumane et al., <https://doi.org/10.57745/WVP0FJ>).

## Requirements

- Preprocessing (Step 1): MATLAB R2015a or later with the Image Processing Toolbox; run with MATLAB R2015a for the article.
- Analysis (Steps 2-3): MATLAB R2020a or later with the Statistics and Machine Learning Toolbox (`fitlme`, `fishertest`, `binocdf`, `finv`, `chi2inv`); run with MATLAB R2022a for the article.

No other toolbox is needed. MCR-ALS (non-negativity and unit-length spectra), SIMPLISMA and exact non-negative least squares are implemented in `functions/`.

If PLS_Toolbox is installed, its `dataset` class and `table2dataset` function conflict with `fitlme`. `run_analysis` removes PLS_Toolbox from the path for the current MATLAB session only (it is back after restarting MATLAB). If the mixed models still fail, restart MATLAB, run `restoredefaultpath; rehash toolboxcache`, go to the code folder and run `run_analysis('steps', {'main', 'sensitivity', 'figures'})`.

## Step 1. From the raw images to `alldataandmeta.mat` (`preprocessing/run_preprocessing.m`)

`run_preprocessing` creates `alldataandmeta.mat`, the input of the analysis (Step 2). Set the two default folders at the top of `run_preprocessing.m` (or pass them as options, e.g. `run_preprocessing('raw_dir', 'D:\images', 'work_dir', 'D:\work')`): `raw_dir`, the folder with the downloaded images, and `work_dir`, where a small results file is written for each image. Options `'pots'` and `'dates'` select a subset of the images.

Download from the database the VNIR archives of the dates 11/03 to 29/03/2021 (−1 to 17 dpi; one archive per date, with all pots) and extract each one into `raw_dir`. Each archive contains a folder named after its date (yyyymmdd) with one subfolder per image (`<pot>_<side>`, e.g. `33_D`, with `33_D.hyspex`, `33_D.hdr` and the leaf masks `33_D_leaf1/2/3.csv`), so the images end up in `raw_dir/<date>/<pot>_<side>/`. Only the 20 pots used in the article are processed (variable `pots_used`; set it to `[]` to process all pots).

For each image, the script:

1. reads the raw image (1024 × 1370 × 216, ENVI format) with `hsread`;
2. converts it to reflectance with the Spectralon reference (mean of the first 50 columns, line by line);
3. crops the inoculated area (781 × 1370 pixels);
4. reads the leaf masks of the database (`<pot>_<side>_leaf<k>.csv`: pixel coordinates x = column, y = row in the cropped image), which include the leaves delineated manually and keep the leaf numbering of the database (the same number identifies the same leaf at all dates; an image may have, e.g., leaves 2 and 3 only); images without leaf masks are skipped;
5. builds the cropped cube with the leaf spectra (zeros elsewhere), resizes it with `imresize` (factor 0.4, bicubic interpolation) to 313 × 548 × 216 and maps the leaf pixels to the resized image (coordinates × 0.4, rounded);
6. saves the spectra and pixel indices of each leaf in `work_dir/per_image/`.

The script can be stopped and restarted: images already in `work_dir/per_image/` are not processed again. At the end, all images are organised into `alldataandmeta.mat` (code folder). The images and leaves used in the article are listed in `metadata/leaves_by_date.csv`.

Requirements: 64-bit MATLAB R2015a or later with the Image Processing Toolbox. The images are read with `hsread` (X. Hadoux, 2012; included with the agreement of R. Bendoula, with `hsreadheader`, in `preprocessing/hsi_toolbox/`), which returns the cube as samples × lines × bands. Any other ENVI reader can be used (e.g. `hypercube` in MATLAB or the `spectral` package in Python), taking care of the order of the two spatial dimensions. About 6-8 GB of free memory are needed to process one image (the whole image is held twice in memory).

## Step 2. Input files of the analysis

The code folder must contain `alldataandmeta.mat` and `Data_description.csv`; the model profiles are part of the code:

| File | Content |
|---|---|
| `alldataandmeta.mat` | Leaf pixel spectra (images resized by a factor of 0.4) and image metadata, created in Step 1. |
| `Data_description.csv` | File of the database (download it from <https://doi.org/10.57745/WVP0FJ> and place it in the code folder) giving, for each leaf, genotype, treatment and the confidence level of its infection status (visual evaluation; leaves with a low confidence level are not used). |
| `functions/full_model_profiles.m` | Spectral profiles (`soptTots`, 4 × 106) of the multiset MCR-ALS model fitted on all images with the MCR-ALS graphical interface (Jaumot et al., 2005, Chemom. Intell. Lab. Syst. 76, 101-110; <http://www.mcrals.info/>), with SIMPLISMA initial estimates, non-negativity of contributions and spectra, and spectra normalised to unit length, written in the code with full double precision; and the wavelengths (nm) of the 216 VNIR bands. The profiles are the initial estimates of the refits in each fold, so the results of the article are reproduced exactly with them. |

The spectra in `alldataandmeta.mat` have the 216 bands of the VNIR images; bands 25-130 (476-760 nm) are used.

## Step 3. How to run the analysis

Go to the code folder and run:

```matlab
run_analysis
```

All results are written to `results/`. The script runs, in this order:

| Step | What it does | Article |
|---|---|---|
| `data` | Organises the leaf spectra per pot × date (`revision_data.mat`) | – |
| (quality check) | Excludes the leaf × date observations whose mask lies mostly on the background (fewer than 60% of the pixels above 0.15 reflectance at 760 nm; list in `excluded_leaf_observations.csv`) | Section 2.3 |
| `check` | Pre-inoculation check of the healthy reference (below) | Section 2.6 |
| `main` | Design table; singular values, MCR-ALS with 3-5 components and stability over random initialisations; leave-one-pot-out validation (MCR-ALS, PCA reference model and control limit fitted without the left-out pot; sensitivity, specificity, false-positive rate and balanced accuracy with pot-bootstrap 95% intervals; leaf-level detection date, Fisher test with Holm correction); linear mixed models (genotype × treatment × dpi; random pot and leaf) and divergence dates; percentage of pixels assigned to the infection-associated class; references built from mock-inoculated leaves of the same date or of the two preceding dates; metadata files | Tables 1-3, Tables S1, S2, S4, S5, Figs. S1-S3, Section 3.3 |
| `sensitivity` | Leave-one-pot-out validation with three PCs (≥ 95% variance) monitored with T² and SPE, with SIMPLISMA initial estimates computed in each fold, and with a non-negative projection; summary of all sensitivity analyses (`TableS3_S4_sensitivity_summary.xlsx`) | Tables S3, S4 |
| `figures` | Figures 2-6 (from the saved results) | Figs. 2-6 |

One or more steps can be run alone, e.g. `run_analysis('steps', {'figures'})` to redraw the figures; each step reads the results saved by the previous ones. Other options (output folder, number of bootstrap resamples, etc.) are listed at the top of `run_analysis.m`. Helper functions (MCR-ALS, SIMPLISMA, MSPC, metrics, maps) are in `functions/`.

Figures are saved as 300-dpi PNG files in Times New Roman. Files named `FigS_metrics_…` are check figures and are not in the article.

## Pots

The analysis uses the 20 pots of the original selection of the experiment. They are pots 25-48, in which both sides of the pot were imaged at least until 13 dpi (in pots 1-24, the leaves of one or both sides were sampled for RNA analyses at 4, 7 or 10 dpi; Latchoumane et al., 2025, Data in Brief 59, 111404), except pots 30, 31, 39 and 40, whose infection status was uncertain (low confidence level in the database for pots 30, 31 and 40, probable annotation error for pot 39). Using only pots with both sides keeps the pot, the experimental unit, complete. The pots are set in `default_config` (`cfg.inoc_EPO67`, `cfg.mock_EPO67`, `cfg.inoc_EPO68`, `cfg.mock_EPO68`) and in `pots_used` in `run_preprocessing.m`:

| Genotype | Treatment | Pots |
|---|---|---|
| EPO67 | Inoculated | 32, 33, 43, 44, 45 |
| EPO67 | Mock-inoculated | 34, 35, 36, 46, 47, 48 |
| EPO68 | Inoculated | 25, 26, 27, 37, 38 |
| EPO68 | Mock-inoculated | 28, 29, 41, 42 |

Leaves whose infection status was recorded with a low confidence level in the database (file `Data_description.csv` of the database, column "Confidence level" = 2, visual evaluation of symptoms) are not used, because their status is uncertain: among the 20 pots this concerns the three leaves of pot 42, side G (mock-inoculated) and leaf 2 of pot 44, side D (inoculated). The excluded leaves are listed in `results/metadata/excluded_leaves_low_confidence.csv`.

Mock-inoculated pots whose leaves already differed from the healthy reference before inoculation are excluded (step `check`): a pot is flagged when the proportion of its observations at −1 and 0 dpi above the control limit is significantly higher than 5% (one-sided binomial test, p < 0.05). Because an anomalous pot widens the reference and can mask another one, the most significant pot is removed and the check is repeated until no pot is flagged. The results of each iteration are in `results/all_20_pots/` (first iteration, also the sensitivity analysis with all 20 pots) and `results/reference_check/`, and the excluded pots in `results/metadata/excluded_pots.csv`. In the article, pots 36 and 46 (EPO67, mock-inoculated) were excluded.

## Reproducibility

The random seed is 2026 (pixel subsampling and bootstrap). With the same inputs and seed, all tables and figures are reproduced exactly. The main results of the article are printed on screen and written to `results/` (`Table_metrics_full_initfull_k4.xlsx`, `Table_mixed_models_full_initfull_k4.xlsx`, `Table_pixels_full_initfull_k4.xlsx`, `TableS3_S4_sensitivity_summary.xlsx`).

The folder `metadata/` contains the files written by the analysis (`results/metadata/`): `leaves_by_date.csv` (one row per leaf × date observation analysed: pot, side, leaf, date, dpi, genotype, treatment, number of pixels, image, validation fold), `folds.csv` (validation fold of each pot), `excluded_pots.csv` (pots excluded by the pre-inoculation check) `excluded_leaves_low_confidence.csv` (leaves with a low confidence level of their infection status) and `excluded_leaf_observations.csv` (observations excluded by the quality check of the leaf masks).

## Names

All file, function, variable and column names are in English, except `soptTots` (variable name of the original analysis) and the fields of the image metadata. The input files and their variables are:

| Name | Meaning |
|---|---|
| `alldataandmeta.mat` | leaf pixel spectra and image metadata |
| `soptTots` | spectral profiles of the full multiset MCR-ALS model (`functions/full_model_profiles.m`) |
| `wavelengths` | wavelengths (nm) of the 216 bands (`functions/full_model_profiles.m`) |
| `metadataG`, `metadataD`, fields `fecha`, `tipo`, `planta`, `fila`, `columna`, `filename`, `variables_cargadas`, `faltantes` | date, side, pot, date index, file index, file name, variables found and variables missing for each image |
| `G` / `D` | left / right side of the pot (one image per side) |

## Licence and contact

MIT licence (see `LICENSE`). Contact: Silvia Mas Garcia, INRAE (corresponding author of the article).
