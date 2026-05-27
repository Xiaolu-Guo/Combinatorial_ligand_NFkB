# Combinatorial_ligand_NFkB

Simulation and analysis code for the NF-κB signaling model in the context of
ligand mixtures (combinatorial / dual-ligand stimulation).

This repository accompanies the paper:

> * Macrophages maintain signaling fidelity in response to ligand mixtures
> * V1 (2025).
> https://www.biorxiv.org/content/10.1101/2025.03.27.645368v1.abstract 

The workflow has two stages: first **run the simulations** (the scripts in
`run_on_server/`), which write `.mat` result files; then **reproduce the
figures** by running the relevant blocks in `run_me_main_pub.m`, which load
those `.mat` files and draw the panels.

---

## 1. Requirements

- **MATLAB** R2021a or newer (the code uses `heatmap`, `parula`, and other
  recent graphics functions).
- MATLAB toolboxes: **Statistics and Machine Learning Toolbox** (clustering /
  distance functions used in the Wasserstein-distance and codon analyses).
- **Python 3.9+** with Jupyter, only for the machine-learning notebooks in
  `python/` (TensorFlow/Keras, scikit-learn, pandas, numpy, matplotlib). These
  are optional and independent of the MATLAB pipeline.

---

## 2. Repository layout

| Path | Contents |
|------|----------|
| `run_me_main_pub.m` | Master plotting script. Each block reproduces specific figure panels (see Section 5). |
| `run_on_server/` | Driver scripts that run the simulations and save `.mat` result files (see Section 4). |
| `src/` | Core model and analysis functions (ODE model, single-cell sampling, codon calculation, plotting). |
| `lib/` | Helper utilities (metrics, peak detection, feature lists, initialization, plotting helpers). |
| `bin/` | Additional plotting / distance utilities (`w_distance.m`, RMSD pair plots). |
| `python/` | Jupyter notebooks for the classifier / LSTM analyses and their CSV inputs (`python/data/`). |
| `All_ligand_codon_2023.mat` | Fitted single-cell parameter dataset used to seed the simulations (see Section 3). |
| `LICENSE` | MIT license. |

---

## 3. Required data file: `All_ligand_codon_2023.mat`

The simulations are seeded from a dataset of fitted single-cell parameters
distributed with the paper.

1. Download `All_ligand_codon_2023.mat` from https://datadryad.org/downloads/file_stream/4082734  , data files of the
   paper: https://www.nature.com/articles/s41467-025-60901-3. This is the data file for all single-cell parameters, to run simulation of NFkB signaling dynamics.
2. Place it in the **repository root** (the same folder as `run_me_main_pub.m`).

Because several functions load this file by name only
(e.g. `load('All_ligand_codon_2023.mat')` in
`src/Sim_sc_diff_dose_comb_ligand_2025_CpGpIC_SR.m`), it must be **on the MATLAB
path** — adding the repository root with `addpath(pwd)` (or running from the
repository root) is sufficient.

---

## 4. Running the simulations (`run_on_server/`)

The scripts in `run_on_server/` generate the simulation `.mat` files consumed by
`run_me_main_pub.m`. **Run them from the repository root**, for example:

```matlab
% from the repository root (folder containing run_me_main_pub.m)
addpath(pwd)              % so All_ligand_codon_2023.mat is found
run('run_on_server/run_pair1.m')
```

Each script begins with `addpath('./lib/')`, `addpath('./src/')`,
`addpath('./bin/')`, which resolve correctly when the current directory is the
repository root. Output is written to `../test_data/` (one level **above** the
repository); the scripts create this folder automatically if it does not exist.

### What each script produces

**Population-level, all conditions (31 combinatorial conditions):**

- `Sim_high_dose_all_comb.m` → `Sim2025_all_comb_SR.mat`
  (CpG–polyIC competing for endosomal transport, pIC–TLR3 degradation increased).
- `Sim_high_dose_all_comb_CpG_pIC_no_compete.m` → `Sim2025_all_comb_CpG_PIC_no_compete.mat`
  (same model but **without** CpG–pIC competition).

**Single-cell, dose-matched ligand pairs (synergy / antagonism):**

- `run_pair1.m` … `run_pair10.m` → `Sim1008_<ligand1>_<ligand2>_sc_synergy_different_dose.mat`,
  one file per ligand pair. The 10 pairs are indexed identically in every script:

  1. CpG + PolyIC  2. LPS + PolyIC  3. TNF + PolyIC  4. TNF + LPS  5. LPS + CpG
  6. CpG + Pam3CSK  7. Pam3CSK + PolyIC  8. LPS + Pam3CSK  9. TNF + CpG  10. TNF + Pam3CSK

- `run_pair_LPS_Pam_CD14.m` → `Sim1008_LPS_Pam3CSK_sc_antag_Diff_CD14.mat`
  (LPS+Pam3CSK with varied CD14-related parameters `params29`/`params30`).
- `run_pair_all_dose2_IKKhill_1.m` → `Sim1008_<L1>_<L2>_sc_dose2_IKK_hill1.mat`
  (all pairs, two-dose subset, IKK activation Hill coefficient = 1).
- `run_pair_all_dose2_IKKhill_2.m` → `Sim1008_<L1>_<L2>_sc_dose2_IKK_hill2.mat`
  (all pairs, two-dose subset, IKK activation Hill coefficient = 2, default).

These can be run independently and, on a cluster, in parallel (one job per
`run_pairN.m`).

### Making outputs available to the plotting script

`run_me_main_pub.m` loads simulation results and experimental data from the folder set in its
`data_save_file_path` variable (see Section 5) (The data can be downloaded:
DOI:10.17632/hcg272hzyd.1). After running the simulations,
copy (or move) the generated `.mat` files from `../test_data/` into that folder
so the plotting blocks can find them.

---

## 5. Reproducing the figures (`run_me_main_pub.m`)

`run_me_main_pub.m` is organized as a series of blocks, each guarded by
`if 0 ... end` and labeled with the figure(s) it produces. To run a block,
change its `if 0` to `if 1` and execute the script (or run the block
interactively in the MATLAB editor). Leaving blocks at `if 0` keeps them off so
you can reproduce one panel at a time.

### Before running, edit the two paths at the top of the file

```matlab
data_save_file_path = '<path to the folder containing the .mat result files>';
fig_save_path       = '../subfigures/';   % where output PDFs are written
```

Set `data_save_file_path` to the **repository root** (or a `data/` subfolder of
it) — i.e. wherever you placed `All_ligand_codon_2023.mat` and the simulation
`.mat` files produced in Section 4. The script also runs
`addpath('./lib/')`, `addpath('./src/')`, `addpath('./bin/')`, so launch it from
the repository root.

### Block-to-figure map

| Block label in `run_me_main_pub.m` | Figure(s) | Required `.mat` input |
|------------------------------------|-----------|-----------------------|
| Responder ratio (sim / exp) | Fig. 1, S1 | `Sim2025_all_comb_CpG_PIC_no_compete.mat` (and experimental data) |
| Trajectory heatmaps (sim / exp) | Fig. S1, 2A, 2B, S4 | `Sim2025_all_comb_SR.mat`, `Sim2025_all_comb_CpG_PIC_no_compete.mat`, `All_exp_sim_codon_dual_ligand_202511_nonminmaxscaled.mat` |
| `compare_sample_supriya_202605` | Fig. 2C, 2F, S3 | (loaded inside the function) |
| `bar_diff_exp_old_20251118` | Fig. 2G | — |
| `codon_distrib_comb` | Fig. S1, S4, 3A, 3B, 4E, 5E, EV1, EV2 | `Sim2025_all_comb_SR.mat` / `Sim2025_all_comb_CpG_PIC_no_compete.mat` |
| `plot_hist_Wdis_compete_vs_non_compete` | Fig. 5E | both all-comb `.mat` files |
| Save data to ML format (`codon_save_to_ML_format_*`) | Fig. 3 (classifier inputs) | `Exp_data_metrics_2024.mat`, `Sim2025_all_comb_SR.mat` |
| `single_cell_synergy` | Fig. 4, Fig. 5, supp. | `Sim1008_<pair>_sc_synergy_different_dose.mat` |
| Pairwise heatmaps | Fig. S8 | `Sim1008_<pair>_sc_signergy_different_dose.mat` (see note below) |
| `Fig_compare_synergy_hill1_hill2` | Fig. 4D | `Sim1008_<pair>_sc_dose2_IKK_hill1/2.mat` |
| `draw_heatmaps_hill12` / `draw_W_dist_main_seperate_cal_codon` / `plot_integral_hist_hill1_hill2` | Fig. S11A, 4E/4F, S11E | Hill-coefficient `.mat` files |


(The figure labels above are taken directly from the comments in
`run_me_main_pub.m`; refer to those comments for the authoritative mapping, as a
few panels are produced by helper functions that load their own inputs.)

---

## 6. Python notebooks (`python/`)

The notebooks reproduce the machine-learning analyses (Fig. 3 and related). They
read the CSV files in `python/data/` (codon and trajectory feature matrices
exported from the MATLAB pipeline) and are self-contained — open them in Jupyter
and run top to bottom. Pretrained LSTM weights are provided as
`lstm_model(*).keras`.

---

## 7. Quick start

```matlab
% 1. Download All_ligand_codon_2023.mat (see Section 3) into the repo root.
% 2. From the repository root in MATLAB:
addpath(pwd)

% 3. Run a simulation (writes to ../test_data/):
run('run_on_server/Sim_high_dose_all_comb.m')

% 4. Copy the resulting .mat from ../test_data/ into your data folder, then
%    set data_save_file_path at the top of run_me_main_pub.m to that folder.

% 5. Open run_me_main_pub.m, flip the desired block from `if 0` to `if 1`,
%    and run it to regenerate the figure PDFs in fig_save_path.
```

---

## License

MIT License. See [LICENSE](LICENSE).
