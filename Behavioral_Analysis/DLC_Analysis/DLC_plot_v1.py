"""
DLC plotting utilities.

This script loads DeepLabCut (DLC) summary CSVs and generates a variety of plots
(barplots, speed/duration comparisons, and statistical outputs).
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats

DURATION_FILENAME = "_9000_36000"

def processed_data(genotype, manipulation, rig, age, day, datapath, midbrain, duration):

    duration_filename = DURATION_FILENAME

    if genotype == 'MitoPark':
        if rig == 'SquareOF':
            folder = [filename for filename in os.listdir(datapath) if filename.startswith('Analysis_' + rig + '_' + genotype + '_' + age + duration_filename)][0]
            data_ctrl = pd.read_csv(os.path.join(datapath + folder + '/Analysis_' + rig + '_' + genotype + '_' + age + '_' + day + 'DLC_Control_summary.csv' ))
            data_test = pd.read_csv(os.path.join(datapath + folder + '/Analysis_' + rig + '_' + genotype + '_' + age + '_' + day + 'DLC_Test_summary.csv' ))
            data_ctrl = data_ctrl[data_ctrl['DLCscorer'].str.startswith('Sq')]
            data_test = data_test[data_test['DLCscorer'].str.startswith('Sq')]
        elif rig == 'Cylinder':
            folder = [filename for filename in os.listdir(datapath) if filename.startswith('Analysis_' + rig + '_' + genotype + '_' + age + '_0_27000')][0]
            data_ctrl = pd.read_csv(os.path.join(datapath + folder + '/Analysis_' + rig + '_' + genotype + '_' + age + '_' + day + 'DLC_Control_summary.csv' ))
            data_test = pd.read_csv(os.path.join(datapath + folder + '/Analysis_' + rig + '_' + genotype + '_' + age + '_' + day + 'DLC_Test_summary.csv' ))
            data_ctrl = data_ctrl[data_ctrl['DLCscorer'].str.startswith('Cy')]
            data_test = data_test[data_test['DLCscorer'].str.startswith('Cy')]
        data_ctrl['key'] = 'Ctrl'
        data_test['key'] = 'MitoPark'
        data_ctrl['age'] = age
        data_test['age'] = age
    else:
        if rig == 'SquareOF':
            folder = [filename for filename in os.listdir(datapath) if filename.startswith('Analysis_' + rig + '_' + genotype + '_' + manipulation + duration_filename)][0]

            if genotype == 'DAT':
                data_ctrl = pd.read_csv(os.path.join(datapath + folder + '/Analysis_' + rig + '_' + genotype + '_' + manipulation + '_' + midbrain + '_' + day + 'DLC_Control_summary.csv' ))
                data_test = pd.read_csv(os.path.join(datapath + folder + '/Analysis_' + rig + '_' + genotype + '_' + manipulation + '_' + midbrain + '_' + day + 'DLC_Test_summary.csv' ))
            else:
                data_ctrl = pd.read_csv(os.path.join(datapath + folder + '/Analysis_' + rig + '_' + genotype + '_' + manipulation + '_' + day + 'DLC_Control_summary.csv' ))
                data_test = pd.read_csv(os.path.join(datapath + folder + '/Analysis_' + rig + '_' + genotype + '_' + manipulation + '_' + day + 'DLC_Test_summary.csv' ))

            data_ctrl = data_ctrl[data_ctrl['DLCscorer'].str.startswith('Sq')]
            data_test = data_test[data_test['DLCscorer'].str.startswith('Sq')]

        elif rig == 'Cylinder':
            folder = [filename for filename in os.listdir(datapath) if filename.startswith('Analysis_' + rig + '_' + genotype + '_' + manipulation + '_0_27000')][0]

            if genotype == 'DAT':
                data_ctrl = pd.read_csv(os.path.join(datapath + folder + '/Analysis_' + rig + '_' + genotype + '_' + manipulation + '_' + midbrain + '_' + day + 'DLC_Control_summary.csv' ))
                data_test = pd.read_csv(os.path.join(datapath + folder + '/Analysis_' + rig + '_' + genotype + '_' + manipulation + '_' + midbrain + '_' + day + 'DLC_Test_summary.csv' ))
            else:
                data_ctrl = pd.read_csv(os.path.join(datapath + folder + '/Analysis_' + rig + '_' + genotype + '_' + manipulation + '_' + day + 'DLC_Control_summary.csv' ))
                data_test = pd.read_csv(os.path.join(datapath + folder + '/Analysis_' + rig + '_' + genotype + '_' + manipulation + '_' + day + 'DLC_Test_summary.csv' ))

            data_ctrl = data_ctrl[data_ctrl['DLCscorer'].str.startswith('Cy')]
            data_test = data_test[data_test['DLCscorer'].str.startswith('Cy')]
        data_ctrl['key'] = 'Ctrl'
        data_test['key'] = genotype + '-' + manipulation
    return data_ctrl,data_test

def DLC_histgram(
    genotype, manipulation, dataframe, save_dir,
    move_or_stop, speed_th, xlim_min, xlim_max,
    day, midbrain, duration
):
    """
    Generate histogram(s) of stop/move bout durations comparing Ctrl vs test group.

    Parameters
    ----------
    genotype : str
        'MitoPark' triggers 3 subplots (8/16/24wks); otherwise 1 subplot across ages.
    manipulation : str
        Used in test label: f"{genotype}-{manipulation}" (for non-MitoPark).
    dataframe : list-like
        Expects dataframe[0] to be a pandas DataFrame with columns:
        - 'key' in {'Ctrl','MitoPark', f'{genotype}-{manipulation}'}
        - 'age' in {'8wks','16wks','24wks'}
        - f'{move_or_stop}_frames_len_list' (list-like of bout durations in seconds)
        - f'{move_or_stop}_bout_length(s)' (per-bout mean/summary used for mean line)
    save_dir : str
    move_or_stop : str
        'move' or 'stop' (determines which columns to read).
    speed_th : float
        Lower bound for histogram bins.
    xlim_min, xlim_max : float
        X-axis limits in seconds (also used to form bins).
    day, midbrain, duration : str
        duration in {'0_5min','5_10min','10_15min','15_20min','20_25min','25_30min','5_20min',''}.
    """

    # ----- duration → filename suffix -----
    duration_filename = DURATION_FILENAME

    df0 = dataframe[0].copy()
    ages = ["8wks", "16wks", "24wks"]
    list_col = f"{move_or_stop}_frames_len_list"
    mean_col = f"{move_or_stop}_bout_length(s)"
    bins = np.arange(max(speed_th, xlim_min), xlim_max, 0.1)  # consistent 0.1s step
    colors = {"Ctrl": "black", "MitoPark": "red"}  # default mapping
    test_label = f"{genotype}-{manipulation}"

    def _explode_numeric(df, col):
        if df.empty or col not in df.columns:
            return pd.DataFrame(columns=[col])
        out = df.explode(col, ignore_index=True)
        out[col] = pd.to_numeric(out[col], errors="coerce")
        return out[out[col].notna()]

    def _n(df):
        return int(len(df))

    def _panel(ax, df_left, label_left, df_right, label_right, title_suffix):
        """Draw a single overlay panel and return means."""
        # explode
        long_left  = _explode_numeric(df_left,  list_col)
        long_right = _explode_numeric(df_right, list_col)

        # hist + KDE
        sns.histplot(ax=ax, data=long_left,  x=list_col, bins=bins, stat="density",
                     kde=True, alpha=0.30, line_kws={"lw":3}, color=colors.get(label_left, "black"), legend=False)
        sns.histplot(ax=ax, data=long_right, x=list_col, bins=bins, stat="density",
                     kde=True, alpha=0.30, line_kws={"lw":3}, color=colors.get(label_right, "red"),  legend=False)

        # means
        means = {}
        if mean_col in df_left.columns and not df_left.empty:
            m_left = float(df_left[mean_col].mean())
            ax.axvline(m_left, color=colors.get(label_left, "black"), linestyle="--", lw=1.5)
            means[label_left] = m_left
        if mean_col in df_right.columns and not df_right.empty:
            m_right = float(df_right[mean_col].mean())
            ax.axvline(m_right, color=colors.get(label_right, "red"), linestyle="--", lw=1.5)
            means[label_right] = m_right

        # axes/labels
        ax.set_xlim(xlim_min, xlim_max)
        ax.spines[["right","top"]].set_visible(False)
        ax.set_xlabel(f"{move_or_stop.capitalize()} duration (s)", size=24)
        ax.set_ylabel("Density", size=24)
        ax.set_title(f"{move_or_stop.capitalize()} Duration Distribution {title_suffix}", size=18)
        ax.tick_params(which="both", labelsize=18)

        # legends
        handles = [
            Patch(facecolor=colors.get(label_left, "black"), edgecolor="none", alpha=0.30,
                  label=f"{label_left} (n={_n(df_left)})"),
            Patch(facecolor=colors.get(label_right, "red"), edgecolor="none", alpha=0.30,
                  label=f"{label_right} (n={_n(df_right)})"),
        ]
        leg1 = ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=15)
        ax.add_artist(leg1)
        mean_handle = Line2D([0], [0], linestyle="--", lw=1.5, color="black", label="Mean")

        return means

    os.makedirs(save_dir, exist_ok=True)

    if genotype == "MitoPark":
        # 3 panels: Ctrl vs MitoPark for each age
        fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharex=False, sharey=False)
        all_means = {}

        for j, age in enumerate(ages):
            df_ctrl_age = df0[(df0["key"] == "Ctrl") & (df0["age"] == age)]
            df_mp_age   = df0[(df0["key"] == "MitoPark") & (df0["age"] == age)]
            all_means[age] = _panel(axs[j], df_ctrl_age, "Ctrl", df_mp_age, "MitoPark", age)

        fig.tight_layout()
        out_path = os.path.join(
            save_dir,
            f"{move_or_stop.capitalize()}_duration_CtrlvsMP_8wks_16wks_24wks_{day}_{midbrain}_{duration_filename}.pdf"
        )
        plt.savefig(out_path, dpi=300, bbox_inches="tight")

        return out_path, all_means

    else:
        # 1 panel: Ctrl vs <genotype-manipulation> across ages
        df_ctrl = df0[df0["key"] == "Ctrl"]
        df_test = df0[df0["key"] == test_label]
        # ensure a color for custom label
        colors[test_label] = colors.get(test_label, "red")

        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        means = _panel(ax, df_ctrl, "Ctrl", df_test, test_label, f"{genotype}-{manipulation}")

        fig.tight_layout()
        out_path = os.path.join(
            save_dir,
            f"{move_or_stop.capitalize()}_duration_Ctrlvs{test_label}_{day}_{midbrain}_{duration_filename}.pdf"
        )
        plt.savefig(out_path, dpi=300, bbox_inches="tight")

        return out_path, means


def DLC_histgram_MP_age(
    genotype, dataframe, save_dir,
    move_or_stop, speed_th,
    xlim_min, xlim_max, ylim_min, ylim_max,
    day, midbrain, duration,
    bins=None,          # optional: custom bin edges
    id_col=None         # optional: animal ID column for n=unique animals (fallback: row count)
):
    """
    Plot duration histograms for Ctrl vs MitoPark across ages (8/16/24wks) in two subplots.
    Uses columns:
      - list column: f"{move_or_stop}_frames_len_list"      (list-like of bout durations in seconds)
      - mean column: f"{move_or_stop}_bout_length(s)"       (per-bout summary used for mean lines)

    Saves to:
      <save_dir>/<Move|Stop>_duration_Ctrl_or_MP_8wks_16wks_24wks_<day>_<duration_suffix>.pdf

    Returns
    -------
    out_path : str
        Saved PDF path.
    means : dict
        {'Ctrl': {'8wks': m, '16wks': m, '24wks': m},
         'MitoPark': {'8wks': m, '16wks': m, '24wks': m}}
    """

    # ----- duration → filename suffix -----
    duration_suffix = "_9000_36000"

    df0 = dataframe[0].copy()
    ages = ["8wks", "16wks", "24wks"]
    age_colors = {"8wks": "darkgreen", "16wks": "lightblue", "24wks": "red"}

    list_col = f"{move_or_stop}_frames_len_list"
    mean_col = f"{move_or_stop}_bout_length(s)"

    # ----- default bins (0.1 s resolution unless overridden) -----
    if bins is None:
        start = max(float(speed_th), float(xlim_min))
        bins = np.arange(start, float(xlim_max), 0.1)

    # ----- helpers -----
    def _explode_numeric(df, col):
        """Explode list column to long form and keep numeric values only."""
        if df is None or df.empty or col not in df.columns:
            return pd.DataFrame(columns=[col])
        out = df.explode(col, ignore_index=True)
        out[col] = pd.to_numeric(out[col], errors="coerce")
        return out[out[col].notna()]

    def _n_samples(df):
        """Legend n: unique animals if id_col is provided & present, else row count."""
        if df is None or df.empty:
            return 0
        if id_col and id_col in df.columns:
            return int(df[id_col].nunique(dropna=True))
        return int(len(df))

    def _plot_panel(ax, df_group, panel_label):
        """Overlay 3 ages on one axis for a given group (Ctrl or MitoPark)."""
        means = {}
        for age in ages:
            df_age = df_group[df_group["age"] == age]
            if df_age.empty:
                continue

            # histogram + KDE for each age
            long = _explode_numeric(df_age, list_col)
            sns.histplot(
                ax=ax, data=long, x=list_col, bins=bins,
                kde=True, stat="density", alpha=0.30,
                line_kws={"lw": 3}, color=age_colors[age], legend=False
            )

            # per-age mean (dashed line)
            if mean_col in df_age.columns and not df_age.empty:
                m = float(df_age[mean_col].mean())
                means[age] = m
                ax.axvline(m, color=age_colors[age], linestyle="--", lw=1.5)

        # axis cosmetics
        ax.set_xlim(xlim_min, xlim_max)
        ax.set_ylim(ylim_min, ylim_max)
        ax.spines[["right", "top"]].set_visible(False)
        ax.set_xlabel(f"{move_or_stop.capitalize()} duration (s)", size=24)
        ax.set_ylabel("Density", size=24)
        ax.set_title(f"{move_or_stop.capitalize()} Duration Distribution {panel_label}", size=18)
        ax.tick_params(which="both", labelsize=18)

        # legend: ages with (n=)
        handles = []
        for age in ages:
            df_age = df_group[df_group["age"] == age]
            n = _n_samples(df_age)
            handles.append(Patch(facecolor=age_colors[age], edgecolor="none", alpha=0.30,
                                 label=f"{age} (n={n})"))
        leg1 = ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=15)
        ax.add_artist(leg1)

        # mean style legend
        mean_handle = Line2D([0], [0], linestyle="--", lw=1.5, color="black", label="Mean (per age)")

        return means

    # ----- slice groups -----
    df_ctrl = df0[df0["key"] == "Ctrl"].copy()
    df_mp   = df0[df0["key"] == "MitoPark"].copy()

    # ----- figure -----
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    means_ctrl = _plot_panel(axs[0], df_ctrl, "Ctrl")
    means_mp   = _plot_panel(axs[1], df_mp, "MitoPark")

    # optional context
    fig.suptitle(
        f"{genotype} | {midbrain} | {move_or_stop.capitalize()} | {day} {duration}",
        y=1.02, fontsize=12
    )
    fig.tight_layout()

    # ----- save -----
    os.makedirs(save_dir, exist_ok=True)
    out_path = os.path.join(
        save_dir,
        f"{move_or_stop.capitalize()}_duration_Ctrl_or_MP_8wks_16wks_24wks_{day}_{duration_suffix}.pdf"
    )
    plt.savefig(out_path, dpi=300, bbox_inches="tight")

    return out_path, {"Ctrl": means_ctrl, "MitoPark": means_mp}


def DLC_histgram_speed(
    genotype, manipulation, dataframe, save_dir,
    mobile_or_immobile, threshold,
    xlim_min, xlim_max, ylim_min, ylim_max,
    day, midbrain, duration,
    bins_speed=None  # if None, will default to step=1 for 'mobile' and 0.1 for 'immobile'
):
    """
    Create 3-panel histograms (8wks, 16wks, 24wks) of DLC speeds comparing Ctrl vs MitoPark
    when genotype == 'MitoPark'; otherwise create a single-panel Ctrl vs <genotype-manipulation>
    across all ages. Overlays KDE, draws dashed mean lines, and adds legends with (n=#animals).

    Parameters
    ----------
    genotype : str
    manipulation : str
    dataframe : any
        Expects dataframe[0] to be a pandas.DataFrame with at least:
        columns: ['key' in {'Ctrl','MitoPark', '<genotype>-<manipulation>'}, 'age' in {'8wks','16wks','24wks'},
                  f"{mobile_or_immobile}_speed_list" (list-like), 'average_speed_spine2(cm/s)']
        Ideally also contains an animal ID column (e.g., 'animal_id', 'mouse_id', 'subject_id', ...).
    save_dir : str
    mobile_or_immobile : str
        'mobile' or 'immobile' (determines which *_speed_list to plot)
    threshold : float
        Included in the figure suptitle for reference.
    xlim_min, xlim_max, ylim_min, ylim_max : float
        Axis limits.
    day : str
    midbrain : str
    duration : str
        One of: '0_5min','5_10min','10_15min','15_20min','20_25min','25_30min','5_20min', or ''.
    bins_speed : array-like or None
        If None, defaults to np.arange(0, xlim_max, step=1) for 'mobile'
        and np.arange(0, xlim_max, step=0.1) for 'immobile'.

    Returns
    -------
    out_path : str
        Path to the saved PDF.
    means : dict
        For MitoPark branch: {'8wks': {'Ctrl': m, 'MitoPark': m}, '16wks': {...}, '24wks': {...}}
        For else branch:     {'Ctrl': m, '<genotype>-<manipulation>': m}
    """

    # ----- duration → duration_filename mapping -----
    duration_filename = DURATION_FILENAME

    # ----- defaults for bins -----
    if bins_speed is None:
        step = 1.0 if mobile_or_immobile.lower() == "mobile" else 0.1
        bins_speed = np.arange(0, float(xlim_max), step=step)

    # ----- helpers -----
    def _explode_numeric(df, col):
        """Explode list column to long form and keep numeric values only."""
        if df is None or df.empty or col not in df.columns:
            return pd.DataFrame(columns=[col])
        out = df.explode(col, ignore_index=True)
        out[col] = pd.to_numeric(out[col], errors="coerce")
        return out[out[col].notna()]

    def _infer_id_col(df):
        """Try to find an animal/subject ID column for n counting."""
        preferred = ["animal_id", "mouse_id", "subject_id", "subject", "animal", "mouse", "id"]
        for c in preferred:
            if c in df.columns:
                return c
        for c in df.columns:
            if isinstance(c, str) and c.lower().endswith("id"):
                return c
        for c in df.columns:
            if isinstance(c, str):
                lc = c.lower()
                if any(k in lc for k in ["animal", "mouse", "subject", "name"]):
                    return c
        return None

    def _n_animals(df):
        if df is None or df.empty:
            return 0
        col = _infer_id_col(df)
        return int(df[col].nunique(dropna=True)) if col else int(len(df))

    # ----- setup -----
    os.makedirs(save_dir, exist_ok=True)
    df0 = dataframe[0]  # master DataFrame
    speed_col = f"{mobile_or_immobile}_speed_list"

    if genotype == "MitoPark":
        # =======================
        # 3 panels: Ctrl vs MitoPark per age
        # =======================
        ages = ["8wks", "16wks", "24wks"]
        group_colors = {"Ctrl": "black", "MitoPark": "red"}

        # build lists
        df_ctrl_by_age = [df0[(df0["key"] == "Ctrl") & (df0["age"] == a)].copy() for a in ages]
        df_mp_by_age   = [df0[(df0["key"] == "MitoPark") & (df0["age"] == a)].copy() for a in ages]

        fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharex=False, sharey=False)
        means = {a: {} for a in ages}

        for j, age in enumerate(ages):
            ax = axs[j]
            df_ctrl = df_ctrl_by_age[j]
            df_mp   = df_mp_by_age[j]

            long_ctrl = _explode_numeric(df_ctrl, speed_col)
            long_mp   = _explode_numeric(df_mp, speed_col)

            # Plot (no seaborn legend spam)
            sns.histplot(
                ax=ax, data=long_ctrl, x=speed_col, bins=bins_speed,
                kde=True, stat="density", alpha=0.30, line_kws={"lw": 3},
                color=group_colors["Ctrl"], legend=False
            )
            sns.histplot(
                ax=ax, data=long_mp, x=speed_col, bins=bins_speed,
                kde=True, stat="density", alpha=0.30, line_kws={"lw": 3},
                color=group_colors["MitoPark"], legend=False
            )

            # Means from 'average_speed_spine2(cm/s)'
            if not df_ctrl.empty and "average_speed_spine2(cm/s)" in df_ctrl.columns:
                m_ctrl = float(df_ctrl["average_speed_spine2(cm/s)"].mean())
                ax.axvline(m_ctrl, color=group_colors["Ctrl"], linestyle="--", lw=1.5)
                means[age]["Ctrl"] = m_ctrl
            if not df_mp.empty and "average_speed_spine2(cm/s)" in df_mp.columns:
                m_mp = float(df_mp["average_speed_spine2(cm/s)"].mean())
                ax.axvline(m_mp, color=group_colors["MitoPark"], linestyle="--", lw=1.5)
                means[age]["MitoPark"] = m_mp

            # Cosmetics
            ax.set_xlim(xlim_min, xlim_max)
            ax.set_ylim(ylim_min, ylim_max)
            ax.spines[["right", "top"]].set_visible(False)
            ax.set_xlabel(f"{mobile_or_immobile.capitalize()} speed (cm/s)", size=24)
            if j == 0:
                ax.set_ylabel("Density", size=24)
            ax.set_title(f"{mobile_or_immobile.capitalize()} Speed Distribution {age}", size=18)
            ax.tick_params(which="both", labelsize=18)

            # Legends: group counts + mean style
            n_ctrl = _n_animals(df_ctrl)
            n_mp   = _n_animals(df_mp)
            group_handles = [
                Patch(facecolor=group_colors["Ctrl"], edgecolor="none", alpha=0.30, label=f"Ctrl (n={n_ctrl})"),
                Patch(facecolor=group_colors["MitoPark"], edgecolor="none", alpha=0.30, label=f"MitoPark (n={n_mp})"),
            ]
            leg1 = ax.legend(handles=group_handles, loc="upper right", frameon=False, fontsize=15)
            ax.add_artist(leg1)
            mean_handle = Line2D([0], [0], linestyle="--", lw=1.5, color="black", label="Mean")

        # Suptitle/context
        fig.suptitle(
            f"{genotype} | {manipulation} | {midbrain} | {mobile_or_immobile.capitalize()} "
            f"(th={threshold}) | {day} {duration}",
            y=1.02, fontsize=12
        )
        fig.tight_layout()

        out_path = os.path.join(
            save_dir,
            f"{mobile_or_immobile.capitalize()}_speed_CtrlvsMP_8wks_16wks_24wks_{day}_{midbrain}_{duration_filename}.pdf"
        )
        plt.savefig(out_path, dpi=300, bbox_inches="tight")

        return out_path, means

    else:
        # =======================
        # 1 panel: Ctrl vs <genotype-manipulation> across all ages
        # =======================
        test_label = f"{genotype}-{manipulation}"
        group_colors = {"Ctrl": "black", test_label: "red"}

        df_ctrl = df0[df0["key"] == "Ctrl"].copy()
        df_test = df0[df0["key"] == test_label].copy()

        long_ctrl = _explode_numeric(df_ctrl, speed_col)
        long_test = _explode_numeric(df_test, speed_col)

        fig, ax = plt.subplots(1, 1, figsize=(6, 5))

        sns.histplot(
            ax=ax, data=long_ctrl, x=speed_col, bins=bins_speed,
            kde=True, stat="density", alpha=0.30, line_kws={"lw": 3},
            color=group_colors["Ctrl"], legend=False
        )
        sns.histplot(
            ax=ax, data=long_test, x=speed_col, bins=bins_speed,
            kde=True, stat="density", alpha=0.30, line_kws={"lw": 3},
            color=group_colors[test_label], legend=False
        )

        means = {}
        if not df_ctrl.empty and "average_speed_spine2(cm/s)" in df_ctrl.columns:
            m_ctrl = float(df_ctrl["average_speed_spine2(cm/s)"].mean())
            ax.axvline(m_ctrl, color=group_colors["Ctrl"], linestyle="--", lw=1.5)
            means["Ctrl"] = m_ctrl
        if not df_test.empty and "average_speed_spine2(cm/s)" in df_test.columns:
            m_test = float(df_test["average_speed_spine2(cm/s)"].mean())
            ax.axvline(m_test, color=group_colors[test_label], linestyle="--", lw=1.5)
            means[test_label] = m_test

        ax.set_xlim(xlim_min, xlim_max)
        ax.set_ylim(ylim_min, ylim_max)
        ax.spines[["right", "top"]].set_visible(False)
        ax.set_xlabel(f"{mobile_or_immobile.capitalize()} speed (cm/s)", size=24)
        ax.set_ylabel("Density", size=24)
        ax.set_title(f"{mobile_or_immobile.capitalize()} Speed Distribution {test_label}", size=18)
        ax.tick_params(which="both", labelsize=18)

        n_ctrl = _n_animals(df_ctrl)
        n_test = _n_animals(df_test)
        group_handles = [
            Patch(facecolor=group_colors["Ctrl"], edgecolor="none", alpha=0.30, label=f"Ctrl (n={n_ctrl})"),
            Patch(facecolor=group_colors[test_label], edgecolor="none", alpha=0.30, label=f"{test_label} (n={n_test})"),
        ]
        leg1 = ax.legend(handles=group_handles, loc="upper right", frameon=False, fontsize=15)
        ax.add_artist(leg1)
        mean_handle = Line2D([0], [0], linestyle="--", lw=1.5, color="black", label="Mean")

        fig.suptitle(
            f"{genotype} | {manipulation} | {midbrain} | {mobile_or_immobile.capitalize()} "
            f"(th={threshold}) | {day} {duration}",
            y=1.02, fontsize=12
        )
        fig.tight_layout()

        out_path = os.path.join(
            save_dir,
            f"{mobile_or_immobile.capitalize()}_speed_Ctrlvs{test_label}_{day}_{midbrain}_{duration_filename}.pdf"
        )
        plt.savefig(out_path, dpi=300, bbox_inches="tight")

        return out_path, means

def DLC_histgram_MP_age_speed(
    genotype, dataframe, save_dir,
    mobile_or_immobile, threshold,
    xlim_min, xlim_max, ylim_min, ylim_max,
    day, midbrain, duration,
    *,                       # <-- everything after this is keyword-only (prevents the TypeError)
    bins_speed=None,
    id_col=None             # optional: set to your animal ID column to make n = unique animals
):
    """
    Plot histograms of DLC speeds (Ctrl vs MitoPark) across ages (8/16/24 wks).
    - One figure with 2 subplots: Ctrl (left) and MitoPark (right)
    - Each subplot overlays the 3 ages with distinct colors and dashed mean lines
    - Legends show age with (n=...) per panel, plus a mini legend for the mean line
    - Saves to: <save_dir>/<Mobile|Immobile>_speed_Ctrl_or_MP_8wks_16wks_24wks_<day>_<duration_suffix>.pdf

    Parameters
    ----------
    genotype : str
        Kept for consistency; not used in filtering (we use 'key' == 'Ctrl' / 'MitoPark').
    dataframe : list-like
        dataframe[0] must be a pandas DataFrame with columns:
            'key' in {'Ctrl', 'MitoPark'}
            'age' in {'8wks', '16wks', '24wks'}
            f'{mobile_or_immobile}_speed_list'  (list-like of speeds)
            'average_speed_spine2(cm/s)'        (for mean lines)
            (optionally an animal ID column if you want n = unique animals)
    save_dir : str
    mobile_or_immobile : str
        'mobile' or 'immobile'
    threshold : float
        Included for compatibility / title context (not used directly in plotting).
    xlim_min, xlim_max, ylim_min, ylim_max : float
        Axis limits.
    day, midbrain, duration : str
        duration in {'0_5min','5_10min','10_15min','15_20min','20_25min','25_30min','5_20min',''}.
    bins_speed : array-like, keyword-only
        If None, defaults to np.arange(0, xlim_max, step=1 for 'mobile' else 0.1 for 'immobile').
    id_col : str or None, keyword-only
        If provided and present in the data, n = nunique(id_col); otherwise n = row count.

    Returns
    -------
    out_path : str
        Saved PDF path.
    means : dict
        {'Ctrl': {'8wks': mean, '16wks': mean, '24wks': mean},
         'MitoPark': {'8wks': mean, '16wks': mean, '24wks': mean}}
    """

    # ---------- duration → filename suffix ----------
    duration_suffix = "_9000_36000"

    # ---------- defaults for bins ----------
    if bins_speed is None:
        step = 1.0 if mobile_or_immobile.lower() == "mobile" else 0.1
        bins_speed = np.arange(0, float(xlim_max), step=step)

    # ---------- setup ----------
    os.makedirs(save_dir, exist_ok=True)
    df0 = dataframe[0].copy()
    ages = ["8wks", "16wks", "24wks"]
    age_palette = {"8wks": "darkgreen", "16wks": "lightblue", "24wks": "red"}
    speed_list_col = f"{mobile_or_immobile}_speed_list"

    # ---------- helpers ----------
    def _explode_numeric(df, col):
        if df.empty or col not in df.columns:
            return pd.DataFrame(columns=[col])
        out = df.explode(col, ignore_index=True)
        out[col] = pd.to_numeric(out[col], errors="coerce")
        return out[out[col].notna()]

    def _n_samples(df):
        """Count n for legend. If id_col is provided and present, use nunique(id_col); else row count."""
        if df.empty:
            return 0
        if id_col and id_col in df.columns:
            return int(df[id_col].nunique(dropna=True))
        return int(len(df))

    def _plot_panel(ax, group_key, panel_title):
        """
        Plot one panel (group_key in {'Ctrl','MitoPark'}) overlaying 3 ages.
        Returns means per age.
        """
        means = {}
        for age in ages:
            color = age_palette[age]
            df_age = df0[(df0["key"] == group_key) & (df0["age"] == age)].copy()
            if df_age.empty:
                continue

            # Histogram + KDE
            long = _explode_numeric(df_age, speed_list_col)
            sns.histplot(
                ax=ax, data=long, x=speed_list_col,
                bins=bins_speed, kde=True, stat="density",
                alpha=0.30, line_kws={"lw": 3}, color=color, legend=False
            )

            # Mean line (use column average; df_age is already filtered by this age)
            if "average_speed_spine2(cm/s)" in df_age.columns and not df_age.empty:
                m = float(df_age["average_speed_spine2(cm/s)"].mean())
                means[age] = m
                ax.axvline(m, color=color, linestyle="--", lw=1.5)

        # Axes formatting
        ax.set_xlim(xlim_min, xlim_max)
        ax.set_ylim(ylim_min, ylim_max)
        ax.spines[["right", "top"]].set_visible(False)
        ax.set_xlabel(f"{mobile_or_immobile.capitalize()} speed (cm/s)", size=24)
        ax.set_ylabel("Density", size=24)
        ax.set_title(panel_title, size=18)
        ax.tick_params(which="both", labelsize=18)

        # Legend: Ages with (n=)
        handles = []
        for age in ages:
            df_age = df0[(df0["key"] == group_key) & (df0["age"] == age)]
            n = _n_samples(df_age)
            handles.append(
                Patch(facecolor=age_palette[age], edgecolor="none", alpha=0.30, label=f"{age} (n={n})")
            )
        leg1 = ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=15)
        ax.add_artist(leg1)

        # Mean legend (single dashed line sample)
        mean_handle = Line2D([0], [0], linestyle="--", lw=1.5, color="black", label="Mean (per age)")
        # ax.legend(handles=[mean_handle], loc="upper right", frameon=False, fontsize=15, bbox_to_anchor=(1, 0.75))

        return means

    # ---------- figure ----------
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    means_ctrl = _plot_panel(axs[0], "Ctrl", f"{mobile_or_immobile.capitalize()} Speed Distribution Ctrl")
    means_mp   = _plot_panel(axs[1], "MitoPark", f"{mobile_or_immobile.capitalize()} Speed Distribution MitoPark")

    fig.tight_layout()

    # ---------- save ----------
    out_path = os.path.join(
        save_dir,
        f"{mobile_or_immobile.capitalize()}_speed_Ctrl_or_MP_8wks_16wks_24wks_{day}_{duration_suffix}.pdf"
    )
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    # plt.close(fig)

    means = {"Ctrl": means_ctrl, "MitoPark": means_mp}
    return out_path, means



#################
#################


# ----------------------------
# Shared helpers
# ----------------------------
def _duration_suffix(duration: str) -> str:
    mapping = {
        "5_20min": "_9000_36000",
    }
    return mapping.get(duration, "")

def _explode_numeric(df: pd.DataFrame, col: str) -> pd.DataFrame:
    if df is None or df.empty or col not in df.columns:
        return pd.DataFrame(columns=[col])
    out = df.explode(col, ignore_index=True)
    out[col] = pd.to_numeric(out[col], errors="coerce")
    return out[out[col].notna()]

def _infer_id_col(df: pd.DataFrame):
    preferred = ["animal_id","mouse_id","subject_id","subject","animal","mouse","id"]
    for c in preferred:
        if c in df.columns: return c
    for c in df.columns:
        if isinstance(c, str) and c.lower().endswith("id"): return c
    for c in df.columns:
        if isinstance(c, str) and any(k in c.lower() for k in ["animal","mouse","subject","name"]):
            return c
    return None

def _n_samples(df: pd.DataFrame, id_col: str|None) -> int:
    if df is None or df.empty:
        return 0
    if id_col and id_col in df.columns:
        return int(df[id_col].nunique(dropna=True))
    return int(len(df))

def _age_palette():
    return {"8wks":"darkgreen","16wks":"lightblue","24wks":"red"}

def _group_colors(test_label: str):
    return {"Ctrl":"black", test_label:"red"}

def _finalize_axes(ax, xlim_min, xlim_max, ylim_min, ylim_max, xlabel, title):
    ax.set_xlim(xlim_min, xlim_max)
    if ylim_min is not None and ylim_max is not None:
        ax.set_ylim(ylim_min, ylim_max)
    ax.spines[["right","top"]].set_visible(False)
    ax.set_xlabel(xlabel, size=24)
    ax.set_ylabel("Density", size=24)
    ax.set_title(title, size=18)
    ax.tick_params(which="both", labelsize=18)

def _legend_groups(ax, left_label, right_label, left_color, right_color, n_left, n_right):
    handles = [
        Patch(facecolor=left_color,  edgecolor="none", alpha=0.30, label=f"{left_label} (n={n_left})"),
        Patch(facecolor=right_color, edgecolor="none", alpha=0.30, label=f"{right_label} (n={n_right})"),
    ]
    leg1 = ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=15)
    ax.add_artist(leg1)
    mean_handle = Line2D([0],[0], linestyle="--", lw=1.5, color="black", label="Mean")
    # ax.legend(handles=[mean_handle], loc="upper right", frameon=False, fontsize=15, bbox_to_anchor=(1, 0.75))


# ----------------------------
# Core engine (used by both wrappers)
# ----------------------------
def _dlc_hist_core(
    *,  # keyword-only to keep calls clean
    df0: pd.DataFrame,
    save_dir: str,
    genotype: str,
    manipulation: str,
    measure: str,          # 'speed' or 'duration'
    state: str,            # speed: 'mobile'|'immobile'; duration: 'move'|'stop'
    xlim_min: float, xlim_max: float,
    ylim_min: float|None, ylim_max: float|None,
    day: str, midbrain: str, duration: str,
    layout: str,           # 'per_age' | 'by_group' | 'auto'
    bins: np.ndarray|None,
    id_col: str|None,
    threshold: float|None  # optional for title context
):
    os.makedirs(save_dir, exist_ok=True)
    dur_suf = _duration_suffix(duration)
    ages = ["8wks","16wks","24wks"]

    if measure == "speed":
        list_col = f"{state}_speed_list"
        mean_col = "average_speed_spine2(cm/s)"
        xlabel  = f"{state.capitalize()} speed (cm/s)"
        if bins is None:
            step = 1.0 if state.lower()=="mobile" else 0.1
            bins = np.arange(0, float(xlim_max), step)
        stem = f"{state.capitalize()}_speed"
    elif measure == "duration":
        list_col = f"{state}_frames_len_list"
        mean_col = f"{state}_bout_length(s)"
        xlabel  = f"{state.capitalize()} duration (s)"
        if bins is None:
            bins = np.arange(float(xlim_min), float(xlim_max), 0.1)
        stem = f"{state.capitalize()}_duration"
    else:
        raise ValueError("measure must be 'speed' or 'duration'")

    test_label = f"{genotype}-{manipulation}" if genotype != "MitoPark" else "MitoPark"
    gcolors = _group_colors(test_label)
    apal = _age_palette()

    # Small helper to draw one overlay panel for two groups
    def _panel_groups(ax, df_left, left_label, df_right, right_label, title_suffix):
        long_left  = _explode_numeric(df_left,  list_col)
        long_right = _explode_numeric(df_right, list_col)

        sns.histplot(ax=ax, data=long_left,  x=list_col, bins=bins,
                     kde=True, stat="density", alpha=0.30, line_kws={"lw":3},
                     color=gcolors.get(left_label, "black"), legend=False)
        sns.histplot(ax=ax, data=long_right, x=list_col, bins=bins,
                     kde=True, stat="density", alpha=0.30, line_kws={"lw":3},
                     color=gcolors.get(right_label, "red"), legend=False)

        means = {}
        if mean_col in df_left.columns and not df_left.empty:
            mL = float(df_left[mean_col].mean()); means[left_label]=mL
            ax.axvline(mL, color=gcolors.get(left_label,"black"), linestyle="--", lw=1.5)
        if mean_col in df_right.columns and not df_right.empty:
            mR = float(df_right[mean_col].mean()); means[right_label]=mR
            ax.axvline(mR, color=gcolors.get(right_label,"red"), linestyle="--", lw=1.5)

        _finalize_axes(ax, xlim_min, xlim_max, ylim_min, ylim_max, xlabel,
                       f"{xlabel.split()[0]} {xlabel.split()[1]} Distribution {title_suffix}")
        _legend_groups(ax, left_label, right_label, gcolors.get(left_label,"black"),
                       gcolors.get(right_label,"red"),
                       _n_samples(df_left, id_col), _n_samples(df_right, id_col))
        return means

    # Helper to draw one panel overlaying 3 ages for a single group
    def _panel_ages(ax, df_group, panel_label):
        means = {}
        for age in ages:
            df_age = df_group[df_group["age"]==age]
            if df_age.empty: continue
            long = _explode_numeric(df_age, list_col)
            sns.histplot(ax=ax, data=long, x=list_col, bins=bins,
                         kde=True, stat="density", alpha=0.30, line_kws={"lw":3},
                         color=apal[age], legend=False)
            if mean_col in df_age.columns and not df_age.empty:
                m = float(df_age[mean_col].mean()); means[age]=m
                ax.axvline(m, color=apal[age], linestyle="--", lw=1.5)

        _finalize_axes(ax, xlim_min, xlim_max, ylim_min, ylim_max, xlabel,
                       f"{xlabel.split()[0]} {xlabel.split()[1]} Distribution {panel_label}")
        # Age legend + mean ref
        handles = [Patch(facecolor=apal[a], edgecolor="none", alpha=0.30,
                         label=f"{a} (n={_n_samples(df_group[df_group['age']==a], id_col)})") for a in ages]
        leg1 = ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=15)
        ax.add_artist(leg1)
        mean_handle = Line2D([0],[0], linestyle="--", lw=1.5, color="black", label="Mean (per age)")

        return means

    # Figure routing
    if genotype == "MitoPark":
        if layout == "by_group":
            # 2 panels: Ctrl | MitoPark (each overlays ages)
            df_ctrl = df0[df0["key"]=="Ctrl"].copy()
            df_mp   = df0[df0["key"]=="MitoPark"].copy()
            fig, axs = plt.subplots(1,2, figsize=(12,5), sharey=False)
            means_ctrl = _panel_ages(axs[0], df_ctrl, "Ctrl")
            means_mp   = _panel_ages(axs[1], df_mp, "MitoPark")
            means = {"Ctrl":means_ctrl, "MitoPark":means_mp}
            suffix = f"{stem}_Ctrl_or_MP_8wks_16wks_24wks_{day}_{_duration_suffix(duration)}"
        else:
            # 'per_age' (default): 3 panels with Ctrl vs MP per age
            fig, axs = plt.subplots(1,3, figsize=(15,5), sharex=False, sharey=False)
            means = {}
            for j, age in enumerate(ages):
                df_ctrl_age = df0[(df0["key"]=="Ctrl") & (df0["age"]==age)]
                df_mp_age   = df0[(df0["key"]=="MitoPark") & (df0["age"]==age)]
                means[age] = _panel_groups(axs[j], df_ctrl_age, "Ctrl", df_mp_age, "MitoPark", age)
            suffix = f"{stem}_CtrlvsMP_8wks_16wks_24wks_{day}_{_duration_suffix(duration)}"
    else:
        # Non-MitoPark → 1 panel (Ctrl vs <genotype-manipulation>) across ages
        test_label = f"{genotype}-{manipulation}"
        df_ctrl = df0[df0["key"]=="Ctrl"].copy()
        df_test = df0[df0["key"]==test_label].copy()
        fig, ax = plt.subplots(1,1, figsize=(6,5))
        means = _panel_groups(ax, df_ctrl, "Ctrl", df_test, test_label, test_label)
        suffix = f"{stem}_Ctrlvs{test_label}_{day}_{_duration_suffix(duration)}"

    # Context + save
    fig.suptitle(
        f"{genotype} | {manipulation} | {midbrain} | {xlabel} (th={threshold}) | {day} {duration}",
        y=1.02, fontsize=11
    )
    fig.tight_layout()
    out_path = os.path.join(save_dir, f"{suffix}.pdf")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")

    return out_path, means


# ----------------------------
# Public wrappers
# ----------------------------
def DLC_hist_duration(
    genotype, manipulation, dataframe, save_dir,
    move_or_stop, xlim_min, xlim_max, ylim_min, ylim_max,
    day, midbrain, duration,
    *, bins=None, id_col=None, layout="per_age", threshold=None
):
    """
    Duration histograms using:
      list_col = f"{move_or_stop}_frames_len_list"
      mean_col = f"{move_or_stop}_bout_length(s)"
    layout: 'per_age' (3 panels Ctrl vs MP) or 'by_group' (2 panels Ctrl|MP overlay ages).
    """
    df0 = dataframe[0].copy()
    return _dlc_hist_core(
        df0=df0, save_dir=save_dir, genotype=genotype, manipulation=manipulation,
        measure="duration", state=move_or_stop,
        xlim_min=xlim_min, xlim_max=xlim_max,
        ylim_min=ylim_min, ylim_max=ylim_max,
        day=day, midbrain=midbrain, duration=duration,
        layout=layout, bins=bins, id_col=id_col, threshold=threshold
    )


def DLC_hist_speed(
    genotype, manipulation, dataframe, save_dir,
    mobile_or_immobile, xlim_min, xlim_max, ylim_min, ylim_max,
    day, midbrain, duration,
    *, bins=None, id_col=None, layout="per_age", threshold=None
):
    """
    Speed histograms using:
      list_col = f"{mobile_or_immobile}_speed_list"
      mean_col = "average_speed_spine2(cm/s)"
    layout: 'per_age' (3 panels Ctrl vs MP) or 'by_group' (2 panels Ctrl|MP overlay ages).
    """
    df0 = dataframe[0].copy()
    return _dlc_hist_core(
        df0=df0, save_dir=save_dir, genotype=genotype, manipulation=manipulation,
        measure="speed", state=mobile_or_immobile,
        xlim_min=xlim_min, xlim_max=xlim_max,
        ylim_min=ylim_min, ylim_max=ylim_max,
        day=day, midbrain=midbrain, duration=duration,
        layout=layout, bins=bins, id_col=id_col, threshold=threshold
    )

#################
#################

def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.001:
        return '***'
    elif pvalue <= 0.01:
        return '**'
    elif pvalue <= 0.05:
        return '*'
    return 'ns'


def DLC_barplots_w_statistics(genotype, manipulation, day, dataframe, save_dir, animal, midbrain, duration, stat, ax0_y, ax1_y, ax2_y, ax3_y, ax4_y, ax5_y):

    # These are only for SquareOF
    # NOTE: duration windows removed; always using the supported window.
    duration_filename = DURATION_FILENAME
    cols1 = ['immobile_times(s)','average_speed_spine2(cm/s)','no_of_stops','stop_bout_length(s)', 'move_bout_length(s)']
    cols2 = ['number_of_rearing']
    cols = cols1 + cols2

    if genotype == 'MitoPark':
        animal1 = [animal[0],animal[1],animal[2],animal[3],animal[4],animal[5]]
        animal2 = [animal[6],animal[7],animal[8],animal[9],animal[10],animal[11]]

        if stat == 'Mann-Whitney':
            csv_title = 'Mann-Whitney_test'
            # Mann-Whitney U rank test
            data_stats1 = []
            data_stats2 = []
            for j in range(len(cols1)):
                for i in range(0, len(animal1), 2):
                    data_stats1.append(stats.mannwhitneyu(animal1[i][cols1[j]].dropna(), animal1[i+1][cols1[j]].dropna()))
            for j in range(len(cols2)):
                for i in range(0, len(animal2), 2):
                    data_stats2.append(stats.mannwhitneyu(animal2[i][cols2[j]].dropna(), animal2[i+1][cols2[j]].dropna()))
        elif stat == 'Student':
            csv_title = 'Student_t-test'
            # T-test for the means of two independent samples of scores
            data_stats1 = []
            data_stats2 = []
            for j in range(len(cols1)):
                for i in range(0, len(animal1), 2):
                    data_stats1.append(stats.ttest_ind(animal1[i][cols1[j]].dropna(), animal1[i+1][cols1[j]].dropna()))
            for j in range(len(cols2)):
                for i in range(0, len(animal2), 2):
                    data_stats2.append(stats.ttest_ind(animal2[i][cols2[j]].dropna(), animal2[i+1][cols2[j]].dropna()))

        data_stats = data_stats1 + data_stats2

        # for setting statistical bar/star
        p1y1 = max(dataframe[0][dataframe[0]['age']=='8wks']['immobile_times(s)'])+ax0_y/15 # arbitual y position for statistic star for visualization; automatically assigned by ax0_y
        p1y2 = max(dataframe[0][dataframe[0]['age']=='16wks']['immobile_times(s)'])+ax0_y/15
        p1y3 = max(dataframe[0][dataframe[0]['age']=='24wks']['immobile_times(s)'])+ax0_y/15

        p2y1 = max(dataframe[0][dataframe[0]['age']=='8wks']['average_speed_spine2(cm/s)'])+ax1_y/15
        p2y2 = max(dataframe[0][dataframe[0]['age']=='16wks']['average_speed_spine2(cm/s)'])+ax1_y/15
        p2y3 = max(dataframe[0][dataframe[0]['age']=='24wks']['average_speed_spine2(cm/s)'])+ax1_y/15

        p3y1 = max(dataframe[0][dataframe[0]['age']=='8wks']['no_of_stops'])+ax2_y/15
        p3y2 = max(dataframe[0][dataframe[0]['age']=='16wks']['no_of_stops'])+ax2_y/15
        p3y3 = max(dataframe[0][dataframe[0]['age']=='24wks']['no_of_stops'])+ax2_y/15

        p4y1 = max(dataframe[0][dataframe[0]['age']=='8wks']['stop_bout_length(s)'])+ax3_y/15
        p4y2 = max(dataframe[0][dataframe[0]['age']=='16wks']['stop_bout_length(s)'])+ax3_y/15
        p4y3 = max(dataframe[0][dataframe[0]['age']=='24wks']['stop_bout_length(s)'])+ax3_y/15

        p5y1 = max(dataframe[0][dataframe[0]['age']=='8wks']['move_bout_length(s)'])+ax4_y/15
        p5y2 = max(dataframe[0][dataframe[0]['age']=='16wks']['move_bout_length(s)'])+ax4_y/15
        p5y3 = max(dataframe[0][dataframe[0]['age']=='24wks']['move_bout_length(s)'])+ax4_y/15

        p6y1 = max(dataframe[1][dataframe[1]['age']=='8wks']['number_of_rearing'])+ax5_y/15
        p6y2 = max(dataframe[1][dataframe[1]['age']=='16wks']['number_of_rearing'])+ax5_y/15
        p6y3 = max(dataframe[1][dataframe[1]['age']=='24wks']['number_of_rearing'])+ax5_y/15

        pbar_height = [p1y1,p1y2,p1y3,p2y1,p2y2,p2y3,p3y1,p3y2,p3y3,p4y1,p4y2,p4y3,p5y1,p5y2,p5y3,p6y1,p6y2,p6y3]
        pvalue_axis1 = [-0.25,-0.25,0.25,0.25]
        pvalue_axis2 = [0.75,0.75,1.25,1.25]
        pvalue_axis3 = [1.75,1.75,2.25,2.25]
        plength = [ax0_y*0.03,ax1_y*0.03,ax2_y*0.03,ax3_y*0.03,ax4_y*0.03,ax5_y*0.03]

        # Images
        fig, axs = plt.subplots(1,6,figsize=(25,5))
        color = ['black', 'red']

        g1 = sns.barplot(ax=axs[0], y='immobile_times(s)', x='age', data=dataframe[0], alpha=0.5, err_kws={'color': 'grey'}, palette=color, hue='key')
        g1 = sns.stripplot(ax=axs[0], y='immobile_times(s)', x='age', data=dataframe[0], size=5, alpha=1.0, dodge=True, jitter=True, palette=color, hue='key', legend=False)
        axs[0].set_ylim(0, ax0_y)
        axs[0].spines[['right', 'top']].set_visible(False)
        axs[0].set_ylabel('Immobile time(s)', fontsize=24)
        axs[0].tick_params(axis='both', which='major', labelsize=18)
        axs[0].set(xlabel=None)
        axs[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False, fontsize=18)
        axs[0].plot(pvalue_axis1,[pbar_height[0], pbar_height[0]+plength[0], pbar_height[0]+plength[0], pbar_height[0]], lw=1.5, color = '0.2')
        axs[0].text(x=0, y=pbar_height[0]+plength[0]*1.5, s=convert_pvalue_to_asterisks(data_stats[0][1]), ha='center', size=25, weight='bold',color='0.2')
        axs[0].plot(pvalue_axis2,[pbar_height[1], pbar_height[1]+plength[0], pbar_height[1]+plength[0], pbar_height[1]], lw=1.5, color = '0.2')
        axs[0].text(x=1, y=pbar_height[1]+plength[0]*1.5, s=convert_pvalue_to_asterisks(data_stats[1][1]), ha='center', size=25, weight='bold',color='0.2')
        axs[0].plot(pvalue_axis3,[pbar_height[2], pbar_height[2]+plength[0], pbar_height[2]+plength[0], pbar_height[2]], lw=1.5, color = '0.2')
        axs[0].text(x=2, y=pbar_height[2]+plength[0]*1.5, s=convert_pvalue_to_asterisks(data_stats[2][1]), ha='center', size=25, weight='bold',color='0.2')

        g2 = sns.barplot(ax=axs[1],y='average_speed_spine2(cm/s)', x='age', data=dataframe[0], alpha=0.5, err_kws={'color': 'grey'}, palette=color, hue='key')
        g2 = sns.stripplot(ax=axs[1], y='average_speed_spine2(cm/s)', x='age', data=dataframe[0], size=5, alpha=1.0, dodge=True, jitter=True, palette=color, hue='key', legend=False)
        axs[1].set_ylim(0, ax1_y)
        axs[1].spines[['right', 'top']].set_visible(False)
        axs[1].set_ylabel('Average speed(cm/s)', fontsize=24)
        axs[1].tick_params(axis='both', which='major', labelsize=18)
        axs[1].set(xlabel=None)
        axs[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False, fontsize=18)
        axs[1].plot(pvalue_axis1,[pbar_height[3], pbar_height[3]+plength[1], pbar_height[3]+plength[1], pbar_height[3]], lw=1.5, color = '0.2')
        axs[1].text(x=0, y=pbar_height[3]+plength[1]*1.5, s=convert_pvalue_to_asterisks(data_stats[3][1]), ha='center', size=25, weight='bold',color='0.2')
        axs[1].plot(pvalue_axis2,[pbar_height[4], pbar_height[4]+plength[1], pbar_height[4]+plength[1], pbar_height[4]], lw=1.5, color = '0.2')
        axs[1].text(x=1, y=pbar_height[4]+plength[1]*1.5, s=convert_pvalue_to_asterisks(data_stats[4][1]), ha='center', size=25, weight='bold',color='0.2')
        axs[1].plot(pvalue_axis3,[pbar_height[5], pbar_height[5]+plength[1], pbar_height[5]+plength[1], pbar_height[5]], lw=1.5, color = '0.2')
        axs[1].text(x=2, y=pbar_height[5]+plength[1]*1.5, s=convert_pvalue_to_asterisks(data_stats[5][1]), ha='center', size=25, weight='bold',color='0.2')

        g3 = sns.barplot(ax=axs[2],y='no_of_stops', x='age', data=dataframe[0], alpha=0.5, err_kws={'color': 'grey'}, palette=color, hue='key')
        g3 = sns.stripplot(ax=axs[2], y='no_of_stops', x='age', data=dataframe[0], size=5, alpha=1.0, dodge=True, jitter=True, palette=color, hue='key', legend=False)
        axs[2].set_ylim(0, ax2_y)
        axs[2].spines[['right', 'top']].set_visible(False)
        axs[2].set_ylabel('Number of stops', fontsize=24)
        axs[2].tick_params(axis='both', which='major', labelsize=18)
        axs[2].set(xlabel=None)
        axs[2].legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False, fontsize=18)
        axs[2].plot(pvalue_axis1,[pbar_height[6], pbar_height[6]+plength[2], pbar_height[6]+plength[2], pbar_height[6]], lw=1.5, color = '0.2')
        axs[2].text(x=0, y=pbar_height[6]+plength[2]*1.5, s=convert_pvalue_to_asterisks(data_stats[6][1]), ha='center', size=25, weight='bold',color='0.2')
        axs[2].plot(pvalue_axis2,[pbar_height[7], pbar_height[7]+plength[2], pbar_height[7]+plength[2], pbar_height[7]], lw=1.5, color = '0.2')
        axs[2].text(x=1, y=pbar_height[7]+plength[2]*1.5, s=convert_pvalue_to_asterisks(data_stats[7][1]), ha='center', size=25, weight='bold',color='0.2')
        axs[2].plot(pvalue_axis3,[pbar_height[8], pbar_height[8]+plength[2], pbar_height[8]+plength[2], pbar_height[8]], lw=1.5, color = '0.2')
        axs[2].text(x=2, y=pbar_height[8]+plength[2]*1.5, s=convert_pvalue_to_asterisks(data_stats[8][1]), ha='center', size=25, weight='bold',color='0.2')

        g4 = sns.barplot(ax=axs[3], y='stop_bout_length(s)', x='age', data=dataframe[0], alpha=0.5, err_kws={'color': 'grey'}, palette=color, hue='key')
        g4 = sns.stripplot(ax=axs[3], y='stop_bout_length(s)', x='age', data=dataframe[0], size=5, alpha=1.0, dodge=True, jitter=True, palette=color, hue='key', legend=False)
        axs[3].set_ylim(0, ax3_y)
        axs[3].spines[['right', 'top']].set_visible(False)
        axs[3].set_ylabel('stop bout length(s)', fontsize=24)
        axs[3].tick_params(axis='both', which='major', labelsize=18)
        axs[3].set(xlabel=None)
        axs[3].legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False, fontsize=18)
        axs[3].plot(pvalue_axis1,[pbar_height[9], pbar_height[9]+plength[3], pbar_height[9]+plength[3], pbar_height[9]], lw=1.5, color = '0.2')
        axs[3].text(x=0, y=pbar_height[9]+plength[3]*1.5, s=convert_pvalue_to_asterisks(data_stats[9][1]), ha='center', size=25, weight='bold',color='0.2')
        axs[3].plot(pvalue_axis2,[pbar_height[10], pbar_height[10]+plength[3], pbar_height[10]+plength[3], pbar_height[10]], lw=1.5, color = '0.2')
        axs[3].text(x=1, y=pbar_height[10]+plength[3]*1.5, s=convert_pvalue_to_asterisks(data_stats[10][1]), ha='center', size=25, weight='bold',color='0.2')
        axs[3].plot(pvalue_axis3,[pbar_height[11], pbar_height[11]+plength[3], pbar_height[11]+plength[3], pbar_height[11]], lw=1.5, color = '0.2')
        axs[3].text(x=2, y=pbar_height[11]+plength[3]*1.5, s=convert_pvalue_to_asterisks(data_stats[11][1]), ha='center', size=25, weight='bold',color='0.2')

        g5 = sns.barplot(ax=axs[4], y='move_bout_length(s)', x='age', data=dataframe[0], alpha=0.5, err_kws={'color': 'grey'}, palette=color, hue='key')
        g5 = sns.stripplot(ax=axs[4], y='move_bout_length(s)', x='age', data=dataframe[0], size=5, alpha=1.0, dodge=True, jitter=True, palette=color, hue='key', legend=False)
        axs[4].set_ylim(0, ax4_y)
        axs[4].spines[['right', 'top']].set_visible(False)
        axs[4].set_ylabel('Move bout length(s)', fontsize=24)
        axs[4].tick_params(axis='both', which='major', labelsize=18)
        axs[4].set(xlabel=None)
        axs[4].legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False, fontsize=18)
        axs[4].plot(pvalue_axis1,[pbar_height[12], pbar_height[12]+plength[4], pbar_height[12]+plength[4], pbar_height[12]], lw=1.5, color = '0.2')
        axs[4].text(x=0, y=pbar_height[12]+plength[4]*1.5, s=convert_pvalue_to_asterisks(data_stats[12][1]), ha='center', size=25, weight='bold',color='0.2')
        axs[4].plot(pvalue_axis2,[pbar_height[13], pbar_height[13]+plength[4], pbar_height[13]+plength[4], pbar_height[13]], lw=1.5, color = '0.2')
        axs[4].text(x=1, y=pbar_height[13]+plength[4]*1.5, s=convert_pvalue_to_asterisks(data_stats[13][1]), ha='center', size=25, weight='bold',color='0.2')
        axs[4].plot(pvalue_axis3,[pbar_height[14], pbar_height[14]+plength[4], pbar_height[14]+plength[4], pbar_height[14]], lw=1.5, color = '0.2')
        axs[4].text(x=2, y=pbar_height[14]+plength[4]*1.5, s=convert_pvalue_to_asterisks(data_stats[14][1]), ha='center', size=25, weight='bold',color='0.2')

        g6 = sns.barplot(ax=axs[5], y='number_of_rearing', x='age', data=dataframe[1], alpha=0.5, err_kws={'color': 'grey'}, palette=color, hue='key')
        g6 = sns.stripplot(ax=axs[5], y='number_of_rearing', x='age', data=dataframe[1], size=5, alpha=1.0, dodge=True, jitter=True, palette=color, hue='key', legend=False)
        axs[5].set_ylim(0, ax5_y)
        axs[5].spines[['right', 'top']].set_visible(False)
        axs[5].set_ylabel('Number of rears', fontsize=24)
        axs[5].tick_params(axis='both', which='major', labelsize=18)
        axs[5].set(xlabel=None)
        axs[5].legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False, fontsize=18)
        axs[5].plot(pvalue_axis1,[pbar_height[15], pbar_height[15]+plength[5], pbar_height[15]+plength[5], pbar_height[15]], lw=1.5, color = '0.2')
        axs[5].text(x=0, y=pbar_height[15]+plength[5]*1.5, s=convert_pvalue_to_asterisks(data_stats[15][1]), ha='center', size=25, weight='bold',color='0.2')
        axs[5].plot(pvalue_axis2,[pbar_height[16], pbar_height[16]+plength[5], pbar_height[16]+plength[5], pbar_height[16]], lw=1.5, color = '0.2')
        axs[5].text(x=1, y=pbar_height[16]+plength[5]*1.5, s=convert_pvalue_to_asterisks(data_stats[16][1]), ha='center', size=25, weight='bold',color='0.2')
        axs[5].plot(pvalue_axis3,[pbar_height[17], pbar_height[17]+plength[5], pbar_height[17]+plength[5], pbar_height[17]], lw=1.5, color = '0.2')
        axs[5].text(x=2, y=pbar_height[17]+plength[5]*1.5, s=convert_pvalue_to_asterisks(data_stats[17][1]), ha='center', size=25, weight='bold',color='0.2')

        sns.set(style='white')
        fig.tight_layout()
        plt.savefig(save_dir + '/' + '8wks_16wks_24wks_' + day + '_' + duration_filename + '_barplot_w_statistics' + '_' + csv_title + '.pdf', dpi=300)

        with open(save_dir + '/' + 'Statistics_' + csv_title + '_' + day + '_' + duration_filename + '_' + 'MitoPark_mice.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Animal', 'Rig', 'Analysis','Statistical_test', 'Pvalue_ctrl_test', 'n_Ctrl', 'n_Test'])
            writer.writerow(['MitoPark_8wks', 'SquareOF', 'immobile_times(s)', csv_title, data_stats[0][1],len(animal[0]), len(animal[1])])
            writer.writerow(['MitoPark_16wks', 'SquareOF', 'immobile_times(s)', csv_title, data_stats[1][1],len(animal[2]), len(animal[3])])
            writer.writerow(['MitoPark_24wks', 'SquareOF', 'immobile_times(s)', csv_title, data_stats[2][1],len(animal[4]), len(animal[5])])
            writer.writerow(['MitoPark_8wks', 'SquareOF', 'average_speed_spine2(cm/s)', csv_title, data_stats[3][1],len(animal[0]), len(animal[1])])
            writer.writerow(['MitoPark_16wks', 'SquareOF', 'average_speed_spine2(cm/s)', csv_title, data_stats[4][1],len(animal[2]), len(animal[3])])
            writer.writerow(['MitoPark_24wks', 'SquareOF', 'average_speed_spine2(cm/s)', csv_title, data_stats[5][1],len(animal[4]), len(animal[5])])
            writer.writerow(['MitoPark_8wks', 'SquareOF', 'no_of_stops', csv_title, data_stats[6][1],len(animal[0]), len(animal[1])])
            writer.writerow(['MitoPark_16wks', 'SquareOF', 'no_of_stops', csv_title, data_stats[7][1],len(animal[2]), len(animal[3])])
            writer.writerow(['MitoPark_24wks', 'SquareOF', 'no_of_stops', csv_title, data_stats[8][1],len(animal[4]), len(animal[5])])
            writer.writerow(['MitoPark_8wks', 'SquareOF', 'stop_bout_length(s)', csv_title, data_stats[9][1],len(animal[0]), len(animal[1])])
            writer.writerow(['MitoPark_16wks', 'SquareOF', 'stop_bout_length(s)', csv_title, data_stats[10][1],len(animal[2]), len(animal[3])])
            writer.writerow(['MitoPark_24wks', 'SquareOF', 'stop_bout_length(s)', csv_title, data_stats[11][1],len(animal[4]), len(animal[5])])
            writer.writerow(['MitoPark_8wks', 'SquareOF', 'move_bout_length(s)', csv_title, data_stats[12][1],len(animal[0]), len(animal[1])])
            writer.writerow(['MitoPark_16wks', 'SquareOF', 'move_bout_length(s)', csv_title, data_stats[13][1],len(animal[2]), len(animal[3])])
            writer.writerow(['MitoPark_24wks', 'SquareOF', 'move_bout_length(s)', csv_title, data_stats[14][1],len(animal[4]), len(animal[5])])
            writer.writerow(['MitoPark_8wks', 'Cylinder', 'number_of_rearing', csv_title, data_stats[15][1],len(animal[0]), len(animal[1])])
            writer.writerow(['MitoPark_16wks', 'Cylinder', 'number_of_rearing', csv_title, data_stats[16][1],len(animal[2]), len(animal[3])])
            writer.writerow(['MitoPark_24wks', 'Cylinder', 'number_of_rearing', csv_title, data_stats[17][1],len(animal[4]), len(animal[5])])

    else:
        animal1 = [animal[0],animal[1]]
        animal2 = [animal[2],animal[3]]

        if stat == 'Mann-Whitney':
            csv_title = 'Mann-Whitney_test'
            # Mann-Whitney U rank test
            data_stats1 = []
            data_stats2 = []
            for j in range(len(cols1)):
                for i in range(0, len(animal1), 2):
                    data_stats1.append(stats.mannwhitneyu(animal1[i][cols1[j]].dropna(), animal1[i+1][cols1[j]].dropna()))
            for j in range(len(cols2)):
                for i in range(0, len(animal2), 2):
                    data_stats2.append(stats.mannwhitneyu(animal2[i][cols2[j]].dropna(), animal2[i+1][cols2[j]].dropna()))

        elif stat == 'Student':
            csv_title = 'Student_t-test'
            # T-test for the means of two independent samples of scores
            data_stats1 = []
            data_stats2 = []
            for j in range(len(cols1)):
                for i in range(0, len(animal1), 2):
                    data_stats1.append(stats.ttest_ind(animal1[i][cols1[j]].dropna(), animal1[i+1][cols1[j]].dropna()))
            for j in range(len(cols2)):
                for i in range(0, len(animal2), 2):
                    data_stats2.append(stats.ttest_ind(animal2[i][cols2[j]].dropna(), animal2[i+1][cols2[j]].dropna()))

        data_stats = data_stats1 + data_stats2

        # for setting statistical bar/star
        p1y1 = max(dataframe[0]['immobile_times(s)'])+ax0_y/15 # arbitual y position for statistic star for visualization; automatically assigned by ax0_y
        p2y1 = max(dataframe[0]['average_speed_spine2(cm/s)'])+ax1_y/15
        p3y1 = max(dataframe[0]['no_of_stops'])+ax2_y/15
        p4y1 = max(dataframe[0]['stop_bout_length(s)'])+ax3_y/15
        p5y1 = max(dataframe[0]['move_bout_length(s)'])+ax4_y/15
        p6y1 = max(dataframe[1]['number_of_rearing'])+ax5_y/15

        pbar_height = [p1y1,p2y1,p3y1,p4y1,p5y1,p6y1]
        pvalue_axis1 = [0,0,1,1]
        plength = [ax0_y*0.03,ax1_y*0.03,ax2_y*0.03,ax3_y*0.03,ax4_y*0.03,ax5_y*0.03]

        # Images
        fig, axs = plt.subplots(1,6,figsize=(20,5))
        color = ['black', 'red']

        g1 = sns.barplot(ax=axs[0], y='immobile_times(s)', x='key', data=dataframe[0], alpha=0.5, dodge=False, err_kws={'color': 'grey'}, palette=color, hue='key')
        g1 = sns.stripplot(ax=axs[0], y='immobile_times(s)', x='key', data=dataframe[0], size=5, alpha=1.0, dodge=False, jitter=True, palette=color, hue='key', legend=False)
        axs[0].set_ylim(0, ax0_y)
        axs[0].spines[['right', 'top']].set_visible(False)
        axs[0].set_ylabel('Immobile time(s)', fontsize=24)
        axs[0].tick_params(axis='both', which='major', labelsize=18)
        axs[0].set(xlabel=None)
        axs[0].plot(pvalue_axis1,[pbar_height[0], pbar_height[0]+plength[0], pbar_height[0]+plength[0], pbar_height[0]], lw=1.5, color = '0.2')
        axs[0].text(x=0.5, y=pbar_height[0]+plength[0]*1.5, s=convert_pvalue_to_asterisks(data_stats[0][1]), ha='center', size=25, weight='bold',color='0.2')

        g2 = sns.barplot(ax=axs[1], y='average_speed_spine2(cm/s)', x='key', data=dataframe[0], alpha=0.5, dodge=False, err_kws={'color': 'grey'}, palette=color, hue='key')
        g2 = sns.stripplot(ax=axs[1], y='average_speed_spine2(cm/s)', x='key', data=dataframe[0], size=5, alpha=1.0, dodge=False, jitter=True, palette=color, hue='key', legend=False)
        axs[1].set_ylim(0, ax1_y)
        axs[1].spines[['right', 'top']].set_visible(False)
        axs[1].set_ylabel('Average speed(cm/s)', fontsize=24)
        axs[1].tick_params(axis='both', which='major', labelsize=18)
        axs[1].set(xlabel=None)
        axs[1].plot(pvalue_axis1,[pbar_height[1], pbar_height[1]+plength[1], pbar_height[1]+plength[1], pbar_height[1]], lw=1.5, color = '0.2')
        axs[1].text(x=0.5, y=pbar_height[1]+plength[1]*1.5, s=convert_pvalue_to_asterisks(data_stats[1][1]), ha='center', size=25, weight='bold',color='0.2')

        g3 = sns.barplot(ax=axs[2], y='no_of_stops', x='key', data=dataframe[0], alpha=0.5, dodge=False, err_kws={'color': 'grey'}, palette=color, hue='key')
        g3 = sns.stripplot(ax=axs[2], y='no_of_stops', x='key', data=dataframe[0], size=5, alpha=1.0, dodge=False, jitter=True, palette=color, hue='key', legend=False)
        axs[2].set_ylim(0, ax2_y)
        axs[2].spines[['right', 'top']].set_visible(False)
        axs[2].set_ylabel('Number of stops', fontsize=24)
        axs[2].tick_params(axis='both', which='major', labelsize=18)
        axs[2].set(xlabel=None)
        axs[2].plot(pvalue_axis1,[pbar_height[2], pbar_height[2]+plength[2], pbar_height[2]+plength[2], pbar_height[2]], lw=1.5, color = '0.2')
        axs[2].text(x=0.5, y=pbar_height[2]+plength[2]*1.5, s=convert_pvalue_to_asterisks(data_stats[2][1]), ha='center', size=25, weight='bold',color='0.2')

        g4 = sns.barplot(ax=axs[3], y='stop_bout_length(s)', x='key', data=dataframe[0], alpha=0.5, dodge=False, err_kws={'color': 'grey'}, palette=color, hue='key')
        g4 = sns.stripplot(ax=axs[3], y='stop_bout_length(s)', x='key', data=dataframe[0], size=5, alpha=1.0, dodge=False, jitter=True, palette=color, hue='key', legend=False)
        axs[3].set_ylim(0, ax3_y)
        axs[3].spines[['right', 'top']].set_visible(False)
        axs[3].set_ylabel('Stop bout length(s)', fontsize=24)
        axs[3].tick_params(axis='both', which='major', labelsize=18)
        axs[3].set(xlabel=None)
        axs[3].plot(pvalue_axis1,[pbar_height[3], pbar_height[3]+plength[3], pbar_height[3]+plength[3], pbar_height[3]], lw=1.5, color = '0.2')
        axs[3].text(x=0.5, y=pbar_height[3]+plength[3]*1.5, s=convert_pvalue_to_asterisks(data_stats[3][1]), ha='center', size=25, weight='bold',color='0.2')

        g5 = sns.barplot(ax=axs[4], y='move_bout_length(s)', x='key', data=dataframe[0], alpha=0.5, dodge=False, err_kws={'color': 'grey'}, palette=color, hue='key')
        g5 = sns.stripplot(ax=axs[4], y='move_bout_length(s)', x='key', data=dataframe[0], size=5, alpha=1.0, dodge=False, jitter=True, palette=color, hue='key', legend=False)
        axs[4].set_ylim(0, ax4_y)
        axs[4].spines[['right', 'top']].set_visible(False)
        axs[4].set_ylabel('Move bout length(s)', fontsize=24)
        axs[4].tick_params(axis='both', which='major', labelsize=18)
        axs[4].set(xlabel=None)
        axs[4].plot(pvalue_axis1,[pbar_height[4], pbar_height[4]+plength[4], pbar_height[4]+plength[4], pbar_height[4]], lw=1.5, color = '0.2')
        axs[4].text(x=0.5, y=pbar_height[4]+plength[4]*1.5, s=convert_pvalue_to_asterisks(data_stats[4][1]), ha='center', size=25, weight='bold',color='0.2')

        g6 = sns.barplot(ax=axs[5], y='number_of_rearing', x='key', data=dataframe[1], alpha=0.5, dodge=False, err_kws={'color': 'grey'}, palette=color, hue='key')
        g6 = sns.stripplot(ax=axs[5], y='number_of_rearing', x='key', data=dataframe[1], size=5, alpha=1.0, dodge=False, jitter=True, palette=color, hue='key', legend=False)
        axs[5].set_ylim(0, ax5_y)
        axs[5].spines[['right', 'top']].set_visible(False)
        axs[5].set_ylabel('Number of rears', fontsize=24)
        axs[5].tick_params(axis='both', which='major', labelsize=18)
        axs[5].set(xlabel=None)
        axs[5].plot(pvalue_axis1,[pbar_height[5], pbar_height[5]+plength[5], pbar_height[5]+plength[5], pbar_height[5]], lw=1.5, color = '0.2')
        axs[5].text(x=0.5, y=pbar_height[5]+plength[5]*1.5, s=convert_pvalue_to_asterisks(data_stats[5][1]), ha='center', size=25, weight='bold',color='0.2')

        sns.set(style='white')
        fig.tight_layout()
        plt.savefig(save_dir + '/' + genotype + '_' + manipulation + '_' + day + '_' + midbrain + duration_filename + '_barplot_w_statistics' + '_' + csv_title + '.pdf', dpi=300)

        with open(save_dir + '/' + 'Statistics_' + csv_title + '_' + genotype + '_' + manipulation + '_' + day + '_' + midbrain + duration_filename + '_mice.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Animal', 'Rig', 'Analysis','Statistical_test', 'Pvalue_ctrl_test', 'n_Ctrl', 'n_Test'])
            writer.writerow([genotype + '-' + manipulation, 'SquareOF', 'immobile_times(s)', csv_title, data_stats[0][1],len(animal[0]), len(animal[1])])
            writer.writerow([genotype + '-' + manipulation, 'SquareOF', 'average_speed_spine2(cm/s)', csv_title, data_stats[1][1],len(animal[0]), len(animal[1])])
            writer.writerow([genotype + '-' + manipulation, 'SquareOF', 'no_of_stops', csv_title, data_stats[2][1],len(animal[0]), len(animal[1])])
            writer.writerow([genotype + '-' + manipulation, 'SquareOF', 'stop_bout_length(s)', csv_title, data_stats[3][1],len(animal[0]), len(animal[1])])
            writer.writerow([genotype + '-' + manipulation, 'SquareOF', 'move_bout_length(s)', csv_title, data_stats[4][1],len(animal[0]), len(animal[1])])
            writer.writerow([genotype + '-' + manipulation, 'Cylinder', 'number_of_rearing', csv_title, data_stats[5][1],len(animal[2]), len(animal[3])])

