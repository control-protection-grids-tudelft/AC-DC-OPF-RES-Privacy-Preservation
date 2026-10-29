from pathlib import Path
import scipy.io as sio
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False
})


script_dir = Path(__file__).resolve().parent


res_list = ["0.1", "0.5"]

colors = ["0.10", "0.55"]
fill_colors = ["0.25", "0.80"]
markers = ["x", "o"]
linestyles = ["-", "-."]
marker_sizes = [5, 9]

fig, ax = plt.subplots(figsize=(7, 5))

# =========================================================
# Original feasible region: semicircle
# P: 0--60 MW, Q: -60--60 MVAr
# =========================================================
P_ori = np.linspace(0, 60, 300)
Q_ori_upper = np.sqrt(60**2 - P_ori**2)
Q_ori_lower = -Q_ori_upper

P_ori_env = np.concatenate([P_ori, P_ori[::-1]])
Q_ori_env = np.concatenate([Q_ori_upper, Q_ori_lower[::-1]])

ax.fill(
    P_ori_env,
    Q_ori_env,
    color="0.90",
    alpha=0.35,
    edgecolor="none",
    label=r'$\mathrm{By\ Capacity\ Constraint}$'
)

ax.plot(
    P_ori,
    Q_ori_upper,
    linestyle=":",
    color="0.35",
    linewidth=1.5
)

ax.plot(
    P_ori,
    Q_ori_lower,
    linestyle=":",
    color="0.35",
    linewidth=1.5
)



# =========================================================
# Resolution-based feasible regions
# =========================================================
for i, res in enumerate(res_list):

    P_data = sio.loadmat(script_dir / f"P_{res}.mat")
    Q_data = sio.loadmat(script_dir / f"Q_{res}.mat")

    P_env2 = np.squeeze(P_data["P_env2"])
    Q_env2 = np.squeeze(Q_data["Q_env2"])

    ax.fill(
        P_env2,
        Q_env2,
        color=fill_colors[i],
        alpha=0.22,
        edgecolor="none"
    )

    ax.plot(
        P_env2,
        Q_env2,
        marker=markers[i],
        linestyle=linestyles[i],
        color=colors[i],
        linewidth=2.0,
        markersize=marker_sizes[i],
        markerfacecolor='none' if markers[i] == 'o' else colors[i],
        markeredgewidth=1.6,

        label=(
            r'$\mathrm{By\ ECM}\ '
            r'(\Delta p_{r2v}^{\mathrm{res}} := 1\mathrm{MW})$'
            if res == "0.1"
            else
            r'$\mathrm{By\ ECM}\ '
            r'(\Delta p_{r2v}^{\mathrm{res}} := 5\mathrm{MW})$'
        )
    )

# =========================================================
# Adjust
# =========================================================
ax.set_xlabel(
    r'$\mathrm{RES\ Real\ Power\ Output}\ p_{r2v}^{\mathrm{res}}\ [\mathrm{MW}]$',
    fontsize=14
)

ax.set_ylabel(
    r'$\mathrm{RES\ Reactive\ Power\ Output}\ q_{r2v}^{\mathrm{res}}\ [\mathrm{MVAr}]$',
    fontsize=14
)


ax.xaxis.set_major_formatter(
    FuncFormatter(lambda x, pos: rf'${x:g}$')
)

ax.yaxis.set_major_formatter(
    FuncFormatter(lambda y, pos: rf'${y:g}$')
)

ax.tick_params(axis='both', labelsize=14)


ax.set_xlim(0, 60)
ax.set_ylim(-60, 60)


ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)


ax.grid(
    True,
    color="0.92",
    linewidth=0.7
)

ax.legend(
    frameon=False,
    fontsize=12,
    loc="best"
)

plt.tight_layout()
plt.show()