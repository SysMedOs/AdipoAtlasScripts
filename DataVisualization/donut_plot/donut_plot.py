# -*- coding: utf-8 -*-
# Copyright (C) 2016-2021  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
# LipidHunter is Dual-licensed
#     For academic and non-commercial use: `GPLv2 License` Please read more information by the following link:
#         [The GNU General Public License version 2] (https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
#     For commercial use:
#         please contact the SysMedOs_team by email.
# Please cite our publication in an appropriate form.
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import matplotlib.pyplot as plt
import pandas as pd

in_file = r"data/donut_plot_data.csv"
f_name = "donut_plot_data"
df = pd.read_csv(in_file)
# Make data: I have 3 groups and 7 subgroups
subgroup_names = df["CLASS"].values.tolist()
subgroup_size = df["COUNT"].values.tolist()
subgroup_dct = dict(zip(subgroup_names, subgroup_size))

groups = {
    "CAR": ["CAR"],
    "lysoPL": [
        "LPC",
        "LPE",
        "LPI",
    ],
    "PL": [
        "PC",
        "O-PC",
        "P-PC",
        "PE",
        "O-PE",
        "P-PE",
        "PS",
        "PG",
        "PI",
    ],
    "SP": [
        "SM",
        "dhCer",
        "Cer",
        "Hex(n)Cer",
        "sphingadienineCer",
        "Hex(n)sphingadienineCer",
        "dh-deoxyCer",
        "deoxyCer",
        "phytoCer",
        "Hex(n)phytoCer",
    ],
}

csv_group_names = ["SP", "PL", "lysoPL", "CAR"]
group_size = []
defined_colors = [plt.cm.RdPu, plt.cm.PuBu, plt.cm.YlGn, plt.cm.Greys, plt.cm.PuBuGn]
group_names = []
group_colors = []
subgroup_colors = []
subgroup_size = []
subgroup_names = []
counter = 0

for group in csv_group_names:
    tmp_subgroups = groups[group]
    tmp_group_df = df[df["CLASS"].isin(tmp_subgroups)]
    tmp_group_size = tmp_group_df["COUNT"].sum()
    group_names.append(f"{group}: {tmp_group_size}")
    group_size.append(tmp_group_size)
    group_colors.append(defined_colors[counter](0.95))
    sub_counter = 1
    for i, r in tmp_group_df.iterrows():
        subgroup_colors.append(defined_colors[counter](0.9 - 0.075 * sub_counter))
        subgroup_size.append(r["COUNT"])
        subgroup_names.append(f'{r["CLASS"]}: {r["COUNT"]}')
        sub_counter += 1
    counter += 1

# Create colors

# # First Ring (inside)
fig, ax = plt.subplots(figsize=(20, 20))
ax.axis("equal")
mypie, labels1 = ax.pie(
    group_size,
    radius=0.44,
    labels=group_names,
    colors=group_colors,
    labeldistance=0.71,
    textprops={"fontsize": 15, "fontweight": "bold"},
    startangle=296,
)
plt.setp(mypie, width=0.27, edgecolor="white", alpha=0.5)
# do the rotation of the labels
for ea, eb in zip(mypie, labels1):
    mean_angle = (ea.theta1 + ea.theta2) / 2.0  # get mean_angle of the wedge

    eb.set_rotation(mean_angle)  # rotate the label by (mean_angle + 270)
    # eb.set_va("center")
    eb.set_ha("center")

# Second Ring (outside)
mypie2, labels2 = ax.pie(
    subgroup_size,
    radius=0.78,
    labels=subgroup_names,
    labeldistance=0.785,
    colors=subgroup_colors,
    textprops={"fontsize": 14, "fontweight": "bold", "ha": "right"},
    startangle=296,
)
plt.setp(mypie2, width=0.34, edgecolor="white", alpha=0.4)

# do the rotation of the labels
for ea, eb in zip(mypie2, labels2):
    mean_angle = (ea.theta1 + ea.theta2) / 2.0  # get mean_angle of the wedge

    eb.set_rotation(mean_angle)  # rotate the label by (mean_angle + 270)
    eb.set_va("center")
    eb.set_ha("center")
    lb = str(eb.get_text())
    if "Cer" in lb and len(lb) > 7:
        eb.set_fontsize(9)
    # elif "LPC" in lb or "LPE" in lb or "LPI" in lb:
    #     eb.set_fontsize(10)

# TG Ring (outside)

count_tg = subgroup_dct.get("TG")
count_dg = subgroup_dct.get("DG")
count_ce = subgroup_dct.get("CE")

count_polar = 0
for k in subgroup_dct:
    if k not in ["TG", "DG", "CE"]:
        count_polar += subgroup_dct.get(k)

mypie3, labels3 = ax.pie(
    [count_tg, count_dg, count_ce, count_polar],
    radius=0.92,
    labels=[
        f"TG: {count_tg}",
        f"DG:  {count_dg}",
        f"CE:  {count_ce}",
        f"Polar Lipids: {count_polar}",
    ],
    colors=["orangered", "orange", "coral", "teal"],
    labeldistance=0.925,
    textprops={"fontsize": 18, "fontweight": "bold"},
)

for ea, eb in zip(mypie3, labels3):
    mean_angle = (ea.theta1 + ea.theta2) / 2.0  # get mean_angle of the wedge
    eb.set_rotation(mean_angle + 270)  # rotate the label by (mean_angle + 270)
    eb.set_va("center")
    eb.set_ha("center")

plt.setp(mypie3, width=0.12, edgecolor="white", alpha=0.45)
plt.margins(0, 0)
plt.savefig(f"img/{f_name}.svg", type="svg", dpi=600)
plt.savefig(f"img/{f_name}.png", type="png", dpi=600)


## generate again without labels
# # First Ring (inside)
fig, ax = plt.subplots(figsize=(20, 20))
ax.axis("equal")
mypie, labels1 = ax.pie(
    group_size,
    radius=0.44,
    colors=group_colors,
    startangle=296,
)
plt.setp(mypie, width=0.27, edgecolor="white", alpha=0.5)
# do the rotation of the labels
for ea, eb in zip(mypie, labels1):
    mean_angle = (ea.theta1 + ea.theta2) / 2.0  # get mean_angle of the wedge

    eb.set_rotation(mean_angle)  # rotate the label by (mean_angle + 270)
    # eb.set_va("center")
    eb.set_ha("center")

# Second Ring (outside)
mypie2, labels2 = ax.pie(
    subgroup_size,
    radius=0.78,
    colors=subgroup_colors,
    startangle=296,
)
plt.setp(mypie2, width=0.34, edgecolor="white", alpha=0.4)

# TG Ring (outside)

count_tg = subgroup_dct.get("TG")
count_dg = subgroup_dct.get("DG")
count_ce = subgroup_dct.get("CE")

count_polar = 0
for k in subgroup_dct:
    if k not in ["TG", "DG", "CE"]:
        count_polar += subgroup_dct.get(k)

mypie3, labels3 = ax.pie(
    [count_tg, count_dg, count_ce, count_polar],
    radius=0.92,
    colors=["orangered", "orange", "coral", "teal"],
)

plt.setp(mypie3, width=0.12, edgecolor="white", alpha=0.45)
plt.margins(0, 0)
plt.savefig(f"img/{f_name}_no-label.svg", type="svg", dpi=300)
plt.savefig(f"img/{f_name}_no-label.png", type="png", dpi=300)

print("Finished.")
