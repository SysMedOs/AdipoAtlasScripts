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

import pandas as pd
import plotly.graph_objects as go
import natsort

lipid_file = r"data/Cer_Sankey_Plot.xlsx"
cer_df = pd.read_excel(lipid_file, na_values=0)

source_lst = cer_df["From"].to_list()
target_lst = cer_df["To"].to_list()
val_lst = cer_df["Value"].to_list()
row_lst = list(zip(source_lst, target_lst, val_lst))

raw_node_lst = source_lst + target_lst

fa_lst = natsort.natsorted(
    list(set([i for i in cer_df["FA_X"].to_list() if isinstance(i, str)]))
)

spb_order = [
    "SPB 18:0,O2",
    "SPB 16:1;O2",
    "SPB 17:1;O2",
    "SPB 18:1;O2",
    "SPB 19:1;O2",
    "SPB 20:1;O2",
    "SPB 18:2;O2",
    "SPB 17:0;O",
    "SPB 18:0;O",
    "SPB 19:0;O",
    "SPB 18:1;O",
    "SPB 18:2;O",
    "SPB 18:0;O3",
    "SPB 20:0;O3",
]

class_order = [
    "Dihydro-Cer",
    "Cer",
    "Hex-(n)-Cer",
    "Sphingadienine-Cer",
    "Hex-(n)-Sphingadienine-Cer",
    "Dihydro-Deoxy-Cer",
    "Deoxy-Cer",
    "Phyto-Cer",
    "Hex-(n)-Phyto-Cer",
]

node_lst = spb_order + class_order + fa_lst

link_dct = {}

source_idx_lst = []
target_idx_lst = []
link_val_lst = []

f_lst = spb_order + class_order

for s in f_lst:
    for r in row_lst:
        r_s = r[0]
        if r_s == s:
            r_t = r[1]
            r_v = r[2]
            source_idx_lst.append(node_lst.index(r_s))
            target_idx_lst.append(node_lst.index(r_t))
            link_val_lst.append(r_v)

fig = go.Figure(
    data=[
        go.Sankey(
            arrangement="snap",
            node={
                "label": node_lst,
                "pad": 15,
            },
            link=dict(
                source=source_idx_lst,
                target=target_idx_lst,
                value=link_val_lst,
            ),
        )
    ]
)

fig.update_layout(title_text="Plot", font_size=10)

js2 = fig.to_json()
f_name = "Sankey_Cer"
with open(f"data/{f_name}.json", "w") as jsf:
    jsf.write(js2)

print("File saved as:", f"data/{f_name}.json")
print("Run step2 to generate images.")
print("Edit the generated json file and run step2 again to change colors.")
print("Finished.")

# examples of the color codes
# "node": {
#     "label": [
#         "SPB 18:0,O2",
#         "SPB 16:1;O2",
#         "SPB 17:1;O2",
#         "SPB 18:1;O2",
#         "SPB 19:1;O2",
#         "SPB 20:1;O2",
#         "SPB 18:2;O2",
#         "SPB 17:0;O",
#         "SPB 18:0;O",
#         "SPB 19:0;O",
#         "SPB 18:1;O",
#         "SPB 18:2;O",
#         "SPB 18:0;O3",
#         "SPB 20:0;O3",
#         "Dihydro-Cer",
#         "Cer",
#         "Hex-(n)-Cer",
#         "Sphingadienine-Cer",
#         "Hex-(n)-Sphingadienine-Cer",
#         "Dihydro-Deoxy-Cer",
#         "Deoxy-Cer",
#         "Phyto-Cer",
#         "Hex-(n)-Phyto-Cer",
#         "FA 12:0",
#         "FA 14:0",
#         "FA 15:0",
#         "FA 16:0",
#         "FA 16:0;O",
#         "FA 17:0",
#         "FA 18:0",
#         "FA 18:1",
#         "FA 19:0",
#         "FA 20:0",
#         "FA 21:0",
#         "FA 22:0",
#         "FA 22:0;O",
#         "FA 22:1",
#         "FA 23:0",
#         "FA 23:0;O",
#         "FA 23:1",
#         "FA 24:0",
#         "FA 24:0;O",
#         "FA 24:1",
#         "FA 24:1;O",
#         "FA 24:2",
#         "FA 25:0",
#         "FA 25:1",
#         "FA 26:0;O",
#         "FA 26:1"
#     ],
#     "color": [
#           "salmon",
#           "darksalmon",
#           "darksalmon",
#           "darksalmon",
#           "darksalmon",
#           "darksalmon",
#           "lightblue",
#           "khaki",
#           "khaki",
#           "khaki",
#           "teal",
#           "olive",
#           "olive",
#           "purple",
#           "rosybrown",
#           "lightcoral",
#           "indianred",
#           "lightskyblue",
#           "mediumturquoise",
#           "mediumorchid",
#           "plum",
#           "mediumpurple",
#           "mediumpurple",
#           "mediumseagreen",
#           "mediumseagreen",
#           "mediumseagreen",
#           "mediumseagreen",
#           "cornflowerblue",
#           "mediumseagreen",
#           "mediumseagreen",
#           "burlywood",
#           "mediumseagreen",
#           "mediumseagreen",
#           "mediumseagreen",
#           "mediumseagreen",
#           "cornflowerblue",
#           "burlywood",
#           "mediumseagreen",
#           "cornflowerblue",
#           "burlywood",
#           "mediumseagreen",
#           "cornflowerblue",
#           "burlywood",
#           "cornflowerblue",
#           "goldenrod",
#           "mediumseagreen",
#           "burlywood",
#           "cornflowerblue",
#           "burlywood"
#         ],
#     "pad": 15
# }
