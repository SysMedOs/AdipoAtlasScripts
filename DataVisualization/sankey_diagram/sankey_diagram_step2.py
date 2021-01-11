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

import json
import plotly.graph_objects as go


# f_name = "Sankey_Cer"
f_name = "Sankey_Cer_colors"  # if color is modified

with open(f"data/{f_name}.json", "r") as jsf2:
    fig2 = go.Figure(data=json.load(jsf2))
    fig2.show()
    fig2.write_image(f"img/{f_name}.png", width=1000, height=1600)
    fig2.write_image(f"img/{f_name}.svg", width=1000, height=1600)

print("To change order, use the interactive interface in the web browser that pop up automatically.")
print("To save the image after rearrangement, just use screenshot function.")
print("Finished.")
