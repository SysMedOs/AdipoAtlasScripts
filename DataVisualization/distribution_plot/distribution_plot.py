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
import matplotlib.patches as patches
import matplotlib.pyplot as plt

# fmol/ug

plt.rcParams["font.family"] = "Times New Roman"


class DistributionPlot(object):
    def __init__(self, in_file):
        in_df = pd.read_excel(in_file).fillna(0)
        q_df = in_df.groupby(["Class", "Quant_Bulk"])
        groups = q_df.groups
        out_df = pd.DataFrame()
        for g in groups:
            df = q_df.get_group(g)
            if g[-1] > 0:
                n_df = df.copy()
                n_df["C"] = df["C"].div(df.shape[0])
                n_df["Class"] = df["Class"]
                n_df["Lipid"] = df["Lipid"]
                n_df["Quant_Bulk"] = df["Quant_Bulk"]
                print(df.head())
                out_df = out_df.append(n_df)
            else:
                out_df = out_df.append(df)
        out_df = out_df.sort_index()
        out_df.to_excel(f"{in_file}_AVG.xlsx")
        self.df = out_df

        self.lipid_class_lst = self.df["Class"].unique().tolist()

        self.colors = {
            "Acylcarnitines": "forestgreen",
            "Diacylglycerols": "firebrick",
            "Cholesteryl Esters": "orangered",
            "Phosphatidylcholine": "cornflowerblue",
            "Ether/plasmalogen-PC": "cornflowerblue",
            "Phosphatidylethanolamine": "cornflowerblue",
            "Ether/plasmalogen-PE": "cornflowerblue",
            "Phosphatidylserine": "cornflowerblue",
            "Phosphatidylglycerol": "cornflowerblue",
            "Phosphatidylinositol": "cornflowerblue",
            "Lyso-PC": "teal",
            "Lyso-PE": "teal",
            "Sphingomyelin": "saddlebrown",
            "dhCer": "saddlebrown",
            "Ceramides": "saddlebrown",
            "dh(deoxy)": "saddlebrown",
            "Cer": "saddlebrown",
            "dh(deoxy)Cer": "saddlebrown",
            "(deoxy)Ceramides": "saddlebrown",
            "phytoCer": "saddlebrown",
            "MonohexosylCer": "saddlebrown",
            "DihexosylCer": "saddlebrown",
            "TrihexosylCer": "saddlebrown",
            "Triacylglycerols": "firebrick",
            "DihydroCer": "saddlebrown",
            "Hex(n)Cer": "saddlebrown",
            "DhDeoxyCer": "saddlebrown",
            "DeoxyCer": "saddlebrown",
            "Sphingadiene-Cer": "saddlebrown",
            "Sphingadienine-Cer": "saddlebrown",
            "Hex(n)-Sphingadiene-Cer": "saddlebrown",
            "PhytoCer": "saddlebrown",
            "Dihydro-Cer": "saddlebrown",
            "Dihydro-Deoxy-Cer": "saddlebrown",
            "Deoxy-Cer": "saddlebrown",
            "Phyto-Cer": "saddlebrown",
            "Hex-(n)-Cer": "saddlebrown",
            "Hex-(n)-Phyto-Cer": "saddlebrown",
            "Hex-(n)-Sphingadiene-Cer": "saddlebrown",
            "Hex-(n)-Sphingadienine-Cer": "saddlebrown",
        }

    def plot(self):

        fig, ax = plt.subplots(figsize=(15, 10))
        colors = []
        data = []
        sum_data = []
        sum_data_all = []
        label_data = []
        min_lst = []
        max_lst = []
        for lipid_class in self.lipid_class_lst:
            sub_df = self.df[self.df["Class"] == lipid_class]
            c_lst = sub_df["C"].tolist()
            # c_lst.append(sum(c_lst))
            data.append(c_lst)
            sum_data.append([sum(c_lst)])
            sum_data_all.append(sum(c_lst))
            label_data.append(lipid_class)
            colors.append(self.colors[lipid_class.strip("")])
            max_lst.append([max(c_lst), sub_df.at[sub_df["C"].idxmax(), "Lipid"]])
            min_lst.append([min(c_lst), sub_df.at[sub_df["C"].idxmin(), "Lipid"]])

        colors = list(reversed(colors))
        data = list(reversed(data))
        sum_data = list(reversed(sum_data))
        sum_data_all = list(reversed(sum_data_all))
        label_data = list(reversed(label_data))
        max_lst = list(reversed(max_lst))
        min_lst = list(reversed(min_lst))

        ax.eventplot(
            data, lineoffsets=2, linelengths=1, linewidths=0.9, colors=colors, alpha=0.9
        )
        ax.eventplot(
            sum_data, lineoffsets=2, linelengths=1, linewidths=1.5, colors=colors, alpha=1
        )
        sum_all = list(zip(sum_data_all, colors))
        # for sp in sum_all:
        #     ax.scatter(
        #         [sp[0]],
        #         [-2],
        #         color=sp[1],
        #         alpha=0.9,
        #         marker="|",
        #         edgecolors=sp[1],
        #         s=228,
        #         linewidths=0.5,
        #     )
        #     # ax.text(sp[0], -3.5, f'{sp[0]:.1E}', color=sp[1], horizontalalignment='center', rotation=90, fontsize=5)

        y = -0.5
        idx = 0
        for l_lst in data:
            bg_box = patches.Rectangle(
                (min(l_lst), y),
                max(l_lst) - min(l_lst),
                1,
                facecolor=colors[idx],
                edgecolor="none",
                alpha=0.3,
                )
            plt.text(
                min(l_lst) * 0.9,
                y + 0.22,
                f"{label_data[idx]}  {len(l_lst)} ",
                horizontalalignment="right",
                )
            plt.text(
                sum_data[idx][0],
                y + 1.15,
                "sum",
                size=5,
                color=colors[idx],
                horizontalalignment="center",
                )
            # plt.text(sum_data[idx][0] * 1.15, y+0.25, f'{len(l_lst)}', horizontalalignment='left', size=8, color=colors[idx])
            plt.text(
                min_lst[idx][0] * 0.9,
                y + 1.15,
                min_lst[idx][1],
                size=5,
                color=colors[idx],
                horizontalalignment="left",
                )
            plt.text(
                max_lst[idx][0] * 1.1,
                y + 1.15,
                max_lst[idx][1],
                size=5,
                color=colors[idx],
                horizontalalignment="right",
                )
            y += 2
            idx += 1
            ax.add_patch(bg_box)
        ax.get_yaxis().set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        plt.xscale("log")
        plt.grid(which="both", linestyle=":", linewidth=0.55, alpha=0.75)

        plt.savefig("img/distribution_plot.svg", dpi=300)
        plt.savefig("img/distribution_plot.png", dpi=300)
        plt.show()


if __name__ == "__main__":
    f = r"data/distribution_plot_data.xlsx"
    d = DistributionPlot(in_file=f)
    d.plot()
