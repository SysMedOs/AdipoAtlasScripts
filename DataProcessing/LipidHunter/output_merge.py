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

#%%

merged_output = r'output/LipidHunter_PL_sum3.xlsx'

f_info_dct = {
    'Polar_Neg_PL': r'output/LipidHunter_Polar_Neg_PL3.xlsx',
    'Polar_Neg_PL_IT': r'output/LipidHunter_Polar_Neg_PL_IT3.xlsx',
}

main_col_lst = [
    'Class', 'Bulk', 'Discrete', 'Formula', 'Formula_charged',
    'Charge', 'Lib_mz', 'AVG_RT', 'MIN_RT', 'MAX_RT',
]

m_df = pd.DataFrame()

for file in f_info_dct:
    print(file)
    s_df = pd.read_excel(f_info_dct[file])
    s_df = s_df[main_col_lst]
    if m_df.empty:
        m_df = s_df
    else:
        m_df = m_df.merge(s_df, how='outer', on=['Class', 'Bulk', 'Discrete', 'Formula', 'Formula_charged', 'Charge'],
                          suffixes=['', f'_{file}'])

m_df.to_excel(merged_output)


#%%

merged_output = r'output/LipidHunter_GL_sum3.xlsx'

f_info_dct = {
    'Polar': r'output/GL/LipidHunter_TGDG_polar3.xlsx',
    'Polar_IT': r'output/GL/LipidHunter_TGDG_polar_IT3.xlsx',
    'Unpolar': r'output/GL/LipidHunter_TGDG_unpolar3.xlsx',
    'AquireX': r'output/GL/LipidHunter_TGDG_AquireX3.xlsx',
}

main_col_lst = [
    'Class', 'Bulk', 'Discrete', 'Formula', 'Formula_charged',
    'Charge', 'Lib_mz', 'AVG_RT', 'MIN_RT', 'MAX_RT',
]

m_df = pd.DataFrame()

for file in f_info_dct:
    print(file)
    s_df = pd.read_excel(f_info_dct[file])
    s_df = s_df[main_col_lst]
    if m_df.empty:
        m_df = s_df
    else:
        m_df = m_df.merge(s_df, how='outer', on=['Class', 'Bulk', 'Discrete', 'Formula', 'Formula_charged', 'Charge'],
                          suffixes=['', f'_{file}'])

m_df.to_excel(merged_output)
