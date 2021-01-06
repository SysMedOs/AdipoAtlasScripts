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

ls_tg_aq = r'Data/Lipostar/Lipostar_TG_AquireX.xlsx'
ls_tg_po = r'Data/Lipostar/Lipostar_TG_polar.xlsx'
ls_tg_up = r'Data/Lipostar/Lipostar_TG_unpolar.xlsx'
ls_tg_merged = r'Data/data/Lipostar_TGDG_all2.xlsx'

header = ['Class', 'Bulk', 'Discrete', 'Formula', 'Charge']
# sum_header = header.copy()

rn_col_lst = ['AVG_RT', 'MIN_RT', 'MAX_RT', 'SIR', 'SIS', 'VIR', 'VIS']
groups_lst = ['SIR', 'SIS', 'VIR', 'VIS']

file_dct = {'P': ls_tg_po, 'U': ls_tg_up, 'A': ls_tg_aq}
f_col_lst = []
m_df = pd.DataFrame()

for f in file_dct:
    df = pd.read_excel(file_dct[f])
    # df = df[df['[M+NH4]+'] == 'T']
    df.at[:, 'Charge'] = '[M+NH4]+'
    df = df[header+ rn_col_lst]
    tmp_rn_col_dct = {}
    for col in rn_col_lst:
        tmp_rn_col_dct[col] = f'{col}_{f}'
        f_col_lst.append(f'{col}_{f}')

    df.rename(columns=tmp_rn_col_dct, inplace=True)

    if not m_df.empty:
        m_df = m_df.merge(df, how='outer', on=['Class', 'Bulk', 'Discrete', 'Formula', 'Charge'])
    else:
        m_df = df.copy()
for col in groups_lst:
    m_df[col] = ''

m_df = m_df[header + f_col_lst + groups_lst].fillna(0)
m_df['SUM_IDENT'] = ''
for g in groups_lst:
    for f in file_dct:
        m_df[f'{g}_{f}'] = m_df[f'{g}_{f}'].astype(int)
        m_df[f'{g}_{f}_T'] = m_df[f'{g}_{f}'].astype(str)
        m_df[f'{g}_{f}_T'] = m_df[f'{g}_{f}_T'].str.replace(r'0', '')
        m_df[f'{g}_{f}_T'] = m_df[f'{g}_{f}_T'].str.replace(r'\d', f)
        m_df[g] = m_df[g] + m_df[f'{g}_{f}_T']
    m_df[g] = m_df[g].apply(lambda x: ''.join(sorted(list(x))))
    m_df['SUM_IDENT'] = m_df['SUM_IDENT'] + m_df[g]

m_df['SUM_IDENT'] = m_df['SUM_IDENT'].apply(lambda x: ''.join(sorted(list(set(list(x))))))

m_df = m_df[header + f_col_lst + groups_lst + ['SUM_IDENT']]

m_df.to_excel(ls_tg_merged)

print(m_df.shape)
print('FIN')

