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

# %%
import re

from natsort import natsorted
import pandas as pd


# %%

# l_file = r'./data/PL/PC_manual.xlsx'
l_file = r'./data/PL/PS_manual.xlsx'
h_file = r'./data/PL/LipidHunter_PL_sum3.xlsx'
s_file = r'./data/PL/Lipostar_PL_sum2.xlsx'

pl_class = 'PS'

out_file = r'./output/{pl}_sum4.xlsx'.format(pl=pl_class)


l_df = pd.read_excel(l_file)
h_df = pd.read_excel(h_file)
s_df = pd.read_excel(s_file)
l_df['LipidSearch'] = 'LipidSearch'
h_df['LipidHunter'] = 'LipidHunter'
s_df['Lipostar'] = 'Lipostar'

m_df = l_df

m_df = m_df.merge(h_df, how='outer', on=['Class', 'Bulk', 'Discrete'],
                  suffixes=['_Search', f'_Hunter'])
m_df = m_df.merge(s_df, how='outer', on=['Class', 'Bulk', 'Discrete'],
                  suffixes=['', f'_Star'])

m_df = m_df[m_df['Class'] == pl_class]
m_df.to_excel(out_file)

print('FIN')

# %%

l_file = r'./sum/data/DG/DG_manual.xlsx'
h_file = r'./sum/data/DG/LipidHunter_GL_sum3.xlsx'
s_file = r'./sum/data/DG/Lipostar_DG.xlsx'

lp_class = 'DG'

out_file = r'./sum/output/{lp}_sum6.xlsx'.format(lp=lp_class)


l_df = pd.read_excel(l_file)
h_df = pd.read_excel(h_file)
s_df = pd.read_excel(s_file)
l_df['LipidSearch'] = 'LipidSearch'
h_df['LipidHunter'] = 'LipidHunter'
s_df['Lipostar'] = 'Lipostar'

m_df = l_df

m_df = m_df.merge(h_df, how='outer', on=['Class', 'Bulk', 'Discrete'],
                  suffixes=['_Search', f'_Hunter'])
m_df = m_df.merge(s_df, how='outer', on=['Class', 'Bulk', 'Discrete'],
                  suffixes=['', f'_Star'])

m_df = m_df[m_df['Class'] == lp_class]
m_df.to_excel(out_file)

print('FIN')


# %%

l_file = r'./data/TG/TG_manual2.xlsx'
h_file = r'./data/TG/LipidHunter_TG.xlsx'
s_file = r'./data/TG/Lipostar_GL_sum4.xlsx'

lp_class = 'TG'

out_file = r'./output/{c}_sum4.xlsx'.format(c=lp_class)


l_df = pd.read_excel(l_file)
h_df = pd.read_excel(h_file)
s_df = pd.read_excel(s_file)

for i, r in l_df.iterrows():
    u = r['Discrete']
    if isinstance(u, str):
        if u.startswith('TG'):
            print(u)
            u_core = u[2:].strip('()')
            fa = '_'.join(natsorted(u_core.split('_')))
            u_fa = f'TG({fa})'
            if u != u_fa:
                print('L', u, '->', u_fa)
                l_df.at[i, 'Discrete'] = u_fa

for i, r in h_df.iterrows():
    u = r['Discrete']
    if isinstance(u, str):
        if u.startswith('TG'):
            print(u)
            u_core = u[2:].strip('()')
            fa = '_'.join(natsorted(u_core.split('_')))
            u_fa = f'TG({fa})'
            if u != u_fa:
                print('H', u, '->', u_fa)
                h_df.at[i, 'Discrete'] = u_fa

for i, r in s_df.iterrows():
    u = r['Discrete']
    if isinstance(u, str):
        if u.startswith('TG'):
            print(u)
            u_core = u[2:].strip('()')
            fa = '_'.join(natsorted(u_core.split('_')))
            u_fa = f'TG({fa})'
            if u != u_fa:
                print('S', u, '->', u_fa)
                s_df.at[i, 'Discrete'] = u_fa

print(l_df.head())
print(h_df.head())
print(s_df.head())

l_df['LipidSearch'] = 'LipidSearch'
h_df['LipidHunter'] = 'LipidHunter'
s_df['Lipostar'] = 'Lipostar'

m_df = l_df

m_df = m_df.merge(h_df, how='outer', on=['Class', 'Bulk', 'Discrete'],
                  suffixes=['_Search', f'_Hunter'])
m_df = m_df.merge(s_df, how='outer', on=['Class', 'Bulk', 'Discrete'],
                  suffixes=['', f'_Star'])

m_df = m_df[m_df['Class'] == lp_class]
m_df.to_excel(out_file)

print('FIN')
