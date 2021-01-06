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

import re

import pandas as pd
from natsort import natsorted

out_file = r'./output/DG_sum6_sorted_RT2.xlsx'

sum_file = r'./output/DG_sum6_sorted.xlsx'
l_file = r'./output/dg_summary.xlsx'

sum_df = pd.read_excel(sum_file)
l_df = pd.read_excel(l_file)

gl_rgx = re.compile(r'(?P<GL>\w{2,3})\((?P<FAs>.*)\)')

s_l_list = sum_df['Discrete'].unique().tolist()

l_info_dct = {
    'Leipzig': {
        'polar': 'Polar_L',
        'unpolar': 'Unpolar_L',
    },
    'Bremen': {
        'polar': 'Polar_IT_L',
        'Aquire': 'AquireX_L',
        'all_extract': 'Unpolar_all_extract_L',
        'conc': 'Unpolar_all_extract_conc_L',
    }
}


for li, lr in l_df.iterrows():
    l_p = lr['Discrete']
    l_l = lr['location']
    l_f = lr['s_fraction']
    l_rt = lr['AVG_RT']
    # l_mz = lr['Lib_mz']

    l_rt = round(l_rt, 2)
    # l_mz = f'.2f{l_mz}'

    l_match = re.match(gl_rgx, l_p)
    if l_match:
        h_info_dct = l_match.groupdict()
        h = f"{h_info_dct['GL']}({'_'.join(natsorted(h_info_dct['FAs'].split('_')))})"
        if h in s_l_list:
            print(h, l_rt, l_l, l_f)
            sum_df[l_info_dct[l_l][l_f]].loc[sum_df['Discrete'] == h] = l_rt

elem_rgx = re.compile(r'C(?P<c>\d\d)H(?P<h>\d\d\d?)NO(?P<o>\d)\+')

rt_gp_dct = {
    'Polar': ['Polar_L', 'Polar_H', 'Polar_S'],
    'Polar_IT': ['Polar_IT_L', 'Polar_IT_H', 'Polar_IT_S'],
    'Unpolar': ['Unpolar_L', 'Unpolar_H', 'Unpolar_S'],
    'AquireX': ['AquireX_L', 'AquireX_H', 'AquireX_S'],
    'all_extract': ['Unpolar_all_extract_L', 'Unpolar_all_extract_conc_L'],
}

for rt in rt_gp_dct:
    sum_df[f'RT_{rt}'] = sum_df[rt_gp_dct[rt]].mean(axis=1)
    sum_df[f'{rt}_delta'] = 0.0
    # sum_df[f'RT_{rt}'].round(2)

for i, r in sum_df.iterrows():
    formula = r['Formula_Charged']
    f_match = re.match(elem_rgx, formula)

    if f_match:
        f_dct = f_match.groupdict()
        # print(f_dct)
        mz = int(f_dct['c']) * 12 + int(f_dct['h']) * 1.0078250321 + 14.0030740052 + int(f_dct['o']) * 15.9949146221
        sum_df.at[i, 'Lib_mz'] = mz
        for rt in rt_gp_dct:
            rt_avg = r[f'RT_{rt}']

            if rt_avg > 0:
                rt_d = 0
                rt_col_lst = rt_gp_dct[rt]
                for col in rt_col_lst:
                    rt_dx = abs(r[col] - rt_avg)
                    if rt_dx > rt_d:
                        rt_d = rt_dx
                sum_df.at[i, f'{rt}_delta'] = rt_d

    l_p = r['Discrete']
    if '-' in l_p:
        l_match = re.match(gl_rgx, l_p)
        if l_match:
            h_info_dct = l_match.groupdict()
            fa_lst = h_info_dct['FAs'].split('_')
            fa_a_lst = []
            fa_op = ''
            for fa in fa_lst:
                if '-' in fa:
                    fa_op = fa
                else:
                    fa_a_lst.append(fa)
            h = f"{h_info_dct['GL']}({fa_op}_{'_'.join(natsorted(fa_a_lst))})"
            sum_df.at[i, 'Discrete'] = h

sum_df['Lib_mz'] = sum_df['Lib_mz'].round(decimals=5)
sum_df.drop_duplicates(inplace=True, keep='first')

idx_lst =[]
out_df = pd.DataFrame()
for i, r in sum_df.iterrows():

    d_lst = []

    for d in ['Polar_delta', 'Polar_IT_delta', 'Unpolar_delta', 'AquireX_delta']:
        if r[d] > 0:

            d_lst.append(r[d])
    if d_lst:
        if min(d_lst) < 0.3:
            # print(d_lst)
            # out_df = out_df.append(sum_df, sort=False)
            idx_lst.append(i)
    else:
        # print(d_lst)
        # print(r)
        idx_lst.append(i)
        # out_df = out_df.append(sum_df, sort=False)

out_df = sum_df.loc[idx_lst, :]

out_df.to_excel(out_file)
print('FIN')
