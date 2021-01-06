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

import os
import re

from natsort import natsorted
import numpy as np
import pandas as pd


class LipidSearchParser:

    def __init__(self, lab_lst, fraction_lst, charge_lst, header_info, file_groups, fa_check=True):
        self.lab_lst, self.fraction_lst, self.charge_lst = lab_lst, fraction_lst, charge_lst
        self.charge_dict = {'M-H': '[M-H]-', 'M+HCOO': '[M+HCOO]-',
                            'M+H': '[M+H]+', 'M+NH4': '[M+NH4]+', 'M+Na': '[M+Na]+'}

        header_df = pd.read_excel(header_info)
        self.header_lst = header_df['Unified'].values.tolist()
        self.header_dct = dict(zip(header_df['LipidSearch'].values.tolist(), self.header_lst))

        file_df = pd.read_excel(file_groups, index_col=0)
        self.file_dct = file_df.to_dict(orient='index')

        self.file_abbr_lst = []
        self.abbr_group_dct = {}
        self.group_abbr_dct = {}
        for f in self.file_dct:
            self.file_abbr_lst.append(self.file_dct[f]['ABBR'])
            self.abbr_group_dct[self.file_dct[f]['ABBR']] = self.file_dct[f]['GROUP']
            if self.file_dct[f]['GROUP'] in self.group_abbr_dct:
                self.group_abbr_dct[self.file_dct[f]['GROUP']].append(self.file_dct[f]['ABBR'])
            else:
                self.group_abbr_dct[self.file_dct[f]['GROUP']] = [self.file_dct[f]['ABBR']]
        self.headers = ['Class', 'Bulk', 'Discrete', 'Formula', 'Formula_charged',
                        'Charge', 'Lib_mz', 'AVG_M_SCORE', 'AVG_T_SCORE', 'AVG_ppm', 'AVG_RT', 'MIN_RT',
                        'MAX_RT', 'IDENT_PEAK', 'IDENT_HG', 'IDENT_COUNT']
        self.headers.extend(sorted(list(self.group_abbr_dct.keys())))
        self.headers.extend(sorted(list(self.abbr_group_dct.keys())))
        self.fa_check = fa_check

    @staticmethod
    def merge_inputs(data_folder):

        file_lst = [os.path.join(data_folder, f) for f in os.listdir(data_folder)
                    if os.path.isfile(os.path.join(data_folder, f))]

        x_lst = [x for x in file_lst if x.endswith('.xlsx')]

        m_df = pd.DataFrame()
        for x in x_lst:
            tmp_df = pd.read_excel(x)
            m_df = m_df.append(tmp_df)

        return m_df

    def format_df(self, m_df: pd.DataFrame):

        m_df = m_df.replace(self.charge_dict)

        m_df.rename(columns=self.header_dct, inplace=True)
        m_header = m_df.columns.tolist()
        rest_header_lst = [h for h in m_header if h not in self.headers]
        header_lst = self.headers + sorted(rest_header_lst)

        for h in header_lst:
            if h not in m_header:
                m_df[h] = ''
        try:
            m_df = m_df[m_df['location'].isin(self.lab_lst)]
        except Exception as e:
            print(e)
        m_df = m_df[m_df['s_fraction'].isin(self.fraction_lst)]
        m_df = m_df[m_df['Charge'].isin(self.charge_lst)]

        for f in self.fraction_lst:
            m_df[f] = ''

        m_df.reset_index(drop=True, inplace=True)

        m_df['Class'] = m_df['Bulk'].str.replace(r'\(.*', '')
        m_df['Formula'] = m_df['Formula'].str.replace(r'N1', 'N')
        m_df['Formula'] = m_df['Formula'].str.replace(r'P1', 'P')
        m_df['Formula'] = m_df['Formula'].str.replace(r'N0', '')
        m_df['Formula'] = m_df['Formula'].str.replace(r'P0', '')
        m_df['Formula'] = m_df['Formula'].str.replace(r' ', '')
        m_df['Formula_charged'] = m_df['Formula_charged'].str.replace(r'N1', 'N')
        m_df['Formula_charged'] = m_df['Formula_charged'].str.replace(r'P1', 'P')
        m_df['Formula_charged'] = m_df['Formula_charged'].str.replace(r'N0', '')
        m_df['Formula_charged'] = m_df['Formula_charged'].str.replace(r'P0', '')

        try:
            m_df['Discrete'] = m_df['LipidIon'].str.replace(r'\).*', ')')
        except KeyError:
            print(m_df.head())

        m_df = m_df[header_lst]

        decimal_dct = {'Lib_mz': 2, 'AVG_RT': 2, 'MIN_RT': 2, 'MAX_RT': 2}
        m_df = m_df.round(decimal_dct)

        unique_df = pd.DataFrame()

        for i, r in m_df.iterrows():
            u = r['Discrete']
            if isinstance(u, str):
                if u.startswith('TG'):
                    print(u)
                    u_core = u[2:].strip('()')
                    fa = '_'.join(natsorted(u_core.split('_')))
                    u_fa = f'TG({fa})'
                    if u != u_fa:
                        print(u, '->', u_fa)
                        m_df.at[i, 'Discrete'] = u_fa

        unique_discrete_lst = m_df['Discrete'].unique().tolist()
        print(unique_discrete_lst)
        for u in unique_discrete_lst:
            tmp_df = m_df[m_df['Discrete'] == u]
            if not tmp_df.empty:
                try:
                    tmp_df['IDENT_PEAK'] = tmp_df['ProductIon'].str.count('FA')
                except KeyError:
                    pass
                r_df = tmp_df.head(1)
                try:
                    file_lst = tmp_df['s_name'].unique().tolist()
                    for file in file_lst:
                        if file in self.file_abbr_lst:
                            r_df.at[:, file] = 1
                except KeyError:
                    pass

                for gp in self.group_abbr_dct:
                    r_df.at[:, gp] = r_df[self.group_abbr_dct[gp]].sum(axis=1)
                try:
                    r_df.at[:, 'Lib_mz'] = tmp_df['CalcMz'].mean()
                    r_df.at[:, 'AVG_M_SCORE'] = tmp_df['m.Score'].mean()
                    r_df.at[:, 'AVG_T_SCORE'] = tmp_df['t.Score'].mean()
                    r_df.at[:, 'AVG_ppm'] = tmp_df['Delta.PPM.'].mean()
                    r_df.at[:, 'AVG_RT'] = tmp_df['Rt'].mean()
                    r_df.at[:, 'MIN_RT'] = tmp_df['Rt'].min()
                    r_df.at[:, 'MAX_RT'] = tmp_df['Rt'].max()
                    r_df.at[:, 'IDENT_PEAK'] = tmp_df['IDENT_PEAK'].max()
                    r_df.at[:, 'IDENT_COUNT'] = r_df[list(self.group_abbr_dct.keys())].sum(axis=1)
                except KeyError:
                    pass

                found_fraction_lst = tmp_df['s_fraction'].unique().tolist()
                for f in found_fraction_lst:
                    if f in self.fraction_lst:
                        m_df[f] = 'T'

                print(r_df['Discrete'])
                print(r_df[['IDENT_COUNT'] + list(self.group_abbr_dct.keys()) + ['IDENT_PEAK']])
                if self.fa_check:
                    if u.startswith('TG'):
                        print(f'{u} is TG')
                        unique_fa_lst = r_df['fa1'].tolist() + r_df['fa2'].tolist() + r_df['fa3'].tolist()
                        if r_df['IDENT_PEAK'].min() >= len(set(unique_fa_lst)):
                            unique_df = unique_df.append(r_df)
                        else:
                            print(f'WARNING: Not enough peaks identified - {u}')
                    elif u.startswith('DG') or u.startswith('PL'):
                        print(f'{u} is DG or PL')
                        unique_fa_lst = r_df['fa1'].tolist() + r_df['fa2'].tolist()
                        if r_df['IDENT_PEAK'].min() >= len(set(unique_fa_lst)):
                            unique_df = unique_df.append(r_df)
                        else:
                            print(f'WARNING: Not enough peaks identified - {u}')
                    elif u.startswith('MG') or u.startswith('LPL'):
                        if r_df['IDENT_PEAK'].min() >= 1:
                            unique_df = unique_df.append(r_df)
                        else:
                            print(f'WARNING: Not enough peaks identified - {u}')
                    else:
                        print(f'WARNING: Not supported Lipid class - {u}')
                else:
                    unique_df = unique_df.append(r_df)

        # unique_df['IDENT_COUNT'] = unique_df[list(self.group_abbr_dct.keys())].sum(axis=1)
        print(unique_df.head())

        for h in header_lst:
            if re.match(r'\w\w\w\d', h):
                # m_df[h] = m_df[h].astype
                try:
                    unique_df[h] = unique_df[h].replace(range(1, 100), 'T')
                except KeyError:
                    pass

        for idx, r in unique_df.iterrows():
            if r['Discrete'] is None or r['Discrete'] in ['', np.nan]:
                print(r['Bulk'], r['Discrete'])
                unique_df.at[idx, 'Discrete'] = r['Bulk']
        try:
            unique_df = unique_df[header_lst]
        except KeyError:
            pass
        decimal_lst = {'AVG_M_SCORE': 2, 'AVG_T_SCORE': 2, 'AVG_ppm': 2, 'AVG_RT': 2, 'MIN_RT': 2, 'MAX_RT': 2}
        unique_df = unique_df.round(decimal_lst)

        unique_df = unique_df.sort_values(by=['Class', 'Lib_mz', 'Discrete'])
        unique_df.reset_index(inplace=True, drop=True)

        return unique_df

    def merge_info(self, data_folder, output_path):

        m_df = self.merge_inputs(data_folder)
        m_df = self.format_df(m_df)

        m_df.to_excel(output_path)


if __name__ == '__main__':

    # usr_header_info = r'../Configurations/Headers.xlsx'
    # usr_data_folder = r'data/PL'
    # usr_file_groups = r'../Configurations/file_groups_Polar_Neg.xlsx'
    # usr_output_path = r'output/PC_manual.xlsx'
    # usr_lab_lst = ['Leipzig']
    # usr_fraction_lst = ['polar']
    # usr_charge_lst = ['[M-H]-', '[M+HCOO]-', '[M+H]+']

    # usr_output_path = r'Data/LipidSearch/TG_polar.xlsx'
    # usr_fraction_lst = ['polar']
    # usr_output_path = r'Data/LipidSearch/TG_unpolar.xlsx'
    # usr_fraction_lst = ['unpolar']

    # usr_data_folder = r'Data/LipidSearch/TG_AquireX'
    # usr_output_path = r'Data/LipidSearch/TG_AquireX.xlsx'
    # usr_file_groups = r'Configurations/file_groups_AquireX.xlsx'
    # usr_fraction_lst = ['aq']
    # usr_lab_lst = ['Leipzig']

    # usr_header_info = r'../Configurations/Headers.xlsx'
    # usr_data_folder = r'data/TG'
    # usr_file_groups = r'../Configurations/file_groups_Polar_Neg.xlsx'
    # usr_output_path = r'output/TG_manual2.xlsx'
    # usr_lab_lst = ['Leipzig', 'Bremen']
    # usr_fraction_lst = ['polar', 'unpolar', 'Aquire', 'all_extract', 'conc']
    # usr_charge_lst = ['[M+NH4]+']

    # usr_header_info = r'../Configurations/Headers.xlsx'
    # usr_data_folder = r'data/DG'
    # usr_file_groups = r'../Configurations/file_groups_Polar_Neg.xlsx'
    # usr_output_path = r'output/DG_manual.xlsx'
    # usr_lab_lst = ['Leipzig', 'Bremen']
    # usr_fraction_lst = ['polar', 'unpolar', 'Aquire', 'all_extract', 'conc']
    # usr_charge_lst = ['[M+NH4]+']

    # usr_header_info = r'../Configurations/Headers.xlsx'
    # usr_data_folder = r'data/PE'
    # usr_file_groups = r'../Configurations/file_groups_Polar_Neg.xlsx'
    # usr_output_path = r'output/PE_manual.xlsx'
    # usr_lab_lst = ['Leipzig', 'Bremen']
    # usr_fraction_lst = ['polar', 'unpolar']
    # usr_charge_lst = ['[M-H]-', '[M+HCOO]-', '[M+H]+']

    usr_header_info = r'../Configurations/Headers.xlsx'
    usr_data_folder = r'data/PS'
    usr_file_groups = r'../Configurations/file_groups_Polar_Neg.xlsx'
    usr_output_path = r'output/PS_manual.xlsx'
    usr_lab_lst = ['Leipzig', 'Bremen']
    usr_fraction_lst = ['polar', 'unpolar']
    usr_charge_lst = ['[M-H]-', '[M+HCOO]-', '[M+H]+']
    lsp = LipidSearchParser(usr_lab_lst, usr_fraction_lst, usr_charge_lst,
                            usr_header_info, usr_file_groups, fa_check=False)
    lsp.merge_info(usr_data_folder, output_path=usr_output_path)

    usr_header_info = r'../Configurations/Headers.xlsx'
    usr_data_folder = r'data/PG'
    usr_file_groups = r'../Configurations/file_groups_Polar_Neg.xlsx'
    usr_output_path = r'output/PG_manual.xlsx'
    usr_lab_lst = ['Leipzig', 'Bremen']
    usr_fraction_lst = ['polar', 'unpolar']
    usr_charge_lst = ['[M-H]-', '[M+HCOO]-', '[M+H]+']
    lsp = LipidSearchParser(usr_lab_lst, usr_fraction_lst, usr_charge_lst,
                            usr_header_info, usr_file_groups, fa_check=False)
    lsp.merge_info(usr_data_folder, output_path=usr_output_path)

    usr_header_info = r'../Configurations/Headers.xlsx'
    usr_data_folder = r'data/PI'
    usr_file_groups = r'../Configurations/file_groups_Polar_Neg.xlsx'
    usr_output_path = r'output/PI_manual.xlsx'
    usr_lab_lst = ['Leipzig', 'Bremen']
    usr_fraction_lst = ['polar', 'unpolar']
    usr_charge_lst = ['[M-H]-', '[M+HCOO]-', '[M+H]+']
    lsp = LipidSearchParser(usr_lab_lst, usr_fraction_lst, usr_charge_lst,
                            usr_header_info, usr_file_groups, fa_check=False)
    lsp.merge_info(usr_data_folder, output_path=usr_output_path)

    print('FIN')
