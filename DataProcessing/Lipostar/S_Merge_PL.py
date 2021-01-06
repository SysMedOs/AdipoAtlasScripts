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

from natsort import natsorted
import pandas as pd


class LipostarParser:

    def __init__(self, sum_csv, top_csv, dbi_csv, file_groups, header_info, adduct_lst=[]):
        self.sum_df = pd.read_csv(sum_csv)
        self.top_df = pd.read_csv(top_csv)
        self.dbi_df = pd.read_csv(dbi_csv)

        # if adduct_lst:
        #     self.top_df = self.top_df.loc[self.top_df['Adduct'].isin(adduct_lst)]
        #     print(self.top_df.head())

        self.adduct_lst = adduct_lst

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

        self.headers = ['Lipid_Class', 'BULK', 'DISCRETE', 'Adduct', 'Confidence', 'Formula',
                        'MZ', 'RT', 'Score', 'Iso. Pat.Score', 'Num Frag Matches', 'IDENT_COUNT']
        self.file_abbr_lst = sorted(self.file_abbr_lst)
        self.headers.extend(sorted(list(self.group_abbr_dct.keys())))
        self.headers.extend(self.file_abbr_lst)

        header_df = pd.read_excel(header_info)
        self.header_lst = header_df['Unified'].values.tolist()
        self.header_dct = dict(zip(header_df['Lipostar'].values.tolist(), self.header_lst))

        print(self.group_abbr_dct)

    @staticmethod
    def get_sorted_abbr(discrete_abbr):

        print(discrete_abbr)

        bulk_abbr = None

        pl_discrete_rgx = re.compile(r'(?P<SUB_CLASS>\w+)\('  # HG
                                     r'(?P<LINK>O-|P-)?(?P<FA1>\d{1,2}):(?P<DB1>\d)/'  # FA1
                                     r'(?P<FA2>\d{1,2}):(?P<DB2>\d)\)')  # FA2

        tg_discrete_rgx = re.compile(r'(?P<SUB_CLASS>\w+)\('  # class
                                     r'(?P<LINK>O-|P-)?(?P<FA1>\d{1,2}):(?P<DB1>\d)(\([\dEZ,]*\))?/'  # FA1
                                     r'(?P<FA2>\d{1,2}):(?P<DB2>\d)(\([\dEZ,]*\))?/'  # FA2
                                     r'(?P<FA3>\d{1,2}):(?P<DB3>\d)(\([\dEZ,]*\))?\)')  # FA3

        if discrete_abbr.startswith('LP'):
            bulk_abbr = discrete_abbr

        else:

            pl_match = re.match(pl_discrete_rgx, discrete_abbr)
            tg_match = re.match(tg_discrete_rgx, discrete_abbr)

            if pl_match:
                pl_groups = pl_match.groupdict()
                if pl_groups['LINK']:
                    pl_link = pl_groups['LINK']
                    discrete_abbr = '{sub_class}({link}{pl_c1}:{pl_db1}_{pl_c2}:{pl_db2})'.format(
                        sub_class=pl_groups['SUB_CLASS'],
                        link=pl_link,
                        pl_c1=pl_groups['FA1'],
                        pl_db1=pl_groups['DB1'],
                        pl_c2=pl_groups['FA2'],
                        pl_db2=pl_groups['DB2'],
                    )
                else:
                    pl_link = ''
                    fa_lst = [
                        '{c}:{db}'.format(c=pl_groups['FA1'], db=pl_groups['DB1']),
                        '{c}:{db}'.format(c=pl_groups['FA2'], db=pl_groups['DB2'])
                    ]
                    fa_lst = natsorted(fa_lst)
                    discrete_abbr = '{sub_class}({fa1}_{fa2})'.format(sub_class=pl_groups['SUB_CLASS'], fa1=fa_lst[0],
                                                                      fa2=fa_lst[1])
                pl_c = int(pl_groups['FA1']) + int(pl_groups['FA2'])
                pl_db = int(pl_groups['DB1']) + int(pl_groups['DB2'])
                bulk_abbr = '{sub_class}({link}{pl_c}:{pl_db})'.format(sub_class=pl_groups['SUB_CLASS'], link=pl_link,
                                                                       pl_c=pl_c, pl_db=pl_db)
            elif tg_match:
                tg_groups = tg_match.groupdict()
                if tg_groups['LINK']:
                    tg_link = tg_groups['LINK']
                    discrete_abbr = '{sub_class}({link}{tg_c1}:{tg_db1}_{tg_c2}:{tg_db2})'.format(
                        sub_class=tg_groups['SUB_CLASS'],
                        link=tg_link,
                        tg_c1=tg_groups['FA1'],
                        tg_db1=tg_groups['DB1'],
                        tg_c2=tg_groups['FA2'],
                        tg_db2=tg_groups['DB2'],
                        tg_c3=tg_groups['FA3'],
                        tg_db3=tg_groups['DB3'],
                    )
                else:
                    tg_link = ''
                    fa_lst = [
                        '{c}:{db}'.format(c=tg_groups['FA1'], db=tg_groups['DB1']),
                        '{c}:{db}'.format(c=tg_groups['FA2'], db=tg_groups['DB2']),
                        '{c}:{db}'.format(c=tg_groups['FA3'], db=tg_groups['DB3']),
                    ]
                    fa_lst = natsorted(fa_lst)
                    discrete_abbr = '{sub_class}({fa1}_{fa2}_{fa3})'.format(sub_class=tg_groups['SUB_CLASS'],
                                                                            fa1=fa_lst[0], fa2=fa_lst[1], fa3=fa_lst[2])
                tg_c = int(tg_groups['FA1']) + int(tg_groups['FA2']) + int(tg_groups['FA3'])
                tg_db = int(tg_groups['DB1']) + int(tg_groups['DB2']) + int(tg_groups['DB3'])
                bulk_abbr = '{sub_class}({link}{tg_c}:{tg_db})'.format(sub_class=tg_groups['SUB_CLASS'], link=tg_link,
                                                                       tg_c=tg_c, tg_db=tg_db)

            else:
                print(f'Can not parse Lipid: {bulk_abbr} - {discrete_abbr}')

        return bulk_abbr, discrete_abbr

    def get_abbr(self, abbr):

        discrete_abbr = re.sub(r'\((\d{1,2}[ezEZ],?)(,\d{1,2}[ezEZ],?)*\)', '', abbr)

        left_lyso = re.compile(r'\(0:0/')
        right_lyso = re.compile(r'/0:0\)')

        lipid_class = re.sub(r'\(.*\)', '', discrete_abbr)

        left_match = left_lyso.search(discrete_abbr)
        right_match = right_lyso.search(discrete_abbr)
        if left_match:
            discrete_abbr = re.sub(r'\(0:0/', '(', discrete_abbr)
            if discrete_abbr.startswith('Lyso'):
                discrete_abbr = re.sub(r'Lyso', 'L', discrete_abbr)
            else:
                if lipid_class[0] == 'P':
                    discrete_abbr = f'L{discrete_abbr}'
        else:
            if right_match:
                discrete_abbr = re.sub(r'/0:0', '', discrete_abbr)
                if discrete_abbr.startswith('Lyso'):
                    discrete_abbr = re.sub(r'Lyso', 'L', discrete_abbr)
                else:
                    if lipid_class[0] == 'P':
                        discrete_abbr = f'L{discrete_abbr}'
        if discrete_abbr.startswith('Lyso'):
            discrete_abbr = re.sub(r'Lyso', 'L', discrete_abbr)

        (bulk_abbr, discrete_abbr) = self.get_sorted_abbr(discrete_abbr)

        abbr_dct = {'Lipid_Class': lipid_class, 'BULK': bulk_abbr, 'DISCRETE': discrete_abbr}

        return abbr_dct

    def sort_sum(self):

        sum_dct = {}
        current_feature_idx = None

        for idx, row in self.sum_df.iterrows():
            _cmp = row['Compound']
            _sample = row['Sample']
            if _sample == 'Super Sample' and _cmp[1] != r'>':
                _sup_feature_idx = _cmp.strip(r'>')
                _sup_feature = _sup_feature_idx.split(r'@')
                current_feature_idx = _sup_feature_idx
                current_feature = _sup_feature
                sum_dct[current_feature_idx] = {'MZ': float(current_feature[0]), 'RT': float(current_feature[1]),
                                                'Lipid_Class': row['Lipid Class']}
            else:
                if _sample in self.file_dct:
                    if row[r'Has MS/MS'] == 'Y':
                        sum_dct[current_feature_idx][self.file_dct[_sample]['ABBR']] = 1

        return sum_dct

    def sort_top(self, sum_dct):

        top_dct = {}
        current_feature_idx = None

        for idx, row in self.top_df.iterrows():
            _cmp = row['Compound']

            if _cmp in sum_dct:
                current_feature_idx = _cmp
            else:
                pass
            top_dct[idx] = row.to_dict()
            top_dct[idx].update(sum_dct[current_feature_idx])
            top_dct[idx].update(self.get_abbr(top_dct[idx]['Name']))

        return top_dct

    def merge_info(self, output_path, score=60, peaks=2):

        sum_dct = self.sort_sum()
        top_dct = self.sort_top(sum_dct)
        sorted_sum_df = pd.DataFrame(data=top_dct).T

        try:
            sorted_sum_df = sorted_sum_df.loc[sorted_sum_df['Adduct'].isin(self.adduct_lst)]
            sorted_sum_df = sorted_sum_df[sorted_sum_df['Chains'] != '-']
        except Exception as e:
            print(e)

        sum_df_header_lst = sorted_sum_df.columns.tolist()
        for gp in self.group_abbr_dct:
            tmp_header_lst = self.group_abbr_dct[gp].copy()
            h_lst = []
            for h in tmp_header_lst:
                if h in sum_df_header_lst:
                    h_lst.append(h)
                else:
                    print(f'{gp} - {h} not in header list')
                    sorted_sum_df[h] = 0
            sorted_sum_df.at[:, gp] = sorted_sum_df[h_lst].sum(axis=1)

        sorted_sum_df.at[:, 'IDENT_COUNT'] = sorted_sum_df[list(self.group_abbr_dct.keys())].sum(axis=1)

        sorted_sum_df = sorted_sum_df[sorted_sum_df['Fragment Score'] >= score]

        pre1_output_unique_df = sorted_sum_df[self.headers]
        pre1_output_unique_df = pre1_output_unique_df[pre1_output_unique_df['Num Frag Matches'] > 1]
        # decimal_lst = {'AVG_RANK_SCORE': 2, 'AVG_ISOTOPE_SCORE': 2, 'AVG_ppm': 2, 'AVG_RT': 2, 'MAX_RT': 2}
        # pre2_output_unique_df = pre1_output_unique_df.round(decimal_lst)
        pre1_output_unique_df.rename(columns=self.header_dct, inplace=True)  # type: pd.DataFrame
        pre1_output_unique_df.sort_values(by=['Class', 'Bulk', 'Discrete', 'Lib_mz', 'Charge',
                                              'Lib_mz', 'Score', 'Iso. Pat.Score',
                                              'Num Frag Matches', 'IDENT_COUNT'],
                                          ascending=[True, True, True, True, False,
                                                     True, False, False, False, False],
                                          inplace=True)

        pre1_output_unique_df = pre1_output_unique_df[pre1_output_unique_df['Num Frag Matches'] >= peaks]
        drop_idx_lst = []
        for i, r in pre1_output_unique_df.iterrows():
            if r['Class'] in ['PC', 'PE', 'PS'] and r['Num Frag Matches'] < peaks + 1:
                drop_idx_lst.append(i)
        pre1_output_unique_df = pre1_output_unique_df.drop(index=drop_idx_lst)

        unique_discrete_lst = pre1_output_unique_df['Discrete'].unique().tolist()
        print(len(unique_discrete_lst))
        unique_charge_lst = pre1_output_unique_df['Charge'].unique().tolist()
        if len(unique_charge_lst) > 1:
            output_unique_df = pd.DataFrame()
            for charge in unique_charge_lst:
                pre1_output_unique_df[charge] = ''
            for d in unique_discrete_lst:
                tmp_df = pre1_output_unique_df.query(f'Discrete == "{d}"')
                tmp_charge_lst = tmp_df['Charge'].unique().tolist()
                tmp_df['Charge'] = ','.join(natsorted(tmp_charge_lst, reverse=True))
                for chg in tmp_charge_lst:
                    tmp_df[chg] = 'T'
                output_unique_df = output_unique_df.append(tmp_df.head(1))
        else:
            output_unique_df = pre1_output_unique_df
        output_unique_df.drop_duplicates(subset=['Discrete', 'Charge'], keep='first', inplace=True)
        output_unique_df.reset_index(drop=True, inplace=True)

        m_header = output_unique_df.columns.tolist()
        rest_header_lst = [h for h in m_header if h not in self.header_lst]
        header_lst = self.header_lst + sorted(rest_header_lst)

        for h in header_lst:
            if h not in m_header:
                output_unique_df[h] = ''
        output_unique_df['Class'] = output_unique_df['Class'].str.replace(r'Lyso', 'L')
        output_unique_df['Bulk'] = output_unique_df['Bulk'].str.replace(r'Lyso', 'L')
        output_unique_df['Discrete'] = output_unique_df['Discrete'].str.replace(r'Lyso', 'L')
        output_unique_df = output_unique_df[header_lst]

        output_unique_df.to_excel(output_path)


if __name__ == '__main__':

    merged_output = r'../Lipostar/output/PL/Lipostar_PL_sum5.xlsx'

    sum_file_info_dct = {}
    usr_header_info = r'../Configurations/Headers.xlsx'

    # sum_file_info_dct['Polar_Neg_PL'] = {
    #     'file_groups': r'../Configurations/file_groups_Polar_Neg.xlsx',
    #     'sum_csv': r'../Lipostar/data/PL/Polar_Neg_PL_1.csv',
    #     'top_csv': r'../Lipostar/data/PL/Polar_Neg_PL_2.csv',
    #     'dbi_csv': r'../Lipostar/data/PL/Polar_Neg_PL_3.csv',
    #     'output_path': r'../Lipostar/output/PL/Lipostar_Polar_Neg_PL.xlsx',
    #     'addut_lst': ['[M-H]-', '[M+HCOO]-'],
    #     'score': 60, 'peaks': 2,
    # }
    # sum_file_info_dct['Polar_Pos_PL'] = {
    #     'file_groups': r'../Configurations/file_groups_Polar_Pos.xlsx',
    #     'sum_csv': r'../Lipostar/data/PL/Polar_Pos_PL_1.csv',
    #     'top_csv': r'../Lipostar/data/PL/Polar_Pos_PL_2.csv',
    #     'dbi_csv': r'../Lipostar/data/PL/Polar_Pos_PL_3.csv',
    #     'output_path': r'../Lipostar/output/PL/Lipostar_Polar_Pos_PL.xlsx',
    #     'addut_lst': ['[M+H]+'],
    #     'score': 60, 'peaks': 2,
    # }
    sum_file_info_dct['Polar_Neg_PL_IT'] = {
        'file_groups': r'../Configurations/file_groups_Polar_Neg_IT.xlsx',
        'sum_csv': r'../Lipostar/data/PL/Polar_Neg_PL_IT_1.csv',
        'top_csv': r'../Lipostar/data/PL/Polar_Neg_PL_IT_2.csv',
        'dbi_csv': r'../Lipostar/data/PL/Polar_Neg_PL_IT_3.csv',
        'output_path': r'../Lipostar/output/PL/Lipostar_Polar_Neg_PL_IT.xlsx',
        'addut_lst': ['[M-H]-', '[M+HCOO]-'],
        'score': 60, 'peaks': 2,
    }
    # sum_file_info_dct['Polar_Pos_PL_IT'] = {
    #     'file_groups': r'../Configurations/file_groups_Polar_Pos_IT.xlsx',
    #     'sum_csv': r'../Lipostar/data/PL/Polar_Pos_PL_IT_1.csv',
    #     'top_csv': r'../Lipostar/data/PL/Polar_Pos_PL_IT_2.csv',
    #     'dbi_csv': r'../Lipostar/data/PL/Polar_Pos_PL_IT_3.csv',
    #     'output_path': r'../Lipostar/output/PL/Lipostar_Polar_Pos_PL_IT.xlsx',
    #     'addut_lst': ['[M+H]+'],
    #     'score': 60, 'peaks': 2,
    # }
    # sum_file_info_dct['Unpolar_Pos_PL'] = {
    #     'file_groups': r'../Configurations/file_groups_Unpolar_Pos.xlsx',
    #     'sum_csv': r'../Lipostar/data/PL/Unpolar_Pos_PL_1.csv',
    #     'top_csv': r'../Lipostar/data/PL/Unpolar_Pos_PL_2.csv',
    #     'dbi_csv': r'../Lipostar/data/PL/Unpolar_Pos_PL_3.csv',
    #     'output_path': r'../Lipostar/output/PL/Unpolar_Pos_PL.xlsx',
    #     'addut_lst': ['[M+H]+'],
    #     'score': 60, 'peaks': 2,
    # }
    # sum_file_info_dct['Unpolar_Pos_PL_AquireX'] = {
    #     'file_groups': r'../Configurations/file_groups_Unpolar_Pos_AquireX.xlsx',
    #     'sum_csv': r'../Lipostar/data/PL/Unpolar_Pos_PL_AquireX_1.csv',
    #     'top_csv': r'../Lipostar/data/PL/Unpolar_Pos_PL_AquireX_2.csv',
    #     'dbi_csv': r'../Lipostar/data/PL/Unpolar_Pos_PL_AquireX_3.csv',
    #     'output_path': r'../Lipostar/output/PL/Unpolar_Pos_PL_AquireX.xlsx',
    #     'addut_lst': ['[M+H]+'],
    #     'score': 60, 'peaks': 2,
    # }

    m_df = pd.DataFrame()
    main_col_lst = ['Class', 'Bulk', 'Discrete', 'Formula', 'Formula_charged',
                    'Charge', 'Lib_mz', 'AVG_RT', 'MIN_RT', 'MAX_RT']

    for file in sum_file_info_dct:
        info_dct = sum_file_info_dct[file]
        print(file)
        lsp = LipostarParser(info_dct['sum_csv'], info_dct['top_csv'], info_dct['dbi_csv'],
                             info_dct['file_groups'], usr_header_info, info_dct['addut_lst'])
        lsp.merge_info(output_path=info_dct['output_path'], score=info_dct['score'], peaks=info_dct['peaks'])

    for file in sum_file_info_dct:
        info_dct = sum_file_info_dct[file]
        print(file)
        s_df = pd.read_excel(info_dct['output_path'])
        s_df = s_df[main_col_lst]
        if m_df.empty:
            m_df = s_df
        else:
            m_df = m_df.merge(s_df, how='outer', on=['Class', 'Bulk', 'Discrete', 'Formula', 'Formula_charged'],
                              suffixes=['', f'_{file}'])

    m_df.to_excel(merged_output)

    print('FIN')
