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

import logging
import os
import sys

from natsort import natsorted
import pandas as pd
import configparser


class HunterConfig:

    def __init__(self, file_groups, header_cfg, log_level=logging.DEBUG, hg_check=True, fa_check=True):

        logging.basicConfig(format='%(asctime)s-%(levelname)s - %(message)s',
                            datefmt='%b-%d@%H:%M:%S',
                            level=log_level)
        self.logger = logging.getLogger('log')

        self.opt_dct = {
            'lipid_class': 'lipid_class',
            'charge': 'charge_mode',
            'mzml': 'mzml_path_str',
            'xlsx': 'xlsx_output_path_str',
        }

        file_df = pd.read_excel(file_groups, index_col=0)
        self.file_dct = file_df.to_dict(orient='index')
        self.logger.debug(self.file_dct)

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
        self.headers = ['lipid_class', 'Proposed_structures', 'DISCRETE_ABBR', 'Formula_neutral', 'Formula_ion',
                        'Charge', 'Lib_mz', 'AVG_RANK_SCORE', 'AVG_ISOTOPE_SCORE', 'AVG_ppm', 'AVG_RT', 'MIN_RT',
                        'MAX_RT', 'IDENT_FA', 'IDENT_HG', 'IDENT_COUNT']
        self.headers.extend(sorted(list(self.group_abbr_dct.keys())))
        self.headers.extend(sorted(list(self.abbr_group_dct.keys())))

        header_df = pd.read_excel(header_cfg)
        self.header_lst = header_df['Unified'].values.tolist()
        self.header_dct = dict(zip(header_df['LipidHunter'].values.tolist(), self.header_lst))

        self.hg_check = hg_check
        self.fa_check = fa_check

    def check_file(self, file_path):
        is_file = False
        try:
            if os.path.isfile(file_path):
                is_file = True
                self.logger.debug(f'Load file: {file_path}')
            else:
                self.logger.error(f'Failed to load file: {file_path}')
        except IOError:
            self.logger.error(f'Failed to load file: {file_path}')

        return is_file

    def read_cfg(self, cfg):

        config = configparser.ConfigParser()
        config.read(cfg)

        if config.has_section('settings'):
            usr_cfg = 'settings'
        elif config.has_section('parameters'):
            usr_cfg = 'parameters'
        else:
            if config.has_section('default'):
                usr_cfg = 'default'
            else:
                usr_cfg = None

        cfg_dct = {}

        if usr_cfg is not None:
            options = config.options(usr_cfg)
            if self.opt_dct['xlsx'] in options:
                xlsx_path = config.get(usr_cfg, self.opt_dct['xlsx'])

                if self.check_file(xlsx_path):
                    self.logger.debug(xlsx_path)
                    cfg_dct['xlsx_path'] = xlsx_path
                    cfg_dct['xlsx'] = os.path.splitext(os.path.basename(xlsx_path))[0]

                    if self.opt_dct['lipid_class'] in options:
                        cfg_dct['lipid_class'] = config.get(usr_cfg, self.opt_dct['lipid_class'])

                    if self.opt_dct['charge'] in options:
                        cfg_dct['charge'] = config.get(usr_cfg, self.opt_dct['charge'])

                    if self.opt_dct['mzml'] in options:
                        mzml_path = config.get(usr_cfg, self.opt_dct['mzml'])
                        cfg_dct['mzml'] = os.path.splitext(os.path.basename(mzml_path))[0]

                else:
                    self.logger.error(f'Failed to load file: {xlsx_path}')

        if 'charge' in cfg_dct and 'xlsx' in cfg_dct:
            self.logger.info(cfg_dct)
            return cfg_dct

        else:
            return None

    def load_batch_cfg(self, cfg_lst, merge_table, rank_score=40, isotope_score=80):

        sum_cfg_dct = {}

        for cfg in cfg_lst:
            cfg_dct = self.read_cfg(cfg)
            try:
                sum_cfg_dct[cfg_dct['xlsx']] = cfg_dct
            except Exception as err:
                self.logger.error(err)

        cfg_df = pd.DataFrame(sum_cfg_dct).T

        sum_df = self.merge_xlsx(cfg_df)
        self.merge_features(sum_df, merge_table, rank_score=rank_score, isotope_score=isotope_score)

    def merge_xlsx(self, cfg_df):

        sum_df = pd.DataFrame()

        for i, r in cfg_df.iterrows():
            mzml = r['mzml']
            tmp_df = pd.read_excel(r['xlsx_path'])
            tmp_df['xlsx'] = r['xlsx']
            tmp_df['mzml'] = mzml
            if mzml in self.file_dct:
                mzml_abbr = self.file_dct[mzml]['ABBR']
                tmp_df[mzml_abbr] = 1
                tmp_df[self.file_dct[mzml]['GROUP']] = 1
            tmp_df['abs_ppm'] = tmp_df['ppm'].abs()

            tmp_df['lipid_class'] = r['lipid_class']

            if r['lipid_class'].startswith('P'):

                if r['lipid_class'] in ['PC', 'PE']:
                    _op_tmp_df = tmp_df.loc[tmp_df['Proposed_structures'].str.contains('-', regex=False)]
                    _op_tmp_df = _op_tmp_df.fillna(0)
                    op_tmp_df = _op_tmp_df[_op_tmp_df['FA1_[FA-H]-_i'] == 0]
                    op_tmp_df = op_tmp_df[op_tmp_df['FA2_[FA-H]-_i'] > 0]
                    op_tmp_df = op_tmp_df[op_tmp_df['[LPL(FA1)-H]-_i'] + op_tmp_df['[LPL(FA1)-H2O-H]-_i'] > 0]
                    sum_df = sum_df.append(op_tmp_df, sort=False)
                _fa_tmp_df = tmp_df[tmp_df['FA1_[FA-H]-_i'] > 0]
                fa_tmp_df = _fa_tmp_df[_fa_tmp_df['FA2_[FA-H]-_i'] > 0]
                sum_df = sum_df.append(fa_tmp_df, sort=False)

            elif r['lipid_class'].startswith('T'):

                _fa_tmp_df = tmp_df[tmp_df['[M-(FA3)+H]+_i'] > 0]
                _fa_tmp_df = _fa_tmp_df[_fa_tmp_df['[M-(FA2)+H]+_i'] > 0]

                _op_tmp_df = _fa_tmp_df.loc[_fa_tmp_df['Proposed_structures'].str.contains('-', regex=False)]
                _op_tmp_df = _op_tmp_df.fillna(0)
                op_tmp_df = _op_tmp_df[_op_tmp_df['[M-(FA1)+H]+_i'] == 0]
                sum_df = sum_df.append(op_tmp_df, sort=False)

                fa_tmp_df = _fa_tmp_df[_fa_tmp_df['[M-(FA1)+H]+_i'] > 0]
                sum_df = sum_df.append(fa_tmp_df, sort=False)

            elif r['lipid_class'].startswith('D'):

                _fa_tmp_df = tmp_df[tmp_df['Lib_mz'] < 750]
                _fa_tmp_df = _fa_tmp_df[_fa_tmp_df['[MG(FA2)-H2O+H]+_i'] > 0]

                _op_tmp_df = _fa_tmp_df.loc[_fa_tmp_df['Proposed_structures'].str.contains('-', regex=False)]
                _op_tmp_df = _op_tmp_df.fillna(0)
                op_tmp_df = _op_tmp_df[_op_tmp_df['[MG(FA1)-H2O+H]+_i'] == 0]
                sum_df = sum_df.append(op_tmp_df, sort=False)

                fa_tmp_df = _fa_tmp_df[_fa_tmp_df['[MG(FA1)-H2O+H]+_i'] > 0]
                sum_df = sum_df.append(fa_tmp_df, sort=False)

            else:
                sum_df = sum_df.append(tmp_df, sort=False)

            del tmp_df

        sum_df.sort_values(by=['Lib_mz', 'MS1_obs_mz', 'RANK_SCORE', 'ISOTOPE_SCORE', 'lipid_class'], inplace=True)
        sum_df.reset_index(drop=True, inplace=True)
        print(sum_df.head())

        return sum_df

    def merge_features(self, sum_df, merge_table, rank_score=40, isotope_score=80):

        unique_df = self.unique_features(sum_df, rank_score=rank_score, isotope_score=isotope_score)
        unique_df.to_excel(merge_table)

    def unique_features(self, sum_df, rank_score=40, isotope_score=80):

        unique_df = pd.DataFrame()

        unique_discrete_lst = sum_df['DISCRETE_ABBR'].unique().tolist()

        unique_charge_lst = sum_df['Charge'].unique().tolist()
        self.logger.info(unique_charge_lst)
        if len(unique_charge_lst) == 3:
            if '[M+NH4]+' in unique_charge_lst or '[M+Na]+' in unique_charge_lst or '[M+H]+' in unique_charge_lst:
                unique_charge_lst = ['[M+NH4]+', '[M+Na]+', '[M+H]+']
        self.headers = self.headers + unique_charge_lst
        if len(unique_charge_lst) > 1:
            for charge in unique_charge_lst:
                sum_df[charge] = ''

        discrete_count = len(unique_discrete_lst)

        self.logger.info(f'Total Discrete forms: {discrete_count}')

        d_i = 1
        for discrete in unique_discrete_lst:
            tmp_df = sum_df[sum_df['DISCRETE_ABBR'] == discrete]

            tmp_rank_score = rank_score
            if 'DG' in discrete:
                print(discrete, 'O-/P- lipids, lower TH')
                tmp_rank_score = 50
            if '-' in discrete:
                print(discrete, 'O-/P- lipids, lower TH')
                tmp_rank_score -= 10

            tmp1_df = tmp_df.query(f'RANK_SCORE >= {tmp_rank_score} and ISOTOPE_SCORE >= {isotope_score}')
            del tmp_df
            if discrete.startswith('P') or discrete.startswith('LP'):
                try:
                    if self.hg_check:
                        tmp1_df = tmp1_df[tmp1_df['#Specific_peaks'] > 0]
                except KeyError:
                    self.logger.warning('NO Specific_peaks column')
            if discrete.startswith('LP') or discrete.startswith('M'):
                try:
                    if self.fa_check:
                        tmp1_df = tmp1_df[tmp1_df['#Observed_FA'] >= 1]
                except KeyError:
                    self.logger.warning('NO Observed_FA column')
            if discrete.startswith('P') or discrete.startswith('D'):
                try:
                    if self.fa_check:
                        if '-' in discrete:
                            self.logger.debug(f'{discrete} contains O-/P- lipids, Check 1 FA only')
                            tmp1_df = tmp1_df[tmp1_df['#Observed_FA'] >= 1]
                        else:
                            tmp1_df = tmp1_df[tmp1_df['#Observed_FA'] >= 2]
                except KeyError:
                    self.logger.warning('NO Observed_FA column')
            elif discrete.startswith('T'):
                try:
                    if self.fa_check:
                        if '-' in discrete:
                            self.logger.debug(f'{discrete} contains O-/P- lipids, Check 1 FA only')
                            tmp1_df = tmp1_df[tmp1_df['#Observed_FA'] >= 2]
                        else:
                            tmp1_df = tmp1_df[tmp1_df['#Observed_FA'] >= 3]
                except KeyError:
                    self.logger.warning('NO Observed_FA column')
            tmp2_df = tmp1_df.sort_values(by=['RANK_SCORE', 'ISOTOPE_SCORE', 'abs_ppm'],
                                          ascending=[False, False, True])
            # del tmp1_df

            r_df = tmp2_df.head(1)  # type: pd.DataFrame()
            if not r_df.empty:
                self.logger.debug(f'# {d_i} / {discrete_count}: {discrete}')
                r_df.at[:, 'AVG_RANK_SCORE'] = tmp2_df['RANK_SCORE'].mean()
                r_df.at[:, 'AVG_ISOTOPE_SCORE'] = tmp2_df['ISOTOPE_SCORE'].mean()
                r_df.at[:, 'AVG_ppm'] = tmp2_df['abs_ppm'].mean()
                r_df.at[:, 'AVG_RT'] = tmp2_df['MS2_scan_time'].mean()
                r_df.at[:, 'MIN_RT'] = tmp2_df['MS2_scan_time'].min()
                r_df.at[:, 'MAX_RT'] = tmp2_df['MS2_scan_time'].max()
                r_df.at[:, 'IDENT_FA'] = tmp2_df['#Observed_FA'].max()
                if len(unique_charge_lst) >= 1:
                    tmp_charge_lst = tmp2_df['Charge'].unique().tolist()
                    self.logger.debug(tmp_charge_lst)
                    found_chg_lst = unique_charge_lst.copy()
                    for chg in unique_charge_lst:
                        if chg in tmp_charge_lst:
                            r_df.at[:, chg] = 'T'
                        else:
                            found_chg_lst.remove(chg)
                    r_df.at[:, 'Charge'] = ','.join(found_chg_lst)
                else:
                    self.logger.warning(unique_charge_lst)
                    found_chg_lst = unique_charge_lst.copy()

                try:
                    r_df.at[:, 'IDENT_HG'] = tmp2_df['#Specific_peaks'].max()
                except KeyError:
                    r_df.at[:, 'IDENT_HG'] = 0
                i_count = 0
                for file in self.file_abbr_lst:
                    try:
                        if tmp2_df[file].max() > 0:
                            r_df.at[:, file] = 1
                            i_count += 1
                    except KeyError:
                        r_df.at[:, file] = 0
                        #         if self.abbr_group_dct[file] in r_df.columns.tolist():
                #             r_df[self.abbr_group_dct[file]] += 1
                #         else:
                #             r_df.at[:, self.abbr_group_dct[file]] = 1
                for gp in self.group_abbr_dct:
                    r_df.at[:, gp] = r_df[self.group_abbr_dct[gp]].sum(axis=1)

                if discrete.startswith('TG') or discrete.startswith('DG'):
                    if '[M+NH4]+' in found_chg_lst:
                        unique_df = unique_df.append(r_df)
                        self.logger.debug(f'Collected {unique_df.shape[0]} / {d_i} From {discrete_count}')
                    else:
                        pass
                else:
                    unique_df = unique_df.append(r_df)
                    self.logger.debug(f'Collected {unique_df.shape[0]} / {d_i} From {discrete_count}')
                del r_df
                del tmp2_df
            d_i += 1

        unique_df.at[:, 'IDENT_COUNT'] = unique_df[list(self.group_abbr_dct.keys())].sum(axis=1)
        pre1_output_unique_df = unique_df[self.headers]
        del unique_df
        decimal_lst = {'AVG_RANK_SCORE': 2, 'AVG_ISOTOPE_SCORE': 2, 'AVG_ppm': 2, 'AVG_RT': 2, 'MIN_RT': 2, 'MAX_RT': 2}
        pre2_output_unique_df = pre1_output_unique_df.round(decimal_lst)
        output_unique_df = pre2_output_unique_df.sort_values(by=['lipid_class', 'Lib_mz', 'Proposed_structures'])

        output_unique_df.rename(columns=self.header_dct, inplace=True)
        m_header = output_unique_df.columns.tolist()
        rest_header_lst = [h for h in m_header if h not in self.header_lst]
        header_lst = self.header_lst + sorted(rest_header_lst)
        output_unique_df = output_unique_df[header_lst]

        for file in self.file_abbr_lst:
            output_unique_df[file] = output_unique_df[file].replace(range(1, 100), 'T')

        output_unique_df.reset_index(drop=True, inplace=True)

        return output_unique_df


if __name__ == '__main__':

    usr_header_info = r'../Configurations/Headers.xlsx'

    # cfg_folder = r'D:\AdipoAtlas_Project\HunterCfg'
    # usr_file_groups = r'../Configurations/file_groups_Polar_Neg.xlsx'
    # usr_output_path = r'output/LipidHunter_Polar_Neg_PL3.xlsx'

    # cfg_folder = r'D:\AdipoAtlas_Project\HunterCfg\IT_PL'
    # usr_file_groups = r'../Configurations/file_groups_Polar_Neg_IT.xlsx'
    # usr_output_path = r'output/LipidHunter_Polar_Neg_PL_IT3.xlsx'

    # cfg_folder = r'D:\AdipoAtlas_Project\HunterCfg\IT_PL'
    # usr_file_groups = r'Configurations/file_groups_Polar_Neg.xlsx'

    cfg_folder = r'D:\AdipoAtlas_Project\HunterCfg\IT_TGDG\DG'
    usr_file_groups = r'../Configurations/file_groups_Polar_Pos_IT.xlsx'
    usr_output_path = r'output/GL/LipidHunter_DG_polar_IT3.xlsx'

    # cfg_folder = r'D:\AdipoAtlas_Project\HunterCfg\AquireX_TGDG'
    # usr_file_groups = r'../Configurations/file_groups_Unpolar_Pos_AquireX.xlsx'
    # usr_output_path = r'output/GL/LipidHunter_TGDG_AquireX3.xlsx'

    # cfg_folder = r'D:\AdipoAtlas_Project\HunterCfg\IT_TGDG\TG_NH4'
    # usr_file_groups = r'../Configurations/file_groups_Polar_Pos_IT.xlsx'
    # usr_output_path = r'output/GL/LipidHunter_TGDG_polar_IT3.xlsx'

    # cfg_folder = r'D:\AdipoAtlas_Project\HunterCfg\TGDG'
    # usr_file_groups = r'../Configurations/file_groups_Polar_Pos.xlsx'
    # usr_output_path = r'output/GL/LipidHunter_TGDG_polar3.xlsx'

    # cfg_folder = r'D:\AdipoAtlas_Project\HunterCfg\TGDG_unpolar'
    # usr_file_groups = r'../Configurations/file_groups_Unpolar_Pos.xlsx'
    # usr_output_path = r'output/GL/LipidHunter_TGDG_unpolar3.xlsx'

    usr_cfg_lst = [os.path.join(cfg_folder, f) for f in os.listdir(cfg_folder)
                   if os.path.isfile(os.path.join(cfg_folder, f))]

    hunter_cfg = HunterConfig(usr_file_groups, usr_header_info)
    hunter_cfg.load_batch_cfg(usr_cfg_lst, merge_table=usr_output_path, rank_score=60, isotope_score=85)

    print('FIN')
