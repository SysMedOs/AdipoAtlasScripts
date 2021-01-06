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

lh_file = r'LipidHunter_TG_lite.xlsx'
lh_df = pd.read_excel(lh_file)

ls_file = r'LipidSearch_TG_lite.xlsx'
ls_df = pd.read_excel(ls_file)

gl_rgx = re.compile(r'(?P<GL>\w{2,3})\((?P<FAs>.*)\)')

for idx, r in lh_df.iterrows():
    pre_h = r['Discrete']
    l_match = re.match(gl_rgx, pre_h)
    if l_match:
        h_info_dct = l_match.groupdict()
        h = f"{h_info_dct['GL']}({'_'.join(natsorted(h_info_dct['FAs'].split('_')))})"
        lh_df.at[idx, 'Discrete'] = h

m_df = lh_df.merge(ls_df, how='outer', on=['Class', 'Bulk', 'Discrete', 'Charge'])

print(m_df.shape)

m_df.to_excel('LH_LS_TG.xlsx')
print('FIN')
