# -*- coding: utf-8 -*-

'''

    2개의 file (.tsv, .VCF)를 하나로 합쳐주는 script
    합쳐지는 file은 동일한 header와 format을 가져야 한다.

    output은 tsv 형식임

    - VarScan2 output을 somatic, LOH, germline으로 분해해서 filtering하고
    다시 파일로 merging 하기 위해 필요함

    - snv.tsv. indel.tsv를 merge하는데 필요함

'''

import pandas as pd
import vcf2pd as vp
import os

def rbind(filenames, pos):

    if len(filenames) < 2:
        print("2개 이상의 file을 input으로 주어야 합니다.")
        return None

    ext = os.path.splitext(filenames[0])[1]

    if ext == '.vcf':

        merged_df = vp.vcf2pd(filenames[0])

        for f in filenames[1:]:

            df = vp.vcf2pd(f)

            # column이 다르면 에러
            if not (merged_df.columns == df.columns).all():
                print("파일들의 column이 같아야 합니다.")
                return None

            merged_df = pd.concat([merged_df, df])

            # sorting
            merged_df_1_22 = merged_df[merged_df['CHROM'].isin([str(i) for i in range(1, 23)])].copy()
            merged_df_X_Y = merged_df[merged_df['CHROM'].isin(['X', 'Y'])].copy()
            merged_df_others = merged_df[~merged_df['CHROM'].isin([str(i) for i in range(1, 23)] + ['X', 'Y'])].copy()

            merged_df_1_22['CHROM'] = merged_df_1_22['CHROM'].astype(int)

            merged_df_1_22 = merged_df_1_22.sort_values(by=['CHROM', 'POS'])
            merged_df_X_Y = merged_df_X_Y.sort_values(by=['CHROM', 'POS'])
            merged_df_others = merged_df_others.sort_values(by=['CHROM', 'POS'])

            merged_df_1_22['CHROM'] = merged_df_1_22['CHROM'].astype(object)

            merged_df = pd.concat([merged_df_1_22, merged_df_X_Y, merged_df_others])

    elif ext == '.tsv':

        merged_df = pd.read_csv(filenames[0], sep='\t', low_memory=False)

        for f in filenames[1:]:

            df = pd.read_csv(f, sep='\t', low_memory=False)
            df.columns = merged_df.columns

            # .tsv merge는 column이 완전 동일할 필요는 없지만 적어도 개수는 같아야 한다
            # ( snv.tsv와 indel.tsv의 sample 파일 이름이 다를 수 있기 때문 )
            if len(merged_df.columns) != len(df.columns):
                print("파일들의 column이 같아야 합니다.")
                return None

            merged_df = pd.concat([merged_df, df], sort=False)

            # sorting
            merged_df_1_22    = merged_df[merged_df['CHROM'].isin([ str(i) for i in range(1,23) ])].copy()
            merged_df_X_Y     = merged_df[merged_df['CHROM'].isin(['X', 'Y']) ].copy()
            merged_df_others  = merged_df[~merged_df['CHROM'].isin([ str(i) for i in range(1,23) ] + ['X','Y'])].copy()

            merged_df_1_22['CHROM'] = merged_df_1_22['CHROM'].astype(int)

            merged_df_1_22   = merged_df_1_22.sort_values(by=['CHROM', pos])
            merged_df_X_Y    = merged_df_X_Y.sort_values(by=['CHROM', pos])
            merged_df_others = merged_df_others.sort_values(by=['CHROM', pos])

            merged_df_1_22['CHROM'] = merged_df_1_22['CHROM'].astype(object)

            merged_df = pd.concat([merged_df_1_22, merged_df_X_Y, merged_df_others])

    return merged_df

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('files', help="합쳐질 파일들 (.tsv, .VCF)", nargs='+')
    parser.add_argument('-p', '--pos', help="POS 변수 이름", choices=['POS', 'START'], default='POS')
    parser.add_argument('-o', '--out', help="결과(tsv 파일)이 저장될 경로", required=True)
    args = parser.parse_args()

    df = rbind(args.files, args.pos)
    df.to_csv(args.out, sep='\t', index=False)

