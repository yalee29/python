# -*- coding: utf-8 -*-

'''

    Caller들 결과를 합쳐주는 tool

'''

import pandas as pd
import numpy  as np

def merge_call(list_df, list_caller):

    if len(list_df) != len(list_caller):
        print("The lengths of list_df and list_caller must be same!")
        return None

    df = list_df[0]
    df['caller'] = list_caller[0]
    df = df.set_index(['CHROM', 'POS', 'REF', 'ALT']).add_suffix('_' + list_caller[0])

    for i in range(len(list_df)-1):

        list_df[i + 1]['caller'] = list_caller[i + 1]
        list_df[i + 1] = list_df[i + 1].set_index(['CHROM', 'POS', 'REF', 'ALT']).add_suffix('_' + list_caller[i + 1])

        df=df.merge(list_df[i + 1], left_index=True, right_index=True, how='outer')

    # column sorting
    df.sort_index(axis=1, inplace=True)

    # index를 다시 되돌림
    df = df.reset_index()

    # caller 정보 합침
    df['caller']=''
    for c in list_caller:

        df['caller'] = df['caller'].str.cat(df['caller_' + c].str.strip(), sep=' ', na_rep='').str.strip()
        df = df.drop('caller_' + c, axis='columns')


    # sorting
    df_1_22   = df[ df['CHROM'].isin([str(i) for i in range(1, 23)])].copy()
    df_X_Y    = df[ df['CHROM'].isin(['X', 'Y'])].copy()
    df_others = df[~df['CHROM'].isin([str(i) for i in range(1, 23)] + ['X', 'Y'])].copy()

    df_1_22['CHROM'] = df_1_22['CHROM'].astype(int)

    df_1_22   = df_1_22.sort_values(by=['CHROM', 'POS'])
    df_X_Y    = df_X_Y.sort_values(by=['CHROM', 'POS'])
    df_others = df_others.sort_values(by=['CHROM', 'POS'])

    df_1_22['CHROM'] = df_1_22['CHROM'].astype(object)

    df = pd.concat([df_1_22, df_X_Y, df_others])


    return df


def reformat_varscan_alt(df):

    '''
    
        VarScan output의 alt field가 보통의 vcf format과 다르므로 
        vcf format에 맞도록 변경해주는 함수
     
    '''

    df['ALT'] = np.where(df['ALT'].str[0] == '+', df['REF'].str.cat(df['ALT'].str[1:]), df['ALT']) # insertion
    df['REF'] = np.where(df['ALT'].str[0] == '-', df['REF'].str.cat(df['ALT'].str[1:]), df['REF']) # deletion
    df['ALT'] = np.where(df['ALT'].str[0] == '-', df['REF'].str[0], df['ALT'])

    return df


if __name__ == '__main__':

    import argparse, vcf2pd, os, sys

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files",  help="합칠 file list",            nargs='+', required=True)
    parser.add_argument("-c", "--caller", help="각 file들을 만든 caller 정보", nargs='+', required=True)
    parser.add_argument("-o", "--out",    help="Output file 이름",          required=True)
    args = parser.parse_args()

    list_df = []
    for f in args.files:

        ext = f.strip().split('.')[-1]

        if ext == 'vs2': # varscan2 output 이면..

            # merge가 되도록 하기 위해 key의 형태를 맞춰준다.
            df = pd.read_csv(f, sep='\t', low_memory=False).rename(columns={'chrom': 'CHROM', 'position': 'POS', 'ref': 'REF', 'var': 'ALT'})
            df = reformat_varscan_alt(df)

            list_df.append(df)

        elif ext == 'vcf':
            list_df.append(vcf2pd.vcf2pd(f))

        elif ext == 'tsv':
            list_df.append(pd.read_csv(f, sep='\t', low_memory=False).rename(columns={'START':'POS'}))

        else:
            print(f, "는 지원되지 않는 파일입니다")
            sys.exit(1)


    df = merge_call(list_df, args.caller)

    ext = os.path.splitext(args.out)[1]
    if ext == '.csv':
        sep=','
    elif ext == '.tsv':
        sep='\t'

    df.to_csv(args.out, sep=sep, index=False)
