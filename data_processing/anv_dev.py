# -*- coding: utf-8 -*-

'''

    Annovar

'''

import pandas as pd
import numpy  as np


if __name__ == '__main__':

    import argparse, os, subprocess, sys

    parser = argparse.ArgumentParser()
    parser.add_argument('file', help="annotation할 file (CHROM, POS, REF, ALT column이 꼭 있어야 함)", nargs='+')
    parser.add_argument('-p', '--protocol', help="annovar protocol",
                        default='refGene,cytoBand,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,exac03nontcga,esp6500siv2_all,esp6500siv2_aa,esp6500siv2_ea,avsnp150,clinvar_20180603,dbnsfp35c,dbscsnv11,intervar_20180118,cosmic87_coding,cosmic87_noncoding,dgvMerged,targetScanS,oreganno,tfbsConsSites,vistaEnhancers,switchDbTss,laminB1,cpgIslandExt,phastConsElements100way,phastConsElements46way,tRNAs,wgRna,rmsk')
    parser.add_argument('-b', '--buildver', help="Build version", default='hg19')
    parser.add_argument('-o', '--out',    help="Output file", default=None)
    args = parser.parse_args()

    operation = ''
    for o in args.protocol.split(','):
        o = o.strip()
        if o in ['refGene', 'ensGene']:
            operation += 'g,'
        elif o in ['cytoBand', 'dgvMerged', 'phastConsElements46way', 'phastConsElements100way', 'targetScanS', 'oreganno', 'vistaEnhancers', 'switchDbTss', 'laminB1', 'cpgIslandExt', 'tRNAs', 'wgRna', 'rmsk', 'tfbsConsSites']:
            operation += 'r,'
        else:
            operation += 'f,'

    operation = operation.rstrip(',')

    for f in args.file:

        print("\033[0;31mProcessing:\033[0m", f)

        df = pd.read_csv(f, sep='\t', low_memory=False)

        # backup columns
        df['POS_backup'] = df['POS']
        df['ALT_backup'] = df['ALT']
        df['REF_backup'] = df['REF']

        # annovar output과 합치기 위해서는 CHROM, POS, REF, ALT 형식을 바꿔 주어야 함
        df = df.rename(columns={'POS': 'START'})

        # multiallelic site일 경우 가장 앞쪽의 allele만 취한다 (convert2annovar.pl 도 같은 방식으로 작동함)
        df['ALT'] = list(map(lambda x: x[0], df['ALT'].str.split(',')))
        df['ALT'] = df['ALT'].str.strip()

        # snv
        snv = df['REF'].str.len() == df['ALT'].str.len()
        df['END'] = np.where(snv, df['START'] + df['REF'].str.len() - 1, np.nan)

        # insertion (start 위치 = insertion 직전 위치)
        insertion = ( df['REF'].str.len() < df['ALT'].str.len() ) & ~snv # snv locus와 중복수정되는 것을 막기 위함
        df['END'] = np.where(insertion, df['START'], df['END'])
        df['REF'] = np.where(insertion, '-', df['REF'])
        df['ALT'] = np.where(insertion, df['ALT'].str[1:], df['ALT'])

        # deletion (start 위치 = deletion 시작 위치)
        deletion = ( df['REF'].str.len() > df['ALT'].str.len() ) & ~insertion & ~snv # snv, insertion locus와 중복수정 되는 것믈 막기 위
        df['END']   = np.where(deletion, df['START'] + df['REF'].str.len() - 1, df['END'])
        df['START'] = np.where(deletion, df['START'] + 1, df['START'])
        df['REF']   = np.where(deletion, df['REF'].str[1:], df['REF'])
        df['ALT']   = np.where(deletion, '-', df['ALT'])

        if df['END'].isna().any():
            print("Some variants could not be interpreted!")
            print(df[df['END'].isna()])
            sys.exit(1)

        df['END'] = df['END'].astype(int) # float 형태로 저장되는 것을 막기 위해

        df['CHROM'] = np.where(df['CHROM'] == 'chrM', "MT", df['CHROM']) # chrM --> MT로 이름을 변경해주어야함

        # avinput 작성
        if args.out is not None:
            df[['CHROM', 'START', 'END', 'REF', 'ALT']].to_csv(os.path.splitext(args.out)[0] + '.avinput', sep='\t', index=False, header=False)
        else:
            df[['CHROM', 'START', 'END', 'REF', 'ALT']].to_csv(os.path.basename(f).split('.')[0] + '.avinput', sep='\t', index=False, header=False)

        # table_annovar.pl 실행
        if args.out is not None:
            subprocess.run(['table_annovar.pl',
                            os.path.splitext(args.out)[0] + '.avinput',
                            '/home/users/chrono0707/tools/annovar/humandb',
                            '-buildver', args.buildver,
                            '-out', os.path.splitext(args.out)[0],
                            '-remove',
                            '-protocol', args.protocol,
                            '-operation', operation,
                            '--dot2underline'
                            ])
        else:
            subprocess.run(['table_annovar.pl',
                            os.path.basename(f).split('.')[0] + '.avinput',
                            '/home/users/chrono0707/tools/annovar/humandb',
                            '-buildver', args.buildver,
                            '-out', os.path.basename(f).split('.')[0],
                            '-remove',
                            '-protocol', args.protocol,
                            '-operation', operation,
                            '--dot2underline'
                            ])

        # merge
        if args.out is not None:
            df_anv = pd.read_csv(os.path.splitext(args.out)[0] +  '.' + args.buildver + '_multianno.txt', sep='\t', low_memory=False)
        else:
            df_anv = pd.read_csv(os.path.basename(f).split('.')[0] + '.' + args.buildver + '_multianno.txt', sep='\t', low_memory=False)

        df_anv = df_anv.rename(columns={'Chr': 'CHROM', 'Start': 'START', 'End': 'END', 'Ref': 'REF', 'Alt': 'ALT'})
        df = df.merge(df_anv, on=["CHROM", "START", "END", "REF", "ALT"])

        # REF, ALT, START를 원래 형식으로 돌려놓음
        df = df.drop("START", axis=1)
        df = df.drop("END",   axis=1)
        df = df.drop("REF",   axis=1)
        df = df.drop("ALT",   axis=1)

        df = df.rename(columns={"POS_backup": "POS", "REF_backup": "REF", "ALT_backup": "ALT"})

        # CHROM, POS, REF, ALT를 맨 앞으로 보냄
        #print([x for x in df.columns.tolist() if x not in ['CHROM', 'POS', 'REF', 'ALT']])
        df = df[['CHROM', 'POS', 'REF', 'ALT'] + [x for x in df.columns.tolist() if x not in ['CHROM', 'POS', 'REF', 'ALT']]]

        if args.out is not None:
            df.to_csv(args.out, sep='\t', index=False)
        else:
            df.to_csv(os.path.basename(f).split('.')[0] + '_anv.tsv', sep='\t', index=False)
