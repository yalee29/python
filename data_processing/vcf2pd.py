'''
    vcf2를 parsing해서 pandas dataframe으로 만들어주는 함수
    pandas로 만들게 되면 excel, csv 등 다른 format으로 쉽게 변경이 가능.
'''

import argparse, os, re
import pandas as pd
import numpy  as np
from   cyvcf2 import VCF



def split_multiallelic_site(df_mt2):
    '''

        multiallelic site를 분리

    '''

    df_list = list()

    for _, row in df_mt2.iterrows():

        alt = row['ALT'].split(',')

        # multiallelic site라면
        if len(alt) > 1:

            for i, a in enumerate(alt):

                new_row = row.copy() # refenrece하지 않기 위해 dict를 copy

                new_row['ALT'] = a

                for key in new_row.keys():

                    # alt1, alt2, ... 구조
                    if len(str(new_row[key]).split(',')) == len(alt):

                        if new_row[key] is not None:
                            if isinstance(new_row[key], str):
                                new_row[key] = new_row[key].split(',')[i]
                            else:
                                new_row[key] = new_row[key][i]

                    # ref, alt1, alt2, ... 구조
                    elif len(str(new_row[key]).split(',')) == (len(alt) + 1):

                        if new_row[key] is not None:
                            if isinstance(new_row[key], str):
                                new_row[key] = new_row[key].split(',')[0] + ',' + new_row[key].split(',')[i+1]
                            else:
                                new_row[key] = str(new_row[key][0]) + ',' + str(new_row[key][i+1])

                df_list.append(new_row)

        else:

            df_list.append(row)

    df = pd.DataFrame(df_list)

    return df



def vcf2pd(vcf_in):
    '''
    VCF 파일을 읽어서 pandas dataframe 형식으로 return함
    :param vcf_in: VCF 파일 (.vcf/.vcf.gz/.bcf)
    :return: pandas dataframe
    '''

    vcf      = VCF(vcf_in, gts012=True)
    lRows    = []  # row의 list를 만들어서 마지막에 DataFrame으로 변환하는게 가장 빠르다.
    lInfo    = []  # INFO list
    lFormat  = []  # FORMAT list

    samples   = vcf.samples
    n_samples = len(samples)

    # INFO FIELD의 item을 얻는다
    for h in vcf.header_iter():
        if( h['HeaderType'] == 'INFO' ):
            lInfo.append( h['ID'] )

        if( h['HeaderType'] == 'FORMAT' ):
            lFormat.append( h['ID'])

    # sample information이 있는지 찾는다 (mutect2 output의 경우 이 정보가 포함되어있음)
    re_tumor  = re.compile('##tumor_sample=.*')
    re_normal = re.compile('##normal_sample=.*')
    if re_tumor.search(vcf.raw_header) is not None:
        samples[samples.index(re_tumor.search(vcf.raw_header).group().split('=')[1])]='TUMOR'

    if re_normal.search(vcf.raw_header) is not None:
        samples[samples.index(re_normal.search(vcf.raw_header).group().split('=')[1])] = 'NORMAL'

    for v in vcf:

        # 8개의 fixed field를 저장한다.

        dVariant = { 'CHROM': v.CHROM, 'POS': v.POS, 'ID': v.ID, 'REF': v.REF, 'ALT': ','.join(v.ALT),
                     'QUAL': v.QUAL, 'FILTER': v.FILTER}

        if not dVariant['FILTER']:     # cyvcf2에서는 FILTER가 PASS일때 FILTER=None으로 저장하기 때문에 다시 'PASS'로 돌려줌
            dVariant['FILTER'] = 'PASS'

        # INFO field 처리
        for i in lInfo:
            dVariant[i] = v.INFO.get(i)

        # FORMAT field 처리
        for f in lFormat:

            for i in range(n_samples):

                if f == 'GT':
                    dVariant[samples[i] + '_' + f] = v.gt_types[i]
                    # v.format('GT')에는 이상한 형식으로 저장됨.
                    # gt_type = 0 --> hom_ref, gt_type = 1 --> hetero, gt_type = 2 --> hom_alt, gt_type = 3 --> unknown
                else:

                    if v.format(f) is not None: # field가 None이 아니면

                        if isinstance(v.format(f)[i], str):
                            dVariant[samples[i] + '_' + f] = str(v.format(f)[i])  # string일 경우 각 letter들이 comma로 구분되는 것 방지
                        elif np.isnan(v.format(f)[i]).any(): # nan일 경우..
                            dVariant[samples[i] + '_' + f] = None
                        else:
                            dVariant[samples[i] + '_' + f] = ','.join(list(map(str, v.format(f)[i])))

        lRows.append(dVariant)

    cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
    cols.extend(lInfo)
    cols.extend([ s + '_' + f for s in samples for f in lFormat ])

    if lRows:
        df = pd.DataFrame(lRows, columns=cols)
    else:
        df = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER'])


    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help=".vcf/.vcf.gz/.bcf", nargs='+')
    parser.add_argument("-e", "--excel", action="store_true", help="Parsing된 VCF를 excel로 변환")
    parser.add_argument("-c", "--csv",   action="store_true", help="Parsing된 VCF를 csv로 변환")
    parser.add_argument("-t", "--tsv",   action="store_true", help="Parsing된 VCF를 tsv로 변환")
    parser.add_argument("-s", "--split", action="store_true", help="Multiallelic site 분리")
    args = parser.parse_args()


    for f in args.vcf:

        df_vcf = vcf2pd(f)
        if args.split:
            df_vcf = split_multiallelic_site(df_vcf)

        if args.excel:
            df_vcf.to_excel(os.path.basename(f).split('.')[0] + '.xlsx', index=False)

        if args.csv:
            df_vcf.to_csv(os.path.basename(f).split('.')[0] + '.csv', na_rep='', index=False)

        if args.tsv:
            df_vcf.to_csv(os.path.basename(f).split('.')[0] + '.tsv', na_rep='', sep='\t', index=False)
