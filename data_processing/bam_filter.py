# -*- coding: UTF-8 -*-

import pysam
import multiprocessing as mp
import pandas as pd
import numpy as np

chr_length = {'1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276, '5':180915260, '6':171115067, '7':159138663, '8':146364022, '9':141213431, '10':135534747, '11':135006516, '12':133851895,  '13':115169878, '14':107349540, '15':102531392, '16':90354753, '17':81195210, '18':78077248, '19':59128983, '20':63025520, '21':48129895, '22':51304566, 'X':155270560, 'Y':59373566}

if __name__ == '__main__':

    import argparse, os

    parser = argparse.ArgumentParser()
    parser.add_argument('files', help="bam files", nargs = '+')
    parser.add_argument('-o', '--out', help="output file", default=None)
    parser.add_argument('-w', '--window', help="window", default=100000)
    parser.add_argument('-m', '--mapping_quality', help = "mapping quality threshold", default= 15)
    args = parser.parse_args()

    for f in args.files:

        bam_in = pysam.AlignmentFile(f, "rb")

        if args.out is not None:
            bam_out = pysam.AlignmentFile(args.out, "wb", template=bam_in)
        else:
            bam_out = pysam.AlignmentFile(os.path.splitext(f)[0] + '_filtered.bam', "wb", template=bam_in)


##간단한 filtering 조건을 통과하면서 다른 chromosome에 align 되지 않는 read들은 바로 output bamfile에 넣어줌.

        for chr in [str(i) for i in range(1,23)] + ['X', 'Y']:
        #chro = "22"
        for start in list(range(0, chr_length[chro], args.window)): #fetch 별 start

            print(chro+"_"+str(start))

            col_names = ['read_reference_name', 'read_reference_start', 'mate_reference_name', 'mate_reference_start', 'mate_mapping_quality', 'read_query_name']    
            mates_df = pd.DataFrame(columns = col_names)        
            
## 1차 filtering options (mapping quality/mate mapping/read pair)            
            for read in bam_in.fetch(chro, start, start+args.window): 
                if read.mapping_quality < args.mapping_quality: continue
                elif read.mate_is_unmapped: continue
                elif not read.is_paired: continue 
                elif read.mate_is_unmapped: continue    
                elif read.next_reference_name != read.reference_name: #mate가 다른 chromosome에 위치
                    
                    tmp_mates_list = [read.reference_name, read.reference_start, read.next_reference_name, read.next_reference_start, read.mapping_quality, read.query_name]
                    tmp_mates_dict = dict( zip(mates_df.columns, tmp_mates_list))
                    mates_df = mates_df.append(tmp_mates_dict, ignore_index=[True])

                    #print(mates_df)
                else: 
                    bam_out.write(read)

##2차 Filtering (discordant reads 중에서 mate가 같은 chromosome에 있는 것끼리 grouping-> 내림차순 정렬 -> mate의 position이 앞뒤로 1000미만 차이나는 read만 sorting)
            if mates_df.empty == True: continue        
            else:
                top = lambda data,col:data.sort_values(by=col,ascending=[False])
                mates_df = mates_df.groupby("mate_reference_name").apply(top,"mate_reference_start")            
                mates_df = mates_df.assign(diff_f=mates_df[["mate_reference_start"]].groupby("mate_reference_name").diff(periods=-1))
                mates_df = mates_df.reset_index(drop=[True])
                mates_df['diff_b'] = pd.Series([np.nan] + list(mates_df.diff_f[:-1]))
                mates_df['min_diffs'] = mates_df[["diff_f","diff_b"]].min(axis=1)
                mates_df = mates_df[mates_df.min_diffs < 1000]

                #print(mates_df.head(20))
                #print(len(mates_df))

##만들어 놓은 bam_in에서 filtering 된 read들을 뽑아(read_query_name으로 찾음) bam_out에 저장
            for r in bam_in.fetch(chro, start, start+args.window):
                if r.query_name in list(mates_df.read_query_name):
                    bam_out.write(r)
                else: continue        