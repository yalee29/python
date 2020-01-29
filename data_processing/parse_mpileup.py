# -*- coding: utf-8 -*-

'''

    mpileup parser

'''


import re, pandas as pd

def mpileup_parser(df_mp, qual_threshold):

    l = []

    for row in df_mp.itertuples():  # row의 첫번쩨 element는 index가 들어감

        CHROM = row[1]
        POS   = row[2]
        REF   = row[3]

        #depth = row[4]  # pileup에 저장된 depth는 deletion을 모두 포함하므로 진정한 depth라고 할 수 없음

        read  = row[5]
        qual  = row[6]


        # $, ^, * 제거
        read = re.sub('\^.', '', read)
        read = read.replace("$", "")

        # insertion/deletion 처리
        M = dict()
        while re.search(r'[\+\-][0-9]+', read):
            indel = re.search(r'[\+\-][0-9]+', read)
            type  = indel.group()[0]
            n     = int(indel.group()[1:len(indel.group())])
            seq   = type + read[indel.end():(indel.end() + n)]

            if seq in M.keys():  # 이미 있는 key면 count update
                M[seq] = M[seq] + 1
            else: # 처음 등장하는 key면 count를 1로 핳당
                M[seq] = 1

            read  = read[:indel.start()] + read[(indel.end() + n):]

        # filtering based on quality
        if len(read) != len(qual):
            print(len(read), read)
            print(len(qual), qual)
            print("Read and quality string are not in same length at {:s}:{:d}!!".format(CHROM, POS))
            return None

        read = ''.join([read[i] for i in range(len(read)) if ord(qual[i]) >= qual_threshold])

        forward_match = read.count(".")  # forward strand match
        reverse_match = read.count(",")  # reverse strand match
        forward_A     = read.count("A")  # forward strand A
        forward_G     = read.count("G")  # forward strand G
        forward_C     = read.count("C")  # forward strand C
        forward_T     = read.count("T")  # forward strand T
        forward_N     = read.count("N")  # forward strand N
        reverse_A     = read.count("a")  # reverse strand A
        reverse_G     = read.count("g")  # reverse strand G
        reverse_C     = read.count("c")  # reverse strand C
        reverse_T     = read.count("t")  # reverse strand T
        reverse_N     = read.count("n")  # reverse strand n
        #placeholder   = read.count("*") # place holder

        depth = forward_match + reverse_match + forward_A + forward_G + forward_T + forward_C + reverse_A + reverse_G + reverse_T + reverse_C + forward_N + reverse_N

        if REF == 'A':
            forward_A = forward_match
            reverse_A = reverse_match
        elif REF == 'G':
            forward_G = forward_match
            reverse_G = reverse_match
        elif REF == 'C':
            forward_C = forward_match
            reverse_C = reverse_match
        elif REF == 'T':
            forward_T = forward_match
            reverse_T = reverse_match

        l.append({'CHROM': CHROM, 'POS': POS, 'REF': REF, 'depth': depth,
                  'A': forward_A + reverse_A,
                  'G': forward_G + reverse_G,
                  'C': forward_C + reverse_C,
                  'T': forward_T + reverse_T,
                  'N': forward_N + reverse_N,
                  'Indel': M
                  })


    df = pd.DataFrame(l)

    return df[['CHROM', 'POS', 'REF', 'depth', 'A', 'G', 'C', 'T', 'N', 'Indel']]


if __name__ == '__main__':

    import argparse, os, sys

    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='mpileup')
    parser.add_argument('-q', '--quality', help="base quality threshold", default=13)
    parser.add_argument('-o', '--out', help='sample name list', default=None, nargs='*')
    args = parser.parse_args()

    ext = os.path.splitext(args.file)[1]

    if ext == '.gz':
        df = pd.read_csv(args.file, sep='\t', compression='gzip', header=None)
    else:
        df = pd.read_csv(args.file, sep='\t', header=None)

    # sample 수 계산
    n = int((df.shape[1] - 3) / 3)

    for i in range(n):
        print("processing sample {:d}".format(i))
        df_temp = mpileup_parser(df.iloc[:, [0, 1, 2, (i + 1) * 3, (i + 1) * 3 + 1, (i + 1) * 3 + 2]],
                                     qual_threshold=int(args.quality))

        if df_temp is None:
            print("Failed to parse!")
            sys.exit()

        if args.out is not None:
            if len(args.out) != n:
                print("sample list length does not match to sample number!!")
                sys.exit()

            out_name = args.out[i]
        else:
            out_name = 'sample_{:d}.mp'.format(i)

        df_temp.to_csv(out_name, sep='\t', index=False)

