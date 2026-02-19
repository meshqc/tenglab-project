#!/usr/bin/python

from collections import defaultdict
import pysam
import pandas as pd
#import HTSeq
import os,sys
import numpy as np

arg = pd.read_csv(sys.argv[3])
fn=sys.argv[1]
output_prefix=sys.argv[2]

#out=sys.argv[2]
print(fn.split('/')[-1][00:-4])
samA=pysam.AlignmentFile(fn,"rb")


read_seq=pysam.AlignedSegment()
filter_flag=3332
ans1 = []
tag0 = []
tag1 = []
tag2 = []
short_fragment_range = [100,150]
long_fragment_range = [151, 220]
length2count={}
total_reads =0
dic_motif = defaultdict(int)
for i in range(0,arg.shape[0]):
    short_dic = 0
    long_dic = 0
    #print("Processing {0}, start {1}, end {2}".format(arg.iat[i, 0],arg.iat[i, 1],arg.iat[i, 2]))
    for read_seq in samA.fetch(arg.iat[i, 0], arg.iat[i, 1], arg.iat[i, 2]):
        if (not (read_seq.flag & filter_flag) ) and read_seq.mapping_quality > 30:
            if read_seq.isize!=0 and read_seq.is_read1:
                total_reads +=1
                size = abs(read_seq.isize)
                if short_fragment_range[0] <= size <= short_fragment_range[1]:#short judge
                    short_dic += 1
                if long_fragment_range[0] <= size <= long_fragment_range[1]:#long judge
                    long_dic += 1
                if read_seq.is_read1:
                    if 'N' not in read_seq.query_sequence[0:4]:
                      dic_motif[read_seq.query_sequence[0:4]] += 1
                # this might cause problem
                if not read_seq.is_read1:
                    if 'N' not in read_seq.query_sequence[-4:]:
                      dic_motif[read_seq.query_sequence[-4:]] += 1
                #print("the Ratio is {}".format(short_dic / long_dic))
    ans1.append(short_dic / long_dic)
    tag0.append(arg.iat[i, 0])
    tag1.append(arg.iat[i, 1])
    tag2.append(arg.iat[i, 2])
add = sum(dic_motif.values())
temp = sorted(dic_motif.items(), key=lambda x:x[0],reverse = False)       
tag = []
nums = []
for (name, num) in temp:
    tag.append(name)
    nums.append(num / add)
dic_name = {"tag" : tag, "num" : nums}
df = pd.DataFrame(dic_name)

df.to_csv(output_prefix + '_motif_original' + '.csv', encoding='utf-8', index=False, sep='\t')


#print (length2count)
samA.close()

#print(samA)

#print(ans1)
dic_frag = {"chr":tag0, "start":tag1, "end":tag2, "Ratio":ans1}
df = pd.DataFrame(dic_frag)

df.to_csv(output_prefix + '_fragment_original' + '.csv', encoding='utf-8', index=False, sep='\t')


# motif normalization
nums_np = np.array(nums)
nums_afp = nums_np - nums_np.mean()
num_normalization = nums_afp.tolist()
dic_motif_afterprocess = {"tag":tag, "num":num_normalization}
df_after_process = pd.DataFrame(dic_motif_afterprocess)

df_after_process.to_csv(output_prefix + '_motif_normalize' + '.csv', encoding='utf-8', index=False, sep='\t')



#fragment normalization

ans1_np = np.array(ans1)
ans1_afp = ans1_np - ans1_np.mean()
ans1_normalization = ans1_afp.tolist()
dic_frag_after_process = {"chr": tag0, "start":tag1, "end":tag2, "Ratio":ans1_normalization}
df_after_process = pd.DataFrame(dic_frag_after_process)

df_after_process.to_csv(output_prefix + '_fragment_normalize' + '.csv', encoding='utf-8', index=False, sep='\t')

