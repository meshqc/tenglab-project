import sys
import pysam
from collections import defaultdict

from parser import parse_fusion_bam, Segment

from multiprocessing import Pool
import time
from multiprocessing.managers import BaseManager, DictProxy, ListProxy
from collections import defaultdict

class MyManager(BaseManager):
    pass

MyManager.register('defaultdict', defaultdict, DictProxy)

def read_pair_generator(chrom,bam, region_string=None):
    read_dict = defaultdict(lambda: [None, None])
    bam.reset()
    for read in bam.fetch(str(chrom),until_eof=True, region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def count(bam_count,contig,pos):
    #bam=pysam.AlignmentFile(BAM,'rb')
    maxpos=bam_count.get_reference_length(contig)
    linear=0
    circ=0
    same_strand=0
    diff_strand=0
    for read in bam_count.fetch(contig,max(pos-2,0),min(pos+2,maxpos)):
        mapped=[x[1] for x in read.aligned_pairs]
        if (pos-2 in mapped) & (pos+2 in mapped):
            linear+=1
        else:
            circ+=1
        
    for reads in bam_count.fetch(contig,max(0,pos-5000),min(pos+5000,maxpos)):
        if reads.is_reverse == reads.mate_is_reverse:
            same_strand+=1
        
        else:
            diff_strand+=1
    return linear,circ,same_strand,diff_strand
    


def countReads(chrom,BAM,    fusions,mate_reference_length,mate_mapping_quality,l_mapped,r_mapped,ll,lc,ls,ld,rl,rc,rs,rd
):
    temp_fusions=defaultdict(int)
    temp_mate_reference_length=defaultdict(int)
    temp_mate_mapping_quality=defaultdict(int)
    temp_l_mapped=defaultdict(int)
    temp_r_mapped=defaultdict(int)
    samFile_countReads = pysam.AlignmentFile(BAM,'rb')
    for read1_aligns, read2_aligns in read_pair_generator(chrom,samFile_countReads):
        if read1_aligns.is_unmapped or read2_aligns.is_unmapped:  # unmapped reads
            continue
        if read1_aligns.is_supplementary or read2_aligns.is_supplementary:  # supplementary reads
            continue
        if not read1_aligns.has_tag('SA'):#  no SA tag
            if not read2_aligns.has_tag('SA'):
                continue
            else:
                read=read2_aligns
                mate=read1_aligns
        else:
            if not read2_aligns.has_tag('SA'):
                read=read1_aligns
                mate=read2_aligns
            else:
                continue
        if read.reference_id != mate.reference_id: # mate not same chromosome
            continue
        chr1 = samFile_countReads.get_reference_name(read.reference_id)
        strand1 = '+' if not read.is_reverse else '-'
        saInfo = read.get_tag('SA').split(';')[:-1]
        loc = [read.query_alignment_start, read.query_alignment_end, read.reference_start, read.reference_end]
        segments = [loc]
        for sa in saInfo:
            chr2, pos, strand2, cigar = sa.split(',')[:4]
            if chr1 != chr2:  # not same chromosome
                continue
            if strand1 != strand2:  # not same strand
                continue
            segment = Segment(pos=pos, cigar=cigar)
            segments.append([segment.read_start, segment.read_end,
                             segment.ref_start, segment.ref_end])
        segments.sort()
        cov_loc = segments[0][1]
        ref_loc = segments[0][3]
        bflag = 0
        cigar_l, cigar_r = 0, 0
        ref_l, ref_r = 0, 0
        for s in segments[1:]:
            if s[2] < ref_loc and s[0] <= cov_loc:
                bflag += 1
                cigar_l = s[0]
                cigar_r = cov_loc
                ref_l = s[2]
                ref_r = ref_loc
            cov_loc = s[1]
            ref_loc = s[3]
        if bflag == 1:
            fusion_left = str(ref_l)
            fusion_right = str(ref_r - (cigar_r - cigar_l))
            fusion_loc = '\t'.join([chr1, fusion_left, fusion_right])
            temp_fusions[fusion_loc]+= 1
            temp_mate_reference_length[fusion_loc]=mate.reference_length
            temp_mate_mapping_quality[fusion_loc]=mate.mapping_quality
            temp_l_mapped[fusion_loc]=s[1]-s[0]
            temp_r_mapped[fusion_loc]=segments[0][1] - segments[0][0]
    for pos in temp_fusions.keys():
        fusions[pos]=temp_fusions[pos]
        mate_reference_length[pos]=temp_mate_reference_length[pos]
        mate_mapping_quality[pos]=temp_mate_mapping_quality[pos]
        l_mapped[pos]=temp_l_mapped[pos]
        r_mapped[pos]=temp_r_mapped[pos]
        ll[pos],lc[pos],ls[pos],ld[pos]=count(samFile_countReads,chr1,int(fusion_left))
        rl[pos],rc[pos],rs[pos],rd[pos]=count(samFile_countReads,chr1,int(fusion_right))
    print(chrom)
    #print(fusions)
    

if __name__ == "__main__":
    start = time.time()
    input_bam=sys.argv[1]
    cpu=int(sys.argv[2])
    outfile=sys.argv[3]
    samFile_main = pysam.AlignmentFile(input_bam, 'rb')
    chroms=samFile_main.references
    mgr = MyManager()
    mgr.start()
    fusions = mgr.defaultdict(int)
    mate_reference_length=mgr.defaultdict(int)
    mate_mapping_quality=mgr.defaultdict(int)
    l_mapped=mgr.defaultdict(int)
    r_mapped=mgr.defaultdict(int)
    ll=mgr.defaultdict(int)
    lc=mgr.defaultdict(int)
    ls=mgr.defaultdict(int)
    ld=mgr.defaultdict(int)
    rl=mgr.defaultdict(int)
    rc=mgr.defaultdict(int)
    rs=mgr.defaultdict(int)
    rd=mgr.defaultdict(int)



    total_reads=sum([x[3] for x in samFile_main.get_index_statistics()])
    
    pool = Pool(processes=cpu)
    for x in range(len(chroms)):
        pool.apply_async(countReads,(chroms[x],input_bam,fusions,mate_reference_length,mate_mapping_quality,l_mapped,r_mapped,
ll,lc,ls,ld,rl,rc,rs,rd))
    pool.close()
    pool.join()

        
        
    print(fusions)
    total = 0
    with open(outfile, 'w') as outF:
        outF.write("name\tchr\tstart\tend\tcount\tmate_mapped\tmate_mapping_quality\tl_mapped\tr_mapped\t")
        outF.write("left_linear_reads\tleft_circ_reads\tleft_mate_ss\tleft_mate_ds\t")
        outF.write("right_linear_reads\tright_circ_reads\tright_mate_ss\tright_mate_ds\n")
        i=0
        for pos in fusions.keys():
            i+=1
   
            outF.write('FUSIONJUNC_%d\t%s\t%d\t%d\t%d\t%d\t%d\t' % (i,pos, fusions[pos],mate_reference_length[pos],
                                                       mate_mapping_quality[pos],l_mapped[pos],r_mapped[pos]))
            outF.write('%d\t%d\t%d\t%d\t'%(ll[pos],lc[pos],ls[pos],ld[pos]))
            outF.write('%d\t%d\t%d\t%d\n'%(rl[pos],rc[pos],rs[pos],rd[pos]))
            total += fusions[pos]
    outF.close()
    samFile_main.close()
    print('\t Converted %d fusion reads from total %d reads from file %s!' % (total,total_reads,input_bam))
    end = time.time()
    print(end - start)
    
