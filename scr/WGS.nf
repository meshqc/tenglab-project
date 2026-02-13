#!/bin/env nextflow
nextflow.enable.dsl=1

/*
 * Releases
 */

// user: heweihuang
// version: v1.0.1  whhe(2022-06-02) This workflow used to analyze WGS early screening programs
// version: v1.0.2  whhe(2022-09-15) mode advance and merge QC result 
// version: v1.0.3  whhe(2023-07-20) add marker CNV
/*
 * print usage
 */

params.h = false
params.help = false
if(params.help || params.h){
log.info ''
log.info 'WGS.nf'
log.info '=================================================================================================================================='
log.info 'methlylation sample for rawdata fastq.gz'
log.info 'Usage:'
log.info '  nextflow run BS_panel.nf --i /rawdata input path [*fq.gz]/ --t cpu_number --o /outut path/ '
log.info ''
log.info 'Options:'
log.info '      --help                      Show this message and exit.'
log.info '      --i                 <str>   input (fastq.gz/fq.gz) path                              [./]'
log.info '      --id                <str>   the sample name                                          [out]'
log.info '      --o                 <str>   the output dir                                           [./]'
log.info '  QC_options:'
log.info '      --qc                <str>   Whether the rawdata is qc filtered                       [true]'
log.info '      --qc_cpu            <int>   cpu number for each QC job                               [4]'
log.info '      --clean_rate        <int>   threshold of clean data rate                             [80]'
log.info '  Map_options:'
log.info '      --m                 <str>   Whether the reads is mapping the reference database      [true]'
log.info '      --map_cpu           <int>   cpu number for each rm rRNA job                          [8]'
log.info '      --mapping_rate      <int>   threshold of rRNA mapping rate                           [80]'
exit 1
}

outdirAbs              = ""

if( params.d[0] == "." ){ 
    outdirAbs = "${launchDir}/${params.d}"
}else{
    outdirAbs = "${params.d}"
}



Channel
    .fromFilePairs("$params.i/*{_1,_2}*{fq,fastq}*",size:2)
    .ifEmpty { error "Cannot find fq/fq.gz/fastq/fastq.gz files: $params.i" }
    .set{raw_reads}

process QC{
    tag "${libName}"
    cpus "${params.qc_cpu}"
   
    input:
    set val(libNameRaw), file(r_reads) from raw_reads
    output:
    set val(libName), file("*Clean.fastq.gz") into clean_fq1, clean_fq2
    set val(libName), file("*filter_ratio.xls") into fastp_clean_stat, fastp_clean_stat2
    set val(libName), file("*filter_ratio.new.xls") into qc_result 
    script:
       libName = libNameRaw.split(/[_]/)[0]

    """
    ${params.fastp} --thread ${params.qc_cpu} -i ${r_reads[0]} -I ${r_reads[1]} -o ${libName}.read1_Clean.fastq.gz -O ${libName}.read2_Clean.fastq.gz -l ${params.min_len} --detect_adapter_for_pe -w ${params.qc_cpu} --html=${libName}.fastp.html 1>fastp.log 2>fastp.err
    /public/home/weihuang//miniconda3/bin/python3 ${params.fastp_stat} fastp.err ./ ${libName}
    cat ${libName}.filter_ratio.xls|sed 's/,//g'|awk '{if(NR==1){print \$0}if(NR==2){split(\$5,a,"(");split(\$6,b,"(");split(\$1,c,"_");print c[1]"\\t"\$2*2"\\t"\$3*2"\\t"\$4*2"\\t"a[1]*2"("a[2]"\\t"b[1]*2"("b[2]"\\t"\$7*2"\\t"\$8}}' > ${libName}.filter_ratio.new.xls
    md5sum *Clean.fastq.gz > ${libName}.md5.txt
    """ 
}   

process check_QC {
    tag "${libName}"
    cpus "1"
    errorStrategy='finish'

    input:
    set val(libName), file(qc_s) from fastp_clean_stat
    """
    perl ${params.check_cleanR} -i ${qc_s} -threshold ${params.clean_rate}
    """
}

process Map{
    tag "${libName}"
    cpus "${params.mapping_cpu}"
    input:
    set val(libName), file(c_reads2) from clean_fq2
    output:
    set val(libName), file("*sort.mark.bam"), file("*sort.mark.bam.bai") into map_bam1, map_bam2, map_bam3
    set val(libName), file("*.map.xls") into map_result
    """
    export LD_LIBRARY_PATH="/public/home/weihuang/miniconda3/lib/:$LD_LIBRARY_PATH"
    ${params.bwa}  mem -C -R "@RG\\tID:${libName}\\tSM:${libName}\\tPL:Illumina"  -M  -Y -t ${params.mapping_cpu} ${params.fa} ${c_reads2[0]} ${c_reads2[1]} | ${params.sambamba} view -t ${params.mapping_cpu} -f bam -S -h /dev/stdin | ${params.sambamba} sort -m 10G -t ${params.mapping_cpu} --tmpdir=./ -o ${libName}.sort.bam /dev/stdin
    ${params.sinotools} mask_dup_parallel -i ${libName}.sort.bam -o ${libName} -l ${libName}.log -no_clip
    ${params.sinotools} un_map -i ${libName}.sort.bam -o ${libName}_chr0.bam
    ARGV=`ls ${libName}_chr*.bam|sort -V|perl -ne 'chomp;print "\$_ "'`
    ${params.sambamba} merge -t ${params.mapping_cpu} ${libName}.sort.mark.bam \$ARGV
    ${params.samtools} index ${libName}.sort.mark.bam
    ${params.samtools} flagstat ${libName}.sort.bam > ${libName}.flag.txt
    awk -vs=${libName} 'BEGIN{print s"\\tmap-rate"}{if(NR==9){split(\$0,a,/\\(|%/); print s"\\t"a[2]"%"}}' ${libName}.flag.txt > ${libName}.map.xls
    """
}

process Motif{
    tag "${libName}"
    cpus "${params.motif_cpu}"
    input:
    set val(libName), file(bam), file(bai) from map_bam1
    output:
    set val(libName), file("*.motif_original.csv") into out_motif
    set val(libName), file("*.fragment_original.csv") into out_frag
    """
    /public/home/weihuang/miniconda3/bin/python ${params.motif_frag} ${bam} ${libName} ${params.GenRange}
    /public/home/weihuang/miniconda3/bin/python ${params.motif_frag_6bp} ${bam} ${libName}_6bp ${params.GenRange} 
    /public/home/weihuang/miniconda3/bin/python ${params.fea_vis} ${libName}_motif_original.csv ./
    /public/home/weihuang/miniconda3/bin/python ${params.fea_vis} ${libName}_fragment_original.csv ./
    ln -s ${libName}_motif_original.csv ${libName}.motif_original.csv
    ln -s ${libName}_fragment_original.csv ${libName}.fragment_original.csv 
    """
}


workflow.onComplete {
    if( workflow.success ) {
        File f = new File("./pipe.Done")
        f.write("Pipeline completed at: $workflow.complete"+"\n") 
        f.append("Execution status:      ${ workflow.success ? 'OK' : 'failed' }" + "\n")
        f.append("User:                  $workflow.userName" + "\n")
        f.append("Launch time:           ${workflow.start.format('yyyy-MMM-dd HH:mm:ss')}" + "\n")
        f.append("Ending time:           ${workflow.complete.format('yyyy-MMM-dd HH:mm:ss')}" + "\n")
        f.append("Duration:              $workflow.duration" + "\n")
        f.append("Total CPU-Hours:       ${workflow.stats.computeTimeFmt ?: '-'}" + "\n")
        f.append("Tasks stats:           Succeeded ${workflow.stats.succeedCountFmt}; Cached ${workflow.stats.cachedCountFmt}; Ignored ${workflow.stats.ignoredCountFmt}; Failed ${workflow.stats.failedCountFmt}"  + "\n")
    }else{
        File f = new File("./pipe.Failed")
        f.write("Pipeline failed at: $workflow.complete"+"\n")
        f.append("Execution status:   ${ workflow.success ? 'OK' : 'failed' }" + "\n")
        f.append("User:               $workflow.userName" + "\n")
        f.append("Launch time:        ${workflow.start.format('yyyy-MMM-dd HH:mm:ss')}" + "\n")
        f.append("Ending time:        ${workflow.complete.format('yyyy-MMM-dd HH:mm:ss')}" + "\n")
        f.append("Duration:           $workflow.duration" + "\n")
        f.append("Total CPU-Hours:    ${workflow.stats.computeTimeFmt ?: '-'}" + "\n")
        f.append("Tasks stats:        Succeeded ${workflow.stats.succeedCountFmt}; Cached ${workflow.stats.cachedCountFmt}; Ignored ${workflow.stats.ignoredCountFmt}; Failed ${workflow.stats.failedCountFmt}" + "\n")
        f.append("ERROR message:\n" + "  ${workflow.errorMessage}" + "\n")
        f.append("ERROR report:\n" + "  ${workflow.errorReport}" + "\n")
    }
}
