//
// Run Fastp and FastQC on trimmed reads
// Based on nf-core/ViralRecon subworkflow: fastq_trim_fastp_fastqc.nf (https://github.com/nf-core/viralrecon/blob/master/subworkflows/local/fastq_trim_fastp_fastqc.nf)
//

include { FASTP                 } from '../../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_TRIM } from '../../modules/nf-core/fastqc/main'

//
// Function that parses fastp json output file to get total number of reads after trimming
//
import groovy.json.JsonSlurper

def getFastpReadsAfterFiltering(json_file) {
    def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
    return json['after_filtering']['total_reads'].toInteger()
}

workflow FASTP_AND_FASTQC_TRIM {
    take:
    reads             // channel: [ val(meta), [ reads ] ]
    adapter_fasta     //    file: adapters.fasta
    save_trimmed_fail //   value: boolean
    save_merged       //   value: boolean

    main:

    ch_versions = Channel.empty()

    trim_reads        = reads
    trim_json         = Channel.empty()
    trim_html         = Channel.empty()
    trim_log          = Channel.empty()
    trim_reads_fail   = Channel.empty()
    trim_reads_merged = Channel.empty()
    fastqc_trim_html  = Channel.empty()
    fastqc_trim_zip   = Channel.empty()
    FASTP (
        reads,
        adapter_fasta,
        save_trimmed_fail,
        save_merged
    )
    trim_reads        = FASTP.out.reads
    trim_json         = FASTP.out.json
    trim_html         = FASTP.out.html
    trim_log          = FASTP.out.log
    trim_reads_fail   = FASTP.out.reads_fail
    trim_reads_merged = FASTP.out.reads_merged
    ch_versions       = ch_versions.mix(FASTP.out.versions.first())


    //
    // Filter empty FASTQ files after adapter trimming so FastQC doesn't fail
    //
    trim_reads
        .join(trim_json)
        .map {
            meta, reads, json ->
                if (getFastpReadsAfterFiltering(json) > 0) {
                    [ meta, reads ]
                }
        }
        .set { trim_reads }

    FASTQC_TRIM (
        trim_reads
    )
    fastqc_trim_html = FASTQC_TRIM.out.html
    fastqc_trim_zip  = FASTQC_TRIM.out.zip
    ch_versions      = ch_versions.mix(FASTQC_TRIM.out.versions.first())

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]
    trim_json          // channel: [ val(meta), [ json ] ]
    trim_html          // channel: [ val(meta), [ html ] ]
    trim_log           // channel: [ val(meta), [ log ] ]
    trim_reads_fail    // channel: [ val(meta), [ fastq.gz ] ]
    trim_reads_merged  // channel: [ val(meta), [ fastq.gz ] ]

    fastqc_trim_html // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip  // channel: [ val(meta), [ zip ] ]

    versions = ch_versions.ifEmpty(null) // channel: [ version.yml ]
}
