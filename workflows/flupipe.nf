/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowFlupipe.initialise(params, log)


// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.fasta,
    params.adapter_fasta,
    params.kraken2_db,
    params.multiqc_config
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.adapter_fasta) { ch_adapter_fasta = file(params.adapter_fasta) } else { exit 1, 'Input adapter_fasta not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK           } from '../subworkflows/local/input_check'
include { PREPARE_ENVIRONMENT   } from '../subworkflows/local/prepare_environment'
include { FASTP_AND_FASTQC_TRIM } from '../subworkflows/local/fastp_and_fastqc_trim'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_RAW        } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2  } from '../modules/nf-core/kraken2/kraken2/main'
include { UNICYCLER                   } from '../modules/nf-core/unicycler/main'
include { BLAST_BLASTN                } from '../modules/nf-core/blast/blastn/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow FLUPIPE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Prepare all necessary files to run the pipeline
    //
    PREPARE_ENVIRONMENT ()
    ch_versions = ch_versions.mix(PREPARE_ENVIRONMENT.out.versions)

    //
    // MODULE: Run FastQC on raw reads
    //
    if (!params.skip_fastqc) {
        FASTQC_RAW (
            INPUT_CHECK.out.reads
        )
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
    }

    //
    // SUBWORKFLOW: Run Fastp and FastQC on trimmed reads
    //
    FASTP_AND_FASTQC_TRIM (
        INPUT_CHECK.out.reads,
        ch_adapter_fasta,
        params.fastp_save_trimmed_fail,
        params.fastp_save_merged
    )
    ch_trimmed_reads = FASTP_AND_FASTQC_TRIM.out.reads
    ch_versions = ch_versions.mix(FASTP_AND_FASTQC_TRIM.out.versions.first())

    //
    // Filter empty FASTQ files after adapter trimming
    // NOTE: This code is part of nf-core/ViralRecon pipeline.
    //
    ch_fail_reads_multiqc = Channel.empty()
    if (!params.skip_fastp) {
        ch_trimmed_reads
            .join(FASTP_AND_FASTQC_TRIM.out.trim_json)
            .map {
                meta, reads, json ->
                    pass = WorkflowFlupipe.getFastpReadsAfterFiltering(json) > 0
                    [ meta, reads, json, pass ]
            }
            .set { ch_pass_fail_reads }

        ch_pass_fail_reads
            .map { meta, reads, json, pass -> if (pass) [ meta, reads ] }
            .set { ch_variants_fastq }

        ch_pass_fail_reads
            .map {
                meta, reads, json, pass ->
                if (!pass) {
                    fail_mapped_reads[meta.id] = 0
                    num_reads = WorkflowFlupipe.getFastpReadsBeforeFiltering(json)
                    return [ "$meta.id\t$num_reads" ]
                }
            }
            .collect()
            .map {
                tsv_data ->
                    def header = ['Sample', 'Reads before trimming']
                    WorkflowCommons.multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_reads_multiqc }
    }

    //
    // MODULE: Remove host reads using Kraken2 with humanDB
    //
    ch_kraken2_multiqc = Channel.empty()
    KRAKEN2 (
        ch_trimmed_reads,
        PREPARE_ENVIRONMENT.out.kraken2_db,
        params.kraken2_save_output_fastqs,
        params.kraken2_save_reads_assigment
    )
    ch_nonhuman_reads = KRAKEN2.out.unclassified_reads_fastq
    ch_kraken2_multiqc = KRAKEN2.out.report
    ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())

    //
    // MODULE: Preliminary assembly using Unicycler (SPAdes)
    //
    UNICYCLER (
        ch_nonhuman_reads.map { meta, fastq -> [ meta, fastq, [] ] }
    )
    ch_versions = ch_versions.mix(UNICYCLER.out.versions.first())

    //
    // MODULE: Detect hits of previous assembly step using BLASTn with NCBI Influenza Virus Database:
    //
    // TODO

    //
    // MODULE: Run DumpSoftwareVersions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowFlupipe.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowFlupipe.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP_AND_FASTQC_TRIM.out.trim_json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP_AND_FASTQC_TRIM.out.fastqc_trim_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_kraken2_multiqc.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
