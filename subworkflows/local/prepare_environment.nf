//
// Download, uncompress and prepare all necessary files
// Based on nf-core/ViralRecon subworkflow: prepare_genome.nf (https://github.com/nf-core/viralrecon/blob/master/subworkflows/local/prepare_genome.nf)
//

include { UNTAR as UNTAR_KRAKEN2_DB } from '../../modules/nf-core/untar/main'
include { GUNZIP as GUNZIP_FLU_DB   } from '../../modules/nf-core/gunzip/main'
include { KRAKEN2_BUILD             } from '../../modules/local/kraken2_build'
include { BLAST_MAKEBLASTDB         } from '../../modules/nf-core/blast/makeblastdb/main'

workflow PREPARE_ENVIRONMENT {
    main:

    ch_versions = Channel.empty()

    //
    // Prepare files to build humanDB for Kraken2
    //
    ch_kraken2_db = Channel.empty()
    if (params.kraken2_db) {
        if (params.kraken2_db.endsWith('.tar.gz')) {
            UNTAR_KRAKEN2_DB (
                [ [:], params.kraken2_db ]
            )
            ch_kraken2_db = UNTAR_KRAKEN2_DB.out.untar.map { it[1] }
            ch_versions = ch_versions.mix(UNTAR_KRAKEN2_DB.out.versions)
        } else {
            ch_kraken2_db = Channel.value(file(params.kraken2_db))
        }
    } else {
        KRAKEN2_BUILD (
            params.kraken2_db_name
        )
        ch_kraken2_db = KRAKEN2_BUILD.out.db
        ch_versions   = ch_versions.mix(KRAKEN2_BUILD.out.versions)
    }

    //
    // Prepare files to build NCBI Influenza Virus Database for BLASTn
    //
    ch_flu_fna = Channel.empty()
    if (params.flu_fna) {
        if (params.flu_fna.endsWith('.gz')) {
            GUNZIP_FLU_DB (
                [ [:], params.flu_fna ]
            )
            ch_flu_fna = GUNZIP_FLU_DB.out.gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_FLU_DB.out.versions)
        } else {
            ch_flu_fna = Channel.value(file(params.flu_fna))
        }

        PARSE_FLU_DB (
            ch_flu_fna
        )
        ch_versions = ch_versions.mix(PARSE_FLU_DB.out.versions)

        ch_flu_db = Channel.empty()
        BLAST_MAKEBLASTDB (
            PARSE_FLU_DB.out.fna
        )
        ch_flu_db = BLAST_MAKEBLASTDB.out.db
        ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)
    }

    emit:
    kraken2_db    = ch_kraken2_db        // path   : kraken2_db/
    flu_db        = ch_flu_db            // path   : blast_db/
    versions      = ch_versions          // channel: [ versions.yml ]
}


//
// PARSE_FLU_DB: Replace FASTA headers >gi|{gi}|gb|{accession}|{description} with >{accession}|{description} for easier parsing and processing
// Based on https://github.com/peterk87/nf-flu/blob/master/modules/local/misc.nf (process GUNZIP_NCBI_FLU_FASTA)
//
process PARSE_FLU_DB {
    tag "$fasta"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path fasta

    output:
    path "*.fna"       , emit: fna
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    sed -E 's/^>gi\\|[0-9]+\\|gb\\|(\\w+)\\|/>/' $fasta > influenza.parsed.fna
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>1&) | sed 's/^.*(GNU sed) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
