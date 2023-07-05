//
// Download, uncompress and prepare all necessary files
// Based on nf-core/ViralRecon subworkflow: prepare_genome.nf (https://github.com/nf-core/viralrecon/blob/master/subworkflows/local/prepare_genome.nf)
//

include { UNTAR as UNTAR_KRAKEN2_DB     } from '../../modules/nf-core/untar/main'
include { KRAKEN2_BUILD                 } from '../../modules/local/kraken2_build'

workflow PREPARE_ENVIRONMENT {
    main:

    ch_versions = Channel.empty()

    //
    // Prepare files to build humanDB for Kraken2
    //
    ch_kraken2_db = Channel.empty()
    if (!params.skip_kraken2) {
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
            ch_kraken2_db = KRAKEN2_BUILD.out.db.first()
            ch_versions   = ch_versions.mix(KRAKEN2_BUILD.out.versions)
        }
    }

    emit:
    kraken2_db = ch_kraken2_db // path: kraken2_db/
    versions   = ch_versions   // channel: [ versions.yml ]
}
