//
// Preliminary assembly using Unicycler (SPAdes)
//

include { UNICYCLER              } from '../../modules/nf-core/unicycler/main'
include { GUNZIP as GUNZIP_FASTA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFA   } from '../../modules/nf-core/gunzip/main'
include { QUAST                  } from '../../modules/nf-core/quast/main'

workflow assembly {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    ch_fasta_gz = Channel.empty()
    UNICYCLER (
        reads
    )
    ch_fasta_gz = UNICYCLER.out.scaffolds
    ch_gfa_gz   = UNICYCLER.out.gfa
    ch_versions = ch_versions.mix(UNICYCLER.out.versions.first())

    GUNZIP_FASTA (
        ch_fasta_gz
    )
    ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions.first())

    GUNZIP_GFA (
        ch_gfa_gz
    )

    //
    // Filter for empty fasta files
    //
    GUNZIP_FASTA
        .out
        .gunzip
        .filter { meta, fasta -> fasta.size() > 0 }
        .set { ch_fasta }

    GUNZIP_GFA
        .out
        .gunzip
        .filter { meta, gfa -> gfa.size() > 0 }
        .set { ch_gfa }

    emit:
    fasta    = ch_fasta
    gfa      = ch_gfa

    versions = ch_versions
}
