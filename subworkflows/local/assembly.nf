//
// Preliminary assembly using Unicycler (SPAdes)
//

include { UNICYCLER              } from '../../modules/nf-core/unicycler/main'
include { GUNZIP as GUNZIP_FASTA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFA   } from '../../modules/nf-core/gunzip/main'
include { BANDAGE_IMAGE          } from '../modules/nf-core/bandage/image/main'

workflow ASSEMBLY {
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

    //
    // Generate assembly visualisation with Bandage
    //
    ch_bandage_png = Channel.empty()
    ch_bandage_svg = Channel.empty()
    BANDAGE_IMAGE (
        ch_gfa
    )
    ch_bandage_png = BANDAGE_IMAGE.out.png
    ch_bandage_svg = BANDAGE_IMAGE.out.svg
    ch_versions    = ch_versions.mix(BANDAGE_IMAGE.out.versions.first())

    emit:
    fasta    = ch_fasta
    gfa      = ch_gfa

    versions = ch_versions
}
