//
// Run BLASTn on previous assembly step using NCBI Influenza Virus Database
//

include { BLAST_BLASTN                } from '../../modules/nf-core/blast/blastn/main'

workflow BLASTN_AND_PARSE {
    take:
    fasta     // channel: [ val(meta), [ scaffolds ] ]
    blast_db  // channel: /path/to/blast_db/

    main:
    ch_versions = Channel.empty()

    //
    // Run BLASTn os assemblies from Unicycler
    //
    ch_blast_txt        = Channel.empty()
    //ch_blast_filter_txt = Channel.empty()

    BLAST_BLASTN (
        fasta,
        blast_db
    )
    ch_blast_txt = BLAST_BLASTN.out.txt
    ch_versions  = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    emit:

    versions = ch_versions
}
