require 'rbbt'
require 'rbbt/workflow'

require 'rbbt/util/simpleopt'
$cpus ||= SOPT.get("--cpus*")[:cpus]
$cpus = $cpus.to_i if $cpus

Log.info "Loading Structure with #{ $cpus.inspect }" unless $cpus.nil?

module Sequence
  extend Workflow

  input :genomic_mutations, :array, "chr:pos:mutation"
  input :organism, :string, "Organism code", "Hsa"
  input :watson, :boolean, "Mutations all reported on the watson (forward) strand as opposed to the gene strand", true
  task :mutated_isoforms => :tsv do |genomic_mutations,organism,watson|
    exon_info = Sequence.exon_info(organism) 
    exon_info_field_positions = ["Exon Strand", "Exon Chr Start", "Exon Chr End"].collect{|field| exon_info.identify_field field}

    exon_transcript_offsets = Sequence.exon_transcript_offsets(organism) 
    transcript_to_protein = Sequence.transcript_protein(organism)

    tsv = {}
    TSV.traverse genomic_mutations, :cpus => $cpus, :into => tsv do |mutation|
      chr, pos, mut_str = mutation.split(":")
      next if mut_str.nil?
      chr.sub!(/^chr/i,'')
      pos = pos.to_i

      exons = Sequence.exon_position_index(organism, chr)[pos]

      transcript_offsets = []
      exons.each do |exon|
        strand, start, eend = exon_info[exon].values_at *exon_info_field_positions

        if strand == "1"
          offset = pos - start.to_i
        else
          offset = eend.to_i - pos
        end

        Misc.zip_fields(exon_transcript_offsets[exon]).each do |transcript, exon_offset|
          transcript_offsets << [transcript, exon_offset.to_i + offset, strand]
        end if exon_transcript_offsets.include? exon
      end

      isoforms = []
      mis = []

      mut_str.split(',').each do |mut|
        alleles = Sequence.alleles mut

        transcript_offsets.each do |transcript, offset, strand|
          protein = transcript_to_protein[transcript]
          next if protein.nil?

          begin
            codon = Sequence.codon_at_transcript_position(organism, transcript, offset)

            case codon

            when "UTR5", "UTR3"
              mis << [transcript, codon] * ":"

            else # Protein mutation
              triplet, offset, pos = codon.split ":"
              next if not triplet.length == 3
              original = Misc::CODON_TABLE[triplet]
              next if alleles.empty?
              pos = pos.to_i
              alleles.each do |allele|
                change = case allele
                         when "Indel"
                           [original, pos + 1, "Indel"] * ""
                         when "FrameShift"
                           [original, pos + 1, "FrameShift"] * ""
                         else
                           allele = Misc::BASE2COMPLEMENT[allele] if watson and strand == "-1"
                           triplet[offset.to_i] = allele 
                           new = Misc::CODON_TABLE[triplet]
                           [original, pos + 1, new] * ""
                         end
                mis << [protein, change] * ":"
              end
            end
          rescue TranscriptError
            Log.debug{$!.message}
          end
        end
      end
      next if mis.empty?

      [mutation, mis]
    end

    TSV.setup(tsv, :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"], :type => :flat, :namespace => organism)
  end
  export_asynchronous :mutated_isoforms
end

require 'indices'
require 'alleles'
require 'genes'
require 'exons'
require 'transcripts'
require 'mutated_isoforms'
require 'liftover'
