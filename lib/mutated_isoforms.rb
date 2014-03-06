require 'bio'

module Sequence

  desc "Mutated protein isoforms"
  input :organism, :string, "Organism code", "Hsa"
  input :watson, :boolean, "Alleles reported always in the Watson strand (as opposed to the gene's strand)", true
  input :mutations, :array, "Mutation Chr:Position:Mut (e.g. 19:54646887:A). Separator can be ':', space or tab. Extra fields are ignored"
  def self.mutated_isoforms_for_genomic_mutations(organism, watson, mutations)

    log :transcript_offsets, "Finding transcript and offsets for mutations"
    transcript_offsets = transcript_offsets_for_genomic_positions(organism, mutations)
    transcript_to_protein = transcript_protein(organism)

    mutated_isoforms = {}

    log :translating, "Translating changes in transcript offsets into changes in proteins"

    transcript_offsets.monitor = true
    transcript_offsets.to_hash.each do |mutation, list|
      chr, pos, mut_str = mutation.split ":"
      chr.sub!(/chr/,'')

      isoforms = []
      next if mut_str.nil?
      mut_str.split(',').each do |mut|
        case
        when mut.nil?
          alleles = []
        when (mut.length == 1 and mut != '-') #A
          alleles = Misc.IUPAC_to_base(mut) || []
        when (mut[0] == "+"[0] and mut.length % 4 == 0) #+ATG
          alleles = ["Indel"]
        when (mut =~ /^-*$/ and mut.length % 3 == 0) #---
          alleles = ["Indel"]
        when ((_dash = mut.scan('-').length) % 3 == (mut.length - _dash)  % 3) #---
          alleles = ["Indel"]
        when (mut.match(/^[ATCG]+$/) and mut.length > 2 and mut.length % 4 == 0) #GATG where G is the reference
          alleles = ["Indel"]
        when (mut[0] != "-"[0] and mut[1] == "-"[0] and mut.length % 4 == 0) #G---
          alleles = ["Indel"]
        else #+A - GT etc
          alleles = ["FrameShift"]
        end

        list.collect{|t| t.split ":"}.each do |transcript, offset, strand|
          offset = offset.to_i
          begin
            codon = codon_at_transcript_position(organism, transcript, offset)
            case codon
            when "UTR5", "UTR3"
              isoforms << [transcript, codon]
            else
              triplet, offset, pos = codon.split ":"
              next if not triplet.length === 3
              #original = Bio::Sequence::NA.new(triplet).translate
              original = Misc::CODON_TABLE[triplet]
              alleles.each do |allele|
                case allele
                when "Indel"
                  isoforms << [transcript, [original, pos.to_i + 1, "Indel"] * ""]
                when "FrameShift"
                  isoforms << [transcript, [original, pos.to_i + 1, "FrameShift"] * ""]
                else
                  allele = Misc::BASE2COMPLEMENT[allele] if watson and strand.to_i == -1
                  triplet[offset.to_i] = allele 
                  #new = Bio::Sequence::NA.new(triplet).translate
                  new = Misc::CODON_TABLE[triplet]
                  isoforms << [transcript, [original, pos.to_i + 1, new] * ""]
                end
              end
            end
          rescue
            Log.debug{$!.message}
          end
        end
      end

      mutated_isoforms[mutation] = isoforms.collect{|transcript, change| 
        if change =~ /^UTR/
          [transcript, change] * ":"
        else
          protein = transcript_to_protein[transcript]
          next if protein.nil? or protein.empty?
          [protein, change] * ":"
        end
      }.compact
    end

    TSV.setup(mutated_isoforms, :type => :flat, :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"], :namespace => organism, :unnamed => true)
  end
  task :mutated_isoforms_for_genomic_mutations => :tsv
  export_synchronous :mutated_isoforms_for_genomic_mutations






  # STOP














  def self.alleles(mut)
    case
    when (mut.length == 1 and mut != '-') #A
      Misc.IUPAC_to_base(mut) || []
    when (mut[0] == "+"[0] and mut.length % 4 == 0) #+ATG
      ["Indel"]
    when (mut =~ /^-*$/ and mut.length % 3 == 0) #---
      ["Indel"]
    when ((_dash = mut.scan('-').length) % 3 == (mut.length - _dash)  % 3) #---
      ["Indel"]
    when (mut.match(/^[ATCG]+$/) and mut.length > 2 and mut.length % 4 == 0) #GATG where G is the reference
      ["Indel"]
    when (mut[0] != "-"[0] and mut[1] == "-"[0] and mut.length % 4 == 0) #G---
      ["Indel"]
    else #+A - GT etc
      ["FrameShift"]
    end
  end

  desc "Mutated protein isoforms"
  input :organism, :string, "Organism code", "Hsa"
  input :watson, :boolean, "Alleles reported always in the Watson strand (as opposed to the gene's strand)", true
  input :mutations, :array, "Mutation Chr:Position:Mut (e.g. 19:54646887:A). Separator can be ':', space or tab. Extra fields are ignored"
  def self.mutated_isoforms_for_genomic_mutations_pp(organism, watson, mutations)

    log :transcript_offsets, "Finding transcript and offsets for mutations"
    transcript_offsets = transcript_offsets_for_genomic_positions(organism, mutations) 
    transcript_to_protein = transcript_protein(organism)

    mutated_isoforms = {}

    log :translating, "Translating changes in transcript offsets into changes in proteins"
    #transcript_offsets.to_hash.each do |mutation, list|
    transcript_offsets.to_hash.each do |mutation, list|
      chr, pos, mut_str = mutation.split ":"
      chr.sub!(/chr/,'')
      isoforms = []
      next if mut_str.nil?
      mut_str.split(',').each do |mut|
        alleles = self.alleles mut

        list.collect{|t| t.split ":"}.each do |transcript, offset, strand|
          offset = offset.to_i
          begin
            codon = codon_at_transcript_position(organism, transcript, offset)
            case codon
            when "UTR5", "UTR3"
              isoforms << [transcript, codon]
            else
              triplet, offset, pos = codon.split ":"
              next if not triplet.length === 3
              original = Misc::CODON_TABLE[triplet]
              alleles.each do |allele|
                case allele
                when "Indel"
                  isoforms << [transcript, [original, pos.to_i + 1, "Indel"] * ""]
                when "FrameShift"
                  isoforms << [transcript, [original, pos.to_i + 1, "FrameShift"] * ""]
                else
                  allele = Misc::BASE2COMPLEMENT[allele] if watson and strand.to_i == -1
                  triplet[offset.to_i] = allele 
                  #new = Bio::Sequence::NA.new(triplet).translate
                  new = Misc::CODON_TABLE[triplet]
                  isoforms << [transcript, [original, pos.to_i + 1, new] * ""]
                end
              end
            end
          rescue
            Log.debug{$!.message}
          end
        end
      end

      mutated_isoforms[mutation] = isoforms.collect{|transcript, change| 
        if change =~ /^UTR/
          [transcript, change] * ":"
        else
          protein = transcript_to_protein[transcript]
          next if protein.nil? or protein.empty?
          [protein, change] * ":"
        end
      }.compact
    end

    log :saving, "Saving results"

    TSV.setup(mutated_isoforms, :type => :flat, :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"], :namespace => organism, :unnamed => true)
  end
  task :mutated_isoforms_for_genomic_mutations_pp => :tsv
  export_synchronous :mutated_isoforms_for_genomic_mutations_pp
end
