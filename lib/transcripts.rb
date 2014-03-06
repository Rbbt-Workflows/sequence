require 'parallel'
module Sequence

  desc "Transcript offsets of genomic prositions. transcript:offset:strand"
  input :organism, :string, "Organism code", "Hsa"
  input :positions, :array, "Mutation Chr:Position. Separator can be ':', space or tab. Extra fields are ignored"
  def self.transcript_offsets_for_genomic_positions(organism, positions)
    log :exons, "Find exons at genomic positions"
    exons = exons_at_genomic_positions(organism, positions)

    exon_info = exon_info(organism) 
    exon_transcript_offsets = exon_transcript_offsets(organism) 

    field_positions = ["Exon Strand", "Exon Chr Start", "Exon Chr End"].collect{|field| exon_info.identify_field field}

    log :exon_offsets, "Find offsets inside exons"
    exon_offsets = exons.collect do |position, exons|
      chr, pos = position.split(/[\s:\t]/).values_at 0, 1
      chr.sub!(/chr/,'')
      pos = pos.to_i
      list = exons.collect do |exon|
        strand, start, eend = exon_info[exon].values_at *field_positions
        if strand == "1"
          offset = pos - start.to_i
        else
          offset = eend.to_i - pos
        end

        [exon, offset, strand]
      end
      [position, list]
    end

    tsv = {}
    
    exon_transcript_offsets.unnamed = false
    log :transcript_offsets, "Find offsets inside transcripts"
    exon_offsets.each do |position, list|
      next if list.empty?
      offsets = []
      list.each do |exon, offset, strand|
        Misc.zip_fields(exon_transcript_offsets[exon]).each do |transcript, exon_offset|
          offsets << [transcript, exon_offset.to_i + offset, strand] * ":"
        end if exon_transcript_offsets.include? exon
      end
      tsv[position] = offsets
    end

    TSV.setup(tsv, :key_field => "Genomic Position", :fields => ["Ensembl Transcrip ID:Offset:Strand"], :type => :flat, :namespace => organism, :unnamed => true)
  end
  task :transcript_offsets_for_genomic_positions => :tsv
  export_synchronous :transcript_offsets_for_genomic_positions

  #{{{ CODONS

  desc "Return transcript codon at the specified position. Codon:offset:pos"
  input :organism, :string, "Organism code", "Hsa"
  input :transcript, :string, "Ensembl Transcript ID"
  input :offset, :integer, "Offset inside transcript"
  def self.codon_at_transcript_position(organism, transcript, offset)
    transcript_sequence = transcript_sequence(organism)
    transcript_5utr = transcript_5utr(organism) 
    transcript_3utr = transcript_3utr(organism) 
    transcript_phase = transcript_phase(organism)

    utr5 = transcript_5utr[transcript]
      
    if utr5.nil? or utr5 == 0 
      phase = transcript_phase[transcript]
      raise "No UTR5 and no phase for transcript: #{ transcript }" if phase.nil?
      phase = phase.to_i
      raise "No UTR5 but phase is -1: #{ transcript }" if phase == -1
      utr5 = - phase
    else
      utr5 = utr5.to_i
    end

    return "UTR5" if utr5 > offset

    sequence = transcript_sequence[transcript]
    raise "Sequence for transcript was missing: #{ transcript }" if sequence.nil? 

    ccds_offset = offset - utr5
    utr3 = transcript_3utr[transcript].to_i
    # TODO: Check this is ok!
    return "UTR3" if ccds_offset > (sequence.length - utr3 - utr5)

    if utr5 >= 0
      range = (utr5..-1)
      sequence = sequence[range]
    else
      sequence = "N" * utr5.abs << sequence
    end

    codon = ccds_offset / 3
    codon_offset =  ccds_offset % 3

    [sequence[(codon * 3)..((codon + 1) * 3 - 1)], codon_offset, codon] * ":"
  end
  task :codon_at_transcript_position => :string
  export_exec :codon_at_transcript_position
end
