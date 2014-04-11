module Sequence

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

  def self.codon_at_transcript_position(organism, transcript, offset)
    transcript_sequence = Sequence.transcript_sequence(organism)
    transcript_5utr = Sequence.transcript_5utr(organism) 
    transcript_3utr = Sequence.transcript_3utr(organism) 
    transcript_phase = Sequence.transcript_phase(organism)

    utr5 = transcript_5utr[transcript]
      
    if utr5.nil? or utr5 == "0" 
      phase = transcript_phase[transcript]
      raise TranscriptError, "No UTR5 and no phase for transcript: #{ transcript }" if phase.nil?
      phase = phase.to_i
      raise TranscriptError, "No UTR5 but phase is -1: #{ transcript }" if phase == -1
      utr5 = - phase
    else
      utr5 = utr5.to_i
    end

    return "UTR5" if utr5 > offset

    sequence = transcript_sequence[transcript]
    raise "Sequence for transcript was missing: #{ transcript }" if sequence.nil? 

    ccds_offset = offset - utr5
    utr3 = transcript_3utr[transcript].to_i
    utr3 = utr3.to_i

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

  input *POSITIONS_INPUT
  input *ORGANISM_INPUT
  dep :exons
  task :transcript_offsets => :tsv do |positions,organism|
    mutations = step(:genomic_mutations) if step(:genomic_mutations)

    exon_position = Sequence.exon_position(organism)
    exon_transcript_offsets = Sequence.exon_transcript_offsets(organism)

    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Transcript position"], :type => :flat, :namespace => organism
    dumper.init
    
    TSV.traverse step(:exons), :_cpus => 2, :into => dumper, :type => :flat do |position,exons|
      next if position.nil?
      pos = position.split(":")[1]
      next if pos.nil?
      pos = pos.to_i
      if exons.nil?  or exons.empty?
        [position, []]
      else
        transcript_positions = exons.inject([]) do |offsets,exon|
          next offsets unless exon_position.include? exon
          strand, start, eend = exon_position[exon]
          if strand == 1
            offset = pos - start
          else
            offset = eend - pos
          end
          Misc.zip_fields(exon_transcript_offsets[exon]).each do |transcript, exon_offset|
            offsets << [transcript, exon_offset.to_i + offset, strand] * ":"
          end if exon_transcript_offsets.include? exon
          offsets
        end.compact

        [position, transcript_positions]
      end
    end
  end
  export_synchronous :transcript_offsets

  input *MUTATIONS_INPUT
  input *ORGANISM_INPUT
  input *WATSON_INPUT
  input *VCF_INPUT 
  dep do |jobname, options|
    options[:positions] = options[:mutations]
    Sequence.job(:transcript_offsets, jobname, options)
  end
  task :mutated_isoforms => :tsv do |mutations,organism,watson|

    transcript_protein = Sequence.transcript_protein(organism)

    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Mutated Isoform"], :type => :flat, :namespace => organism
    dumper.init
    TSV.traverse step(:transcript_offsets), :_cpus => 2, :into => dumper, :type => :flat do |mutation,transcript_offsets|
      next if mutation.nil?
      chr, pos, mut_str = mutation.split(":")
      next if mut_str.nil?
      chr.sub!(/^chr/i,'')
      pos = pos.to_i

      mis = []
      mut_str.split(',').each do |mut|
        alleles = Sequence.alleles mut

        transcript_offsets.collect{|to| to.split ":" }.each do |transcript, transcript_offset, strand|
          protein = transcript_protein[transcript]
          #eee [:missing_protein, transcript, protein] if protein.nil? or protein.strip.empty?
          next if protein.nil? or protein.empty?

          begin
            codon = Sequence.codon_at_transcript_position(organism, transcript, transcript_offset.to_i);

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
  end
  export_synchronous :mutated_isoforms
end
