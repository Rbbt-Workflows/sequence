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
    transcript_sequence = @@transcript_sequence ||= Sequence.transcript_sequence(organism)
    transcript_5utr = @@transcript_5utr ||= Sequence.transcript_5utr(organism) 
    transcript_3utr = @@transcript_3utr ||= Sequence.transcript_3utr(organism) 
    transcript_phase = @@transcript_phase ||= Sequence.transcript_phase(organism)

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

  dep do |jobname, options|
    options = options.dup
    IndiferentHash.setup options
    options[:positions] = options.delete :mutations 
    if FalseClass === options[:watson]
      Sequence.job(:gene_strand_reference, jobname, options)
    else
      Sequence.job(:reference, jobname, options)
    end
  end
  input *MUTATIONS_INPUT
  input *ORGANISM_INPUT
  input *WATSON_INPUT
  input *VCF_INPUT
  task :type => :tsv do |_muts, _org, watson|
    Misc.consume_stream _muts, true
    reference_job = watson ? step(:reference) : step(:gene_strand_reference)

    reference_job.grace
    organism = reference_job.info[:inputs][:organism]
    mutation_type = TSV::Dumper.new(:key_field => "Genomic Mutation", :fields => ["Mutation type"], :type => :single, :namespace => organism)
    mutation_type.init
    TSV.traverse reference_job, :bar => "Type", :into => mutation_type do |mutation, reference|
      mutation = mutation.first if Array === mutation
      base = mutation.split(":")[2]

      type = case
             when (base.nil? or reference.nil? or base == "?" or reference == "?")
               "unknown"
             when base.index(',')
               "multiple"
             when base == reference
               "none"
             when (base.length > 1 or base == '-')
               "indel"
             when (not %w(A G T C).include? base and not %w(A G T C).include? reference) 
               "unknown"
             when (((Misc::IUPAC2BASE[base] || []) & ["A", "G"]).any? and     ((Misc::IUPAC2BASE[reference] || []) & ["T", "C"]).any?)
               "transversion"
             when (((Misc::IUPAC2BASE[base] || []) & ["T", "C"]).any? and     ((Misc::IUPAC2BASE[reference] || []) & ["A", "G"]).any?)
               "transversion"
             when (((Misc::IUPAC2BASE[base] || []) & ["A", "G"]).any? and     ((Misc::IUPAC2BASE[reference] || [nil]) & ["T", "C", nil]).empty?)
               "transition"
             when (((Misc::IUPAC2BASE[base] || []) & ["T", "C"]).any? and     ((Misc::IUPAC2BASE[reference] || [nil]) & ["A", "G", nil]).empty?)
               "transition"
             else
               "unknown"
             end
      [mutation, type]
    end
  end
  export_asynchronous :type

  dep :exons
  task :transcript_offsets => :tsv do
    mutations = step(:genomic_mutations) if step(:genomic_mutations)
    organism = step(:exons).inputs[1]

    exon_position = Sequence.exon_position(organism)
    exon_transcript_offsets = Sequence.exon_transcript_offsets(organism)

    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Transcript position"], :type => :flat, :namespace => organism
    dumper.init
    
    TSV.traverse step(:exons), :bar => "Transcript offsets", :into => dumper, :type => :flat do |position,exons|
      next if position.nil?
      pos = position.split(":")[1]
      next if pos.nil?
      pos = pos.to_i
      next if exons.nil?  or exons.empty?

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

      next if transcript_positions.empty?

      [position, transcript_positions]
    end
  end
  export_synchronous :transcript_offsets

  input *MUTATIONS_INPUT
  input *ORGANISM_INPUT
  input *WATSON_INPUT
  input *VCF_INPUT
  input *PRINCIPAL_INPUT
  dep do |jobname, options|
    options[:positions] = options[:mutations]
    Sequence.job(:transcript_offsets, jobname, options)
  end
  task :mutated_isoforms => :tsv do |mutations,organism,watson,vcf,principal|
    Misc.consume_stream mutations, true

    transcript_protein = Sequence.transcript_protein(organism)

    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"], :type => :flat, :namespace => organism
    dumper.init
    TSV.traverse step(:transcript_offsets), :bar => "Mutated Isoforms", :_cpus => 2, :into => dumper, :type => :flat do |mutation,transcript_offsets|
      next if mutation.nil?
      chr, pos, mut_str = mutation.split(":")
      next if mut_str.nil?
      chr.sub!(/^chr/i,'')
      pos = pos.to_i

      mis = []
      mut_str.split(',').each do |mut|
        alleles = Sequence.alleles mut

        transcript_offsets.collect{|to| to.split ":" }.each do |transcript, transcript_offset, strand|
          next if principal and not Appris::PRINCIPAL_TRANSCRIPTS.include?(transcript)
          protein = transcript_protein[transcript]
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


  dep do |jobname,options|
    options = options.dup
    IndiferentHash.setup options
    options.merge!(:positions => options[:mutations])
    Sequence.job(:exon_junctions, jobname, options)
  end
  dep :type
  input *MUTATIONS_INPUT
  input *ORGANISM_INPUT
  input *WATSON_INPUT
  input *VCF_INPUT
  task :splicing_mutations => :tsv do |_pos|
    Misc.consume_stream _pos, true
    type = step(:type)
    exon_junctions = step(:exon_junctions)

    type.grace
    exon_junctions.grace

    organism = exon_junctions.info[:inputs][:organism]
    transcript_exons = Sequence.transcript_exons(organism)
    exon_transcripts = Sequence.exon_transcripts(organism)
    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Affected Transcripts"], :namespace => organism, :type => :flat
    dumper.init
    TSV.traverse TSV.paste_streams([type, exon_junctions], :sort => false), :bar => "Splicing Mutations", :type => :array, :into => dumper do |line|
      mutation, type, *exon_junctions = line.split "\t"
      next if type == "none" or type == "unknown"
      next if exon_junctions.empty?
      affected_transcripts = []
      exon_junctions.each do |junction|
        exon, _sep, junction_type = junction.partition ":"
        transcripts = exon_transcripts[exon]
        next if transcripts.nil?

        transcripts.select! do |transcript|
          transcript_info = transcript_exons[transcript]
          next if transcript_info.nil?

          total_exons = transcript_info[0].length
          next if transcript_info[0].index(exon).nil?
          rank = transcript_info[1][transcript_info[0].index(exon)].to_i

          case
          when (rank == 1 and junction_type =~ /acceptor/)
            false
          when (rank == total_exons and junction_type =~ /donor/)
            false
          else
            true
          end
        end
        affected_transcripts.concat transcripts

      end
      next if  affected_transcripts.empty?
      [mutation, affected_transcripts]
    end
  end
  export_asynchronous :splicing_mutations
end
