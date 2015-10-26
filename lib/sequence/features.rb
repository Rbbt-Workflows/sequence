module Sequence

  input *POSITIONS_INPUT
  input *ORGANISM_INPUT
  input *VCF_INPUT
  dep &VCF_CONVERTER
  task :genes => :tsv do |positions,organism|
    begin 
      step(:genomic_mutations)
      Misc.consume_stream positions, true
      positions = step(:genomic_mutations)
    rescue
    end
    raise ParameterException, "No 'positions' specified" if positions.nil?
    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => organism
    dumper.init
    chromosome_files = {}
    TSV.traverse positions, :type => :array, :into => dumper, :bar => "Overlapping genes" do |position|
      chr, pos = position.split(/[\s:\t]+/)
      next if pos.nil?
      chr.sub!(/^chr/i,'')
      index = chromosome_files[chr] ||= Sequence.gene_chromosome_index(organism, chr)
      genes = index[pos.to_i]
      [position, genes]
    end
  end
  export_synchronous :genes

  input *POSITIONS_INPUT
  input *ORGANISM_INPUT
  input *VCF_INPUT
  dep &VCF_CONVERTER
  task :exons => :tsv do |positions,organism|
    begin
      step(:genomic_mutations)
      Misc.consume_stream positions, true
      positions = step(:genomic_mutations)
    rescue
    end
    raise ParameterException, "No 'positions' specified" if positions.nil?
    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Ensembl Exon ID"], :type => :flat, :namespace => organism
    dumper.init
    chromosome_files = {}
    exon_position = Sequence.exon_position(organism)
    TSV.traverse positions, :bar => "Exons", :type => :array, :into => dumper do |position|
      chr, pos = position.split(/[\s:\t]+/)
      next if pos.nil?
      chr.sub!(/^chr/i,'')
      index = chromosome_files[chr] ||= Sequence.exon_chromosome_index(organism, chr)
      exons = index[pos.to_i]
      next if exons.nil? or exons.empty?
      [position, exons]
    end
  end
  export_synchronous :exons

  input *POSITIONS_INPUT
  input *ORGANISM_INPUT
  input *VCF_INPUT
  dep &VCF_CONVERTER
  task :transcripts => :tsv do |positions,organism|
    begin
      step(:genomic_mutations)
      Misc.consume_stream positions, true
      positions = step(:genomic_mutations)
    rescue
    end
    raise ParameterException, "No 'positions' specified" if positions.nil?
    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Ensembl Transcript ID"], :type => :flat, :namespace => organism
    dumper.init
    chromosome_files = {}
    TSV.traverse positions, :type => :array, :into => dumper do |position|
      chr, pos = position.split(/[\s:\t]+/)
      next if pos.nil?
      chr.sub!(/^chr/i,'')
      index = chromosome_files[chr] ||= Sequence.transcript_chromosome_index(organism, chr)
      transcripts = index[pos.to_i]
      [position, transcripts]
    end
  end
  export_synchronous :transcripts

  input *POSITIONS_INPUT
  input *ORGANISM_INPUT
  input *VCF_INPUT
  dep &VCF_CONVERTER
  task :exon_junctions => :tsv do |positions,organism|
    begin
      step(:genomic_mutations)
      Misc.consume_stream positions, true
      positions = step(:genomic_mutations)
    rescue
    end

    raise ParameterException, "No 'positions' specified" if positions.nil?

    exon_position = Sequence.exon_position(organism)
    chromosome_files_start = {}
    chromosome_files_end = {}

    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Ensembl Exon ID"], :type => :flat, :namespace => organism
    dumper.init
    TSV.traverse positions, :bar => "Junctions", :type => :array, :into => dumper do |position|
      chromosome, pos = position.split ":"
      next if pos.nil?
      chromosome.sub!(/^chr/i,'')
      pos = pos.to_i
      junctions = []
      start_index = chromosome_files_start[chromosome] ||= Sequence.exon_start_index(organism, chromosome)
      end_index = chromosome_files_end[chromosome] ||= Sequence.exon_end_index(organism, chromosome)

      end_exons = end_index[pos - 10..pos + 10]
      start_exons = start_index[pos - 10..pos + 10]

      end_exons.each do |exon|
        next unless exon_position.include? exon
        strand, _start, eend = exon_position[exon]
        diff = pos - eend
        case
        when (strand == 1 and diff <= 8 and diff >= -2)
          junctions << exon + ":donor(#{diff})"
        when (strand == -1 and diff <= 8 and diff >= -3)
          junctions << exon + ":acceptor(#{diff})"
        end
      end

      start_exons.each do |exon|
        next unless exon_position.include? exon
        strand, start = exon_position[exon]
        diff = pos - start

        case
        when (strand == 1 and diff >= -8 and diff <= 2)
          junctions << exon + ":acceptor(#{diff})"
        when (strand == -1 and diff >= -8 and diff <= 3)
          junctions << exon + ":donor(#{diff})"
        end
      end
      next if junctions.empty?

      [position, junctions]
    end
  end
  export_synchronous :exon_junctions

  input *RANGES_INPUT
  input *ORGANISM_INPUT
  task :genes_at_ranges => :tsv do |ranges,organism|
    chromosome_files = {}
    dumper = TSV::Dumper.new :key_field => "Chromosome Range", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => organism
    dumper.init
    TSV.traverse ranges, :type => :array, :into => dumper do |range|
      chr, start, eend = range.split(/[\s:\t]+/)
      next if eend.nil?
      chr.sub!(/^chr/i,'')
      index = chromosome_files[chr] ||= Sequence.gene_chromosome_index(organism, chr)
      genes = index[(start.to_i..eend.to_i)]
      [range, genes]
    end
  end
  export_synchronous :genes_at_ranges
end
