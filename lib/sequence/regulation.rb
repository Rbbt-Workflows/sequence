module Sequence

  input *POSITIONS_INPUT
  input *ORGANISM_INPUT
  input *VCF_INPUT
  dep &VCF_CONVERTER
  input :distance, :integer, "Distance to TSS", 1000
  task :TSS => :tsv do |positions,organism,vcf,distance|
    begin
      step(:genomic_mutations)
      Misc.consume_stream positions, true
      positions = step(:genomic_mutations)
    rescue
    end
    raise ParameterException, "No 'positions' specified" if positions.nil?

    gene_position = Sequence.gene_position(organism)
    chromosome_files_start = {}
    chromosome_files_end = {}

    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => organism
    dumper.init
    TSV.traverse positions, :bar => "TSS", :type => :array, :into => dumper do |position|
      chromosome, pos = position.split ":"
      next if pos.nil?
      chromosome.sub!(/^chr/i,'')
      pos = pos.to_i
      overlaps = []
      start_index = chromosome_files_start[chromosome] ||= Sequence.gene_start_index(organism, chromosome)
      end_index = chromosome_files_end[chromosome] ||= Sequence.gene_end_index(organism, chromosome)

      start_genes = start_index[pos-distance..pos + distance]

      start_genes.each do |gene|
        next unless gene_position.include? gene
        strand, start, eend = gene_position[gene]
        diff = pos - start

        case
        when (strand == 1 and pos <= start)
          overlaps << gene 
        when (strand == -1 and pos >= eend)
          overlaps << gene
        end
      end
      next if overlaps.empty?

      [position, overlaps.uniq]
    end
  end
  export_synchronous :TSS

  input *RANGES_INPUT
  input *ORGANISM_INPUT
  input :distance, :integer, "Distance to TSS", 1000
  task :TSS_in_range => :tsv do |ranges,organism,distance|
    raise ParameterException, "No 'ranges' specified" if ranges.nil?

    gene_position = Sequence.gene_position(organism)
    chromosome_files_start = {}
    chromosome_files_end = {}

    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => organism
    dumper.init
    TSV.traverse ranges, :bar => "TSS", :type => :array, :into => dumper do |range|
      chromosome, rstart,rend = range.split ":"
      next if rstart.nil? or rend.nil?
      chromosome.sub!(/^chr/i,'')
      rstart = rstart.to_i
      rend = rend.to_i
      overlaps = []
      start_index = chromosome_files_start[chromosome] ||= Sequence.gene_start_index(organism, chromosome)
      end_index = chromosome_files_end[chromosome] ||= Sequence.gene_end_index(organism, chromosome)

      overlap_genes = start_index[rstart-distance..rend+distance] + end_index[rstart-distance..rend+distance]

      overlap_genes.each do |gene|
        next unless gene_position.include? gene
        strand, start, eend = gene_position[gene]

        case
        when (strand == 1 and rstart <= start)
          overlaps << gene 
        when (strand == -1 and rend >= eend)
          overlaps << gene
        end
      end
      next if overlaps.empty?

      [range, overlaps.uniq]
    end
  end
  export_synchronous :TSS_in_range
end
