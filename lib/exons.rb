module Sequence
  desc "Find exons at particular positions in a chromosome. Multiple values separated by '|'"
  input :organism, :string, "Organism code", "Hsa"
  input :chromosome, :string, "Chromosome name"
  input :positions, :array, "Positions"
  def self.exons_at_chr_positions(organism, chromosome, positions)
    index = exon_position_index(organism, chromosome)
    index.chunked_values_at(positions).collect{|list| list * "|"}
  end
  task :exons_at_chr_positions => :array
  export_exec :exons_at_chr_positions

  desc "Find exons at particular genomic positions. Multiple values separated by '|'"
  input :organism, :string, "Organism code", "Hsa"
  input :positions, :array, "Positions Chr:Position (e.g. 11:533766). Separator can be ':', space or tab. Extra fields are ignored"
  def self.exons_at_genomic_positions(organism, positions)
    chr_positions = {}
    positions = positions.compact.reject{|p| p.empty? }
    positions.each do |position|
      chr, pos = position.split(/[\s:\t]/).values_at 0, 1
      chr.sub!(/chr/,'')
      chr_positions[chr] ||= []
      chr_positions[chr] << pos
    end

    chr_exons = {}
    chr_positions.each do |chr, list|
      chr_exons[chr] = exons_at_chr_positions(organism, chr, list)
    end

    tsv = TSV.setup({}, :key_field => "Genomic Position", :fields => ["Ensembl Exon ID"], :type => :flat, :namespace => organism)
    positions.collect do |position|
      chr, pos = position.split(/[\s:\t]/).values_at 0, 1
      chr.sub!(/chr/,'')
      tsv[position] = chr_exons[chr].shift.split("|")
    end
    tsv
  end
  task :exons_at_genomic_positions => :tsv
  export_synchronous :exons_at_genomic_positions

  #{{{ EXON JUNCTIONS

  desc "Find exon junctions at particular positions in a chromosome. Multiple values separated by '|'"
  input :organism, :string, "Organism code", "Hsa"
  input :chromosome, :string, "Chromosome name"
  input :positions, :array, "Positions"
  def self.exon_junctions_at_chr_positions(organism, chromosome, positions)
    start_index = exon_start_index(organism, chromosome)
    end_index = exon_end_index(organism, chromosome)
    exon_info = exon_info(organism)

    strand_field_pos = exon_info.identify_field "Exon Strand"
    start_field_pos = exon_info.identify_field "Exon Chr Start"
    end_field_pos = exon_info.identify_field "Exon Chr End"
    positions.collect{|pos|
      pos = pos.to_i
      junctions = []

      end_exons = end_index[pos - 10..pos + 10]
      start_exons = start_index[pos - 10..pos + 10]

      end_exons.each do |exon|
        last_exon = exon == end_exons.last

        strand, eend = exon_info[exon].values_at strand_field_pos, end_field_pos
        eend = eend.to_i
        diff = pos - eend
        case
        when (strand == "1" and diff <= 8 and diff >= -2)
          junctions << exon + ":donor(#{diff})"
        when (strand == "-1" and diff <= 8 and diff >= -3)
          junctions << exon + ":acceptor(#{diff})"
        end
      end

      start_exons.each do |exon|
        first_exon = exon == start_exons.first

        strand, start = exon_info[exon].values_at strand_field_pos, start_field_pos
        start = start.to_i
        diff = pos - start

        case
        when (strand == "1" and diff >= -8 and diff <= 2)
          junctions << exon + ":acceptor(#{diff})"
        when (strand == "-1" and diff >= -8 and diff <= 3)
          junctions << exon + ":donor(#{diff})"
        end
      end

      junctions * "|"
    }
  end
  task :exon_junctions_at_chr_positions => :array
  export_exec :exon_junctions_at_chr_positions

  desc "Find exon junctions at particular genomic positions. Multiple values separated by '|'"
  input :organism, :string, "Organism code", "Hsa"
  input :positions, :array, "Positions Chr:Position (e.g. 11:533766). Separator can be ':', space or tab. Extra fields are ignored"
  def self.exon_junctions_at_genomic_positions(organism, positions)
    raise ParameterException, "No 'positions' specified" if positions.nil?

    tsv = TSV.setup({}, :key_field => "Genomic Position", :fields => ["Exon Junction"], :type => :flat, :namespace => organism, :unnamed => true)

    positions = positions.reject{|p| p.nil? or p.empty?}

    chr_positions = {}
    positions.each do |position|
      chr, pos = position.strip.split(/[\s:\t]/).values_at 0, 1
      chr.sub!(/chr/,'')
      chr_positions[chr] ||= []
      chr_positions[chr] << pos
    end

    chr_exon_junctions = {}
    chr_positions.each do |chr, list|
      chr_exon_junctions[chr] = exon_junctions_at_chr_positions(organism, chr, list)
    end

    positions.each do |position|
      chr, pos = position.strip.split(/[\s:\t]/).values_at 0, 1
      chr.sub!(/chr/,'')
      tsv[position] = chr_exon_junctions[chr].shift.split("|")
    end

    tsv.with_unnamed do
      tsv = tsv.select("Exon Junction"){|ej| ej.any?}
    end
    tsv
  end
  task :exon_junctions_at_genomic_positions => :tsv
  export_synchronous :exon_junctions_at_genomic_positions
end
