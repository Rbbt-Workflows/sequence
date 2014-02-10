module Sequence

  #{{{ SNPs

  desc "Identify known SNPs in a chromosome"
  input :organism, :string, "Organism code", "Hsa"
  input :chromosome, :string, "Chromosome name"
  input :positions, :array, "Positions Chr:Position (e.g. 11:533766). Separator can be ':', space or tab. Extra fields are ignored"
  def self.snps_at_chr_positions(organism, chromosome, positions)
    index = snp_position_index(organism, chromosome)
    index.values_at(*positions).collect{|list| list * "|"}
  end
  task :snps_at_chr_positions => :array
  export_exec :snps_at_chr_positions

  def self.somatic_snvs_at_chr_positions(organism, chromosome, positions)
    index = somatic_snv_position_index(organism, chromosome)
    index.values_at(*positions).collect{|list| list * "|"}
  end
  task :somatic_snvs_at_chr_positions => :array
  export_exec :somatic_snvs_at_chr_positions

  desc "Identify known SNPs at genomic positions"
  input :organism, :string, "Organism code", "Hsa"
  input :positions, :array, "Positions Chr:Position (e.g. 11:533766). Separator can be ':', space or tab. Extra fields are ignored"
  def self.snps_at_genomic_positions(organism, positions)
    chr_positions = {}
    positions.each do |position|
      chr, pos = position.split(/[\s:\t]/).values_at 0, 1
      chr_positions[chr] ||= []
      chr_positions[chr] << pos
    end

    chr_snps = {}
    chr_positions.each do |chr, list|
      chr_snps[chr] = snps_at_chr_positions(organism, chr, list)
    end

    tsv = TSV.setup({}, :key_field => "Genomic Position", :fields => ["Germline SNP"], :type => :double, :namespace => organism)
    positions.collect do |position|
      chr, pos = position.split(/[\s:\t]/).values_at 0, 1
      tsv[position] = chr_snps[chr].shift.split("|")
    end
    tsv
  end
  task :snps_at_genomic_positions => :tsv
  export_asynchronous :snps_at_genomic_positions

  desc "Identify known somatic SNVs at genomic positions"
  input :organism, :string, "Organism code", "Hsa"
  input :positions, :array, "Positions Chr:Position (e.g. 11:533766). Separator can be ':', space or tab. Extra fields are ignored"
  def self.somatic_snvs_at_genomic_positions(organism, positions)
    chr_positions = {}
    positions.each do |position|
      chr, pos = position.split(/[\s:\t]/).values_at 0, 1
      chr_positions[chr] ||= []
      chr_positions[chr] << pos
    end

    chr_somatic_snvs = {}
    chr_positions.each do |chr, list|
      chr_somatic_snvs[chr] = somatic_snvs_at_chr_positions(organism, chr, list)
    end

    tsv = TSV.setup({}, :key_field => "Genomic Position", :fields => ["Somatic SNV"], :type => :double, :namespace => organism)
    positions.collect do |position|
      chr, pos = position.split(/[\s:\t]/).values_at 0, 1
      tsv[position] = chr_somatic_snvs[chr].shift.split("|")
    end
    tsv
  end
  task :somatic_snvs_at_genomic_positions => :tsv
  export_asynchronous :somatic_snvs_at_genomic_positions


  #----- GENES

  desc "Find SNPS at particular ranges in a chromosome. Multiple values separated by '|'"
  input :organism, :string, "Organism code", "Hsa"
  input :chromosome, :string, "Chromosome name"
  input :ranges, :array, "Ranges"
  def self.snps_at_chr_ranges(organism, chromosome, ranges)
    index = snp_position_index(organism, chromosome)
    r = ranges.collect{|r| s,e = r.split(":"); (s.to_i..e.to_i)}
    index.values_at(*r).collect{|list| list * "|"}
  end
  task :snps_at_chr_ranges => :array
  export_exec :snps_at_chr_ranges

  desc "Find SNPS at particular genomic ranges. Multiple values separated by '|'"
  input :organism, :string, "Organism code", "Hsa"
  input :ranges, :array, "Positions Chr:Start:End (e.g. 11:533766:553323)"
  def self.snps_at_genomic_ranges(organism, ranges)
    chr_ranges = {}
    ranges.each do |range|
      chr, s, e = range.split(":").values_at 0, 1, 2
      chr.sub!(/chr/,'')
      chr_ranges[chr] ||= []
      chr_ranges[chr] << [s, e] * ":" 
    end

    chr_snps = {}
    chr_ranges.each do |chr, list|
      chr_snps[chr] = snps_at_chr_ranges(organism, chr, list)
    end

    tsv = TSV.setup({}, :key_field => "Genomic Range", :fields => ["SNP ID"], :type => :flat, :namespace => organism)
    ranges.collect do |range|
      chr, s, e = range.split(":").values_at 0, 1, 2
      chr.sub!(/chr/,'')
      tsv[range] = chr_snps[chr].shift.split("|")
    end
    tsv
  end
  task :snps_at_genomic_ranges => :tsv
  export_synchronous :snps_at_genomic_ranges

  desc "Find somatic SNVs at particular ranges in a chromosome. Multiple values separated by '|'"
  input :organism, :string, "Organism code", "Hsa"
  input :chromosome, :string, "Chromosome name"
  input :ranges, :array, "Ranges"
  def self.somatic_snvs_at_chr_ranges(organism, chromosome, ranges)
    index = somatic_snv_position_index(organism, chromosome)
    r = ranges.collect{|r| s,e = r.split(":"); (s.to_i..e.to_i)}
    index.values_at(*r).collect{|list| list * "|"}
  end
  task :somatic_snvs_at_chr_ranges => :array
  export_exec :somatic_snvs_at_chr_ranges

  desc "Find somatic SNVs at particular genomic ranges. Multiple values separated by '|'"
  input :organism, :string, "Organism code", "Hsa"
  input :ranges, :array, "Positions Chr:Start:End (e.g. 11:533766:553323)"
  def self.somatic_snvs_at_genomic_ranges(organism, ranges)
    chr_ranges = {}
    ranges.each do |range|
      next if range.nil? or range.empty?
      chr, s, e = range.split(":").values_at 0, 1, 2
      chr.sub!(/chr/,'')
      chr_ranges[chr] ||= []
      chr_ranges[chr] << [s, e] * ":" 
    end

    chr_somatic_snvs = {}
    chr_ranges.each do |chr, list|
      chr_somatic_snvs[chr] = somatic_snvs_at_chr_ranges(organism, chr, list)
    end

    tsv = TSV.setup({}, :key_field => "Genomic Range", :fields => ["SNP ID"], :type => :flat, :namespace => organism)
    ranges.each do |range|
      next if range.nil? or range.empty?
      chr, s, e = range.split(":").values_at 0, 1, 2
      chr.sub!(/chr/,'')
      tsv[range] = chr_somatic_snvs[chr].shift.split("|")
    end
    tsv
  end
  task :somatic_snvs_at_genomic_ranges => :tsv
  export_synchronous :somatic_snvs_at_genomic_ranges
end
