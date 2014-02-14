module Sequence

  desc "Change mutation reference strand from 'gene' strand to 'watson' strand"
  input :organism, :string, "Organism code", "Hsa"
  input :mutations, :array, "Mutation Chr:Position:Mut (e.g. 19:54646887:A). Separator can be ':', space or tab. Extra fields are ignored"
  def self.to_watson(organism, mutations)
    transcript_offsets = transcript_offsets_for_genomic_positions(organism, mutations)

    fixed = {}
    mutations.each{|mutation| fixed[mutation] = mutation}

    transcript_offsets.each do |mutation, list|
      chr, pos, mut = mutation.split ":"
      next unless Misc::BASE2COMPLEMENT.include? mut
      fixed[mutation] = mutation.sub(mut, Misc::BASE2COMPLEMENT[mut]) if (list.any? and list.first.split(":")[2] == "-1")
    end

    fixed.chunked_values_at mutations
  end
  task :to_watson => :array
  export_synchronous :to_watson

  desc "Guess if mutations are given in watson or gene strand"
  input :organism, :string, "Organism code", "Hsa"
  input :mutations, :array, "Mutation Chr:Position:Mut (e.g. 19:54646887:A). Separator can be ':', space or tab. Extra fields are ignored"
  def self.is_watson(organism, mutations)
    diffs = (mutations - to_watson(organism, mutations)).each

    same = 0
    opposite = 0
    diffs.zip(reference_allele_at_genomic_positions(organism, diffs).values_at *diffs).each do |mutation, reference|
      chr, pos, mut = mutation.split ":"
      same += 1 if mut == reference
      opposite += 1 if mut == Misc::BASE2COMPLEMENT[reference]
    end

    log(:counts, "Opposite: #{ opposite }. Same: #{ same }")
    opposite > same
  end
  task :is_watson => :boolean
  export_synchronous :is_watson

  desc "Reference allele at positions in a chromosome"
  input :organism, :string, "Organism code", "Hsa"
  input :chromosome, :string, "Chromosome name"
  input :positions, :array, "Positions"
  def self.reference_allele_at_chr_positions(organism, chromosome, positions)
    begin
      chromosome_file = Organism[File.join(organism, "chromosome_#{chromosome}")].produce.find
      File.open(chromosome_file) do |f|
        Misc.process_to_hash(positions.sort){|list| list.collect{|position| f.seek(position.to_i - 1); c = f.getc; c.nil? ? nil : c.chr }}.chunked_values_at positions
      end
    rescue
      if $!.message =~ /Fasta file for chromosome not found/i or $!.message =~ /No such file or directory/i
        Log.low{$!.message}
        ["?"] * positions.length
      else
        raise $!
      end
    end
  end
  task :reference_allele_at_chr_positions => :array
  export_exec :reference_allele_at_chr_positions

  desc "Reference allele at genomic positions"
  input :organism, :string, "Organism code", "Hsa"
  input :positions, :array, "Positions Chr:Position (e.g. 19:54646887). Separator can be ':', space or tab. Extra fields are ignored"
  def self.reference_allele_at_genomic_positions(organism, positions)
    chr_positions = {}

    log :parsing, "Parsing positions"
    positions.each do |position|
      chr, pos = position.split(/[\s:\t]/)
      chr.sub!(/chr/,'')
      chr_positions[chr] ||= []
      chr_positions[chr] << pos
    end

    log :processing, "Processing chromosome positions"
    chr_bases = {}
    chr_positions.each do |chr, list|
      chr_bases[chr] = reference_allele_at_chr_positions(organism, chr, list)
    end

    log :loading, "Loading results"
    tsv = TSV.setup({}, :key_field => "Genomic Position", :fields => ["Reference Allele"], :type => :single, :namespace => organism, :unnamed => true)
    positions.collect do |position|
      chr, pos = position.split(/[\s:\t]/).values_at 0, 1
      chr.sub!(/chr/,'')
      tsv[position] = chr_bases[chr].shift
    end

    tsv
  end
  task :reference_allele_at_genomic_positions=> :tsv
  export_synchronous :reference_allele_at_genomic_positions
end
