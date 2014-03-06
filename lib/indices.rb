require 'rbbt/sources/organism'
module Sequence

  CACHE = {}

  def self.gene_position_index(organism, chromosome)
    key = [organism, chromosome]
    CACHE[:gene_position] ||= {}
    if CACHE[:gene_position][key].nil?
      CACHE[:gene_position][key] = TSV.range_index(Organism.gene_positions(organism).produce, "Gene Start", "Gene End", :filters => [["field:Chromosome Name", chromosome]], :persist => true, :data_persist => true, :unnamed => true)
    end
    CACHE[:gene_position][key]
  end

  def self.gene_strand_index(organism)
    key = organism
    CACHE[:gene_strand] ||= {}
    if CACHE[:gene_strand][key].nil?
      CACHE[:gene_strand][key] = Organism.gene_positions(organism).tsv(:fields => ["Strand"], :type => :single, :persist => true, :unnamed => true).to_hash
    end
    CACHE[:gene_strand][key]
  end

  def self.exon_position_index(organism, chromosome)
    key = [organism, chromosome]
    CACHE[:exon_position] ||= {}
    if CACHE[:exon_position][key].nil?
      CACHE[:exon_position][key] = TSV.range_index(Organism.exons(organism).produce, "Exon Chr Start", "Exon Chr End", :filters => [["field:Chromosome Name", chromosome]], :persist => true, :data_persist => true, :unnamed => true)
    end
    CACHE[:exon_position][key]
  end

  def self.exon_start_index(organism, chromosome)
    key = [organism, chromosome]
    CACHE[:exon_start] ||= {}
    if CACHE[:exon_start][key].nil?
      CACHE[:exon_start][key] = TSV.pos_index(Organism.exons(organism).produce, "Exon Chr Start", :filters => [["field:Chromosome Name", chromosome]], :persist => true, :data_persist => true, :unnamed => true)
    end
    CACHE[:exon_start][key]
  end

  def self.exon_end_index(organism, chromosome)
    key = [organism, chromosome]
    CACHE[:exon_end] ||= {}
    if CACHE[:exon_end][key].nil?
      CACHE[:exon_end][key] = TSV.pos_index(Organism.exons(organism).produce, "Exon Chr End", :filters => [["field:Chromosome Name", chromosome]], :persist => true, :data_persist => true, :unnamed => true)
    end
    CACHE[:exon_end][key]
  end


  def self.exon_info(organism)
    key = organism
    CACHE[:exon_info] ||= {}
    if CACHE[:exon_info][key].nil?
      CACHE[:exon_info][key] = Organism.exons(organism).tsv :persist => true, :serializer => :list, :unnamed => true
    end
    CACHE[:exon_info][key]
  end

  def self.exon_transcript_offsets(organism)
    key = organism
    CACHE[:exon_transcript_offsets] ||= {}
    if CACHE[:exon_transcript_offsets][key].nil?
      CACHE[:exon_transcript_offsets][key] = Organism.exon_offsets(organism).tsv :persist => true, :serializer => :double, :unnamed => true
    end
  end

  def self.transcript_sequence(organism)
    key = organism
    CACHE[:transcript_sequence] ||= {}
    if CACHE[:transcript_sequence][key].nil?
      CACHE[:transcript_sequence][key] = Organism.transcript_sequence(organism).tsv(:single, :persist => true, :unnamed => true).to_hash
    end
    CACHE[:transcript_sequence][key]
  end

  def self.transcript_5utr(organism)
    key = organism
    CACHE[:transcript_5utr] ||= {}
    if CACHE[:transcript_5utr][key].nil?
      CACHE[:transcript_5utr][key] = Organism.transcript_5utr(organism).tsv(:single, :persist => true, :unnamed => true).to_hash
    end
    CACHE[:transcript_5utr][key]
  end
   
  def self.transcript_3utr(organism)
    key = organism
    CACHE[:transcript_3utr] ||= {}
    if CACHE[:transcript_3utr][key].nil?
      CACHE[:transcript_3utr][key] = Organism.transcript_3utr(organism).tsv(:single, :persist => true, :unnamed => true).to_hash
    end
    CACHE[:transcript_3utr][key]
  end

  def self.transcript_phase(organism)
    key = organism
    CACHE[:transcript_phase] ||= {}
    if CACHE[:transcript_phase][key].nil?
      CACHE[:transcript_phase][key] = Organism.transcript_phase(organism).tsv(:single, :persist => true, :unnamed => true, :cast => nil).to_hash
    end
    CACHE[:transcript_phase][key]
  end

  def self.transcript_protein(organism)
    key = organism
    CACHE[:transcript_protein] ||= {}
    if CACHE[:transcript_protein][key].nil?
      CACHE[:transcript_protein][key] = Organism.transcripts(organism).tsv(:persist => true, :fields => ["Ensembl Protein ID"], :type => :single, :unnamed => true).to_hash
    end
    CACHE[:transcript_protein][key]
  end

  def self.transcript_info(organism)

    key = organism
    CACHE[:transcript_info] ||= {}
    if CACHE[:transcript_info][key].nil?
      CACHE[:transcript_info][key] = Persist.persist_tsv(Organism.transcripts(organism), "Transcript Info", {:organism => organism},{:persist => true, :serializer => :list}) do |data|
        tsv = Organism.transcript_5utr(organism).tsv(:type => :single, :unnamed => true, :cast => :to_i)
        tsv.attach Organism.transcript_3utr(organism).tsv(:type => :single, :unnamed => true, :cast => :to_i)
        tsv.attach Organism.transcript_phase(organism).tsv(:type => :single, :unnamed => true, :cast => :to_i)

        tsv.annotate data
        data.serializer = :integer_array
        ddd data.serializer
        ddd data.serializer_module
        tsv.each{|k,v| data[k] = v}
        data
      end
    end
    CACHE[:transcript_info][key]
  end

  def self.chromosome_file(organism, chromosome)
    Organism[File.join(organism, "chromosome_#{chromosome}")].produce.find
  end

  def self.snp_position_index(organism, chromosome)
    key = [organism, chromosome]
    CACHE[:snp_position] ||= {}
    CACHE[:germline_variations] ||= Organism.germline_variations(organism).tsv :persist => true, :unnamed => true
    if CACHE[:snp_position][key].nil?
      CACHE[:germline_variations].filter
      CACHE[:germline_variations].add_filter "field:Chromosome Name", chromosome
      CACHE[:snp_position][key] = CACHE[:germline_variations].pos_index("Chromosome Start", :persist => true, :unnamed => true) #TSV.pos_index(Organism.germline_variations(organism), "Chromosome Start", :filters => [["field:Chromosome Name", chromosome]], :persist => true, :data_persist => true, :monitor => true)
    end
    CACHE[:snp_position][key]
  end

  def self.somatic_snv_position_index(organism, chromosome)
    key = [organism, chromosome]
    CACHE[:snv_position] ||= {}
    CACHE[:germline_variations] ||= Organism.somatic_variations(organism).tsv :persist => true, :unnamed => true
    if CACHE[:snv_position][key].nil?
      CACHE[:germline_variations].filter
      CACHE[:germline_variations].add_filter "field:Chromosome Name", chromosome
      CACHE[:snv_position][key] = CACHE[:germline_variations].pos_index("Chromosome Start", :persist => true, :unnamed => true) #TSV.pos_index(Organism.germline_variations(organism), "Chromosome Start", :filters => [["field:Chromosome Name", chromosome]], :persist => true, :data_persist => true, :monitor => true)
    end
    CACHE[:snv_position][key]
  end
end
