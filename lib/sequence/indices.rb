require 'rbbt/sources/organism'

module Sequence

  CACHE = {}

  #{{{ REFERENCE

  def self.chromosome_file(organism, chromosome)
    Persist.memory("chromosome_file", :key => [organism,chromosome]*":", :repo => CACHE) do
      Organism[File.join(organism, "chromosome_#{chromosome}")].produce.open
    end
  end

  #{{{ FEATURES
  
  def self.gene_chromosome_index(organism, chromosome)
    Persist.memory("gene_chromosome_index", :key => [organism,chromosome]*":", :repo => CACHE) do
      TSV.range_index(Organism.gene_positions(organism).produce, "Gene Start", "Gene End", :filters => [["field:Chromosome Name", chromosome]], :persist => true, :data_persist => true, :unnamed => true)
    end
  end

  def self.exon_chromosome_index(organism, chromosome)
    Persist.memory("exon_chromosome_index", :key =>  [organism,chromosome]*":", :repo => CACHE) do
      TSV.range_index(Organism.exons(organism).produce, "Exon Chr Start", "Exon Chr End", :filters => [["field:Chromosome Name", chromosome]], :persist => true, :data_persist => true, :unnamed => true)
    end
  end

  def self.transcript_chromosome_index(organism, chromosome)
    Persist.memory("exon_chromosome_index", :key => [organism,chromosome]*":", :repo => CACHE) do
      TSV.range_index(Organism.transcripts(organism).produce, "Transcript Start (bp)", "Transcript End (bp)", :filters => [["field:Chromosome Name", chromosome]], :persist => true, :data_persist => true, :unnamed => true)
    end
  end

  def self.exon_start_index(organism, chromosome)
    Persist.memory("exon_start_index", :key => [organism,chromosome]*":", :repo => CACHE) do
      TSV.pos_index(Organism.exons(organism).produce, "Exon Chr Start", :filters => [["field:Chromosome Name", chromosome]], :persist => true, :data_persist => true, :unnamed => true)
    end
  end

  def self.exon_end_index(organism, chromosome)
    Persist.memory("exon_end_index", :key => [organism,chromosome]*":", :repo => CACHE) do
      TSV.pos_index(Organism.exons(organism).produce, "Exon Chr Start", :filters => [["field:Chromosome Name", chromosome]], :persist => true, :data_persist => true, :unnamed => true)
    end
  end

  def self.exon_end_index(organism, chromosome)
    key = [organism, chromosome]
    CACHE[:exon_end] ||= {}
    if CACHE[:exon_end][key].nil?
      CACHE[:exon_end][key] = TSV.pos_index(Organism.exons(organism).produce, "Exon Chr End", :filters => [["field:Chromosome Name", chromosome]], :persist => true, :data_persist => true, :unnamed => true)
    end
    CACHE[:exon_end][key]
  end
  def self.exon_position(organism)
    Persist.memory("exon_position", :key => organism, :repo => CACHE) do
      Organism.exons(organism).tsv :persist => true, :fields => ["Exon Strand", "Exon Chr Start", "Exon Chr End"], :cast => :to_i, :unnamed => true
    end
  end

  def self.exon_transcript_offsets(organism)
    Persist.memory("exon_transcript_offsets", :key => organism, :repo => CACHE) do
      Organism.exon_offsets(organism).tsv :persist => true, :serializer => :double, :unnamed => true
    end
  end

  def self.exon_transcripts(organism)
    Persist.memory("exon_transcripts", :key => organism, :repo => CACHE) do
      Organism.transcript_exons(organism).tsv(:key_field => "Ensembl Exon ID", :fields => ["Ensembl Transcript ID"], :type => :flat, :persist => true, :unnamed => true)
    end
  end

  #{{{ TRANSCRIPT CODING

  def self.transcript_exons(organism)
    Persist.memory("transcript_exons", :key => organism, :repo => CACHE) do
      Organism.transcript_exons(organism).tsv(:persist => true, :unnamed => true)
    end
  end

  def self.transcript_sequence(organism)
    Persist.memory("transcript_sequence", :key => organism, :repo => CACHE) do
      Organism.transcript_sequence(organism).tsv(:single, :persist => true, :unnamed => true)
    end
  end

  def self.transcript_5utr(organism)
    Persist.memory("transcript_5utr", :key => organism, :repo => CACHE) do
       Organism.transcript_5utr(organism).tsv(:single, :persist => true, :unnamed => true)
    end
  end
   
  def self.transcript_3utr(organism)
    Persist.memory("transcript_3utr", :key => organism, :repo => CACHE) do
       Organism.transcript_3utr(organism).tsv(:single, :persist => true, :unnamed => true)
    end
  end

  def self.transcript_phase(organism)
    Persist.memory("transcript_phase", :key => organism, :repo => CACHE, :persist => true) do
      Organism.transcript_phase(organism).tsv(:single, :persist => true, :unnamed => true, :cast => nil)
    end
  end

  # TRANSLATION

  def self.transcript_protein(organism)
    Persist.memory("transcript_protein", :key => organism, :repo => CACHE) do
       Organism.transcripts(organism).tsv(:persist => true, :fields => ["Ensembl Protein ID"], :type => :single, :unnamed => true)
    end
  end

end
