$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), 'lib'))

require 'test/unit'
require 'rbbt-util'
require 'rbbt/workflow'

Workflow.require_workflow "Sequence"

class TestSequence < Test::Unit::TestCase
  def test_mutated_isoforms_10_000
    step = Sequence.example_step('mutated_isoforms', '10_000')
    mis = step.recursive_clean.run
    toffs = step.step(:transcript_offsets).load
    organism = step.info[:inputs][:organism]
    missing  = toffs.keys - mis.keys

    index = Organism.transcripts(organism).index :target => "Ensembl Protein ID", :source => "Ensembl Transcript ID", :persist => true
    missing.each do |m|
      proteins = toffs[m].collect do |to|
        t,o,*rest = to.split ":"
        index[t]
      end
      assert proteins.reject{|p| p.empty?}.empty? || m.split(":")[2] == "?"
    end
  end

  def test_mutated_isoforms_10_000_stream
    step = Sequence.example_step('mutated_isoforms', '10_000')
    step = step.recursive_clean.run true
    mis = TSV.open step 
    toffs = step.step(:transcript_offsets).load
    organism = step.info[:inputs][:organism]
    missing  = toffs.keys - mis.keys

    index = Organism.transcripts(organism).index :target => "Ensembl Protein ID", :source => "Ensembl Transcript ID", :persist => true
    missing.each do |m|
      proteins = toffs[m].collect do |to|
        t,o,*rest = to.split ":"
        index[t]
      end
      assert proteins.reject{|p| p.empty?}.empty? || m.split(":")[2] == "?"
    end
  end
end

