require 'rbbt'
require 'rbbt/workflow'

require 'sequence/indices'
module Sequence
  extend Workflow

  class TranscriptError < StandardError; end

  ORGANISM_INPUT = [:organism, :string, "Organism code", "Hsa"]
  MUTATIONS_INPUT = [:mutations, :array, "chr:pos:mutation", nil, :stream => true]
  POSITIONS_INPUT = [:positions, :array, "Positions Chr:Position (e.g. 19:54646887). Separator can be ':', space or tab. Extra fields are ignored", nil, :stream => true]
  WATSON_INPUT = [:watson, :boolean, "Mutations all reported on the watson (forward) strand as opposed to the gene strand", true]

  VCF_INPUT = [:vcf, :boolean, "Is the input an VCF file instead of genomic mutations", false]
  VCF_CONVERTER = Proc.new do |jobname,options|
    if options[:vcf]
      Sequence.job(:genomic_mutations, jobname, options.merge(:vcf_file => options[:mutations] || options[:positions])).run false
    end
  end
end

require 'sequence/reference'
require 'sequence/features'
require 'sequence/coding'
require 'sequence/annotations'
require 'sequence/vcf'
