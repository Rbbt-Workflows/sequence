require 'rbbt/sources/organism'
Workflow.require_workflow "Appris"

require 'sequence/indices'
module Sequence
  extend Workflow

  class TranscriptError < StandardError; end

  ORGANISM_INPUT = [:organism, :string, "Organism code", Organism.default_code("Hsa")]
  MUTATIONS_INPUT = [:mutations, :array, "Mutations Chr:Position:Allele (e.g. 19:54646887:T). Extra fields are ignored but kept", nil, :stream => true]
  POSITIONS_INPUT = [:positions, :array, "Positions Chr:Position (e.g. 19:54646887). Extra fields are ignored but kept", nil, :stream => true]
  RANGES_INPUT = [:ranges, :array, "Ranges Chr:Start:End (e.g. 19:54646887:54647887)", nil, :stream => true]
  WATSON_INPUT = [:watson, :boolean, "Mutations all reported on the watson (forward) strand as opposed to the gene strand", true]
  PRINCIPAL_INPUT = [:principal, :boolean, "Consider only principal isoforms", false]
  CODING_INPUT = [:coding, :boolean, "Consider only transcripts with protein_coding biotype", false]
  NS_INPUT = [:non_synonymous, :boolean, "Report only non_synonymous variants ", false]

  VCF_INPUT = [:vcf, :boolean, "Is the input an VCF file instead of genomic mutations?", false]
  VCF_CONVERTER = Proc.new do |jobname,options|
    if options[:vcf]
      {:task => :genomic_mutations, :jobname => jobname, :inputs => options.merge({:filters => [], :vcf_file => options[:mutations] || options[:positions]})}
    else
      []
    end
  end

  export :genomic_mutations
end

require 'sequence/reference'
require 'sequence/features'
require 'sequence/coding'
require 'sequence/fast'
require 'sequence/annotations'
require 'sequence/vcf'
require 'sequence/significance'
require 'sequence/util'
require 'sequence/regulation'
require 'sequence/bed'
require 'sequence/fusion'
