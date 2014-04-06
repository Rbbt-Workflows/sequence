module Sequence
  dep :mutated_isoforms
  task :mutated_principal_isoforms => :tsv do
    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"], :type => :flat
    dumper.init
    TSV.traverse step(:mutated_isoforms), :into => dumper do |k,mis|
      next if mis.empty?
      mpi = [mis.first]
      [k, mpi]
    end
  end
end

