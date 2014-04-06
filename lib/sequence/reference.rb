module Sequence

  input *POSITIONS_INPUT
  input *ORGANISM_INPUT
  input *VCF_INPUT
  dep &VCF_CONVERTER
  task :reference => :tsv do |positions,organism,vcf|
    positions = step(:genomic_mutations) if step(:genomic_mutations)
    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Reference Allele"], :type => :single, :namespace => organism
    chromosome_files = {}
    TSV.traverse positions, :type => :array, :into => dumper do |position|
      begin
        chr, pos = position.split(/[\s:\t]+/)
        next if pos.nil?
        chr.sub!(/^chr/i,'')
        file = chromosome_files[chr] ||= begin
                                           Sequence.chromosome_file(organism, chr)
                                         rescue
                                           :missing
                                         end
        next if file == :missing

        file.seek pos.to_i - 1
        ref = file.getc
        [position, ref]
      rescue
      end
    end
  end
  export_synchronous :reference
end
