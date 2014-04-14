module Sequence

  input *POSITIONS_INPUT
  input *ORGANISM_INPUT
  input *VCF_INPUT
  dep &VCF_CONVERTER
  task :reference => :tsv do |positions,organism,vcf|
    positions = step(:genomic_mutations) if step(:genomic_mutations)
    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Reference Allele"], :type => :single, :namespace => organism
    dumper.init
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
        Log.exception $!
      end
    end
  end
  export_synchronous :reference

  dep :reference
  dep :exons
  task :gene_strand_reference => :tsv do 
    organism = step(:reference).info[:inputs][:organism]
    exons = step(:exons).join.load
    exon_position = Sequence.exon_position(organism)
    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Gene Strand Reference Allele"], :type => :single, :namespace => organism
    dumper.init
    TSV.traverse step(:reference), :type => :single, :into => dumper do |position,reference|
      begin
        ex = exons[position]
        if ex.nil? or ex.empty?
          next [position, reference] 
        else
          ex_strands = ex.collect{|e| exon_position[e][0] }
          if ex_strands.select{|s| s == -1}.any? and ex_strands.select{|s| s == -1}.empty?
            [position, Misc::BASE2COMPLEMENT[reference]]
          else
            [position, reference]
          end
        end
      rescue
        Log.exception $!
      end
    end
  end
  export_synchronous :gene_strand_reference

end
