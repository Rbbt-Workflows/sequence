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
    pasted = TSV.paste_streams [step(:reference), step(:exons)], :sort => true
    organism = step(:reference).info[:inputs][:organism]
    exon_position = Sequence.exon_position(organism)
    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Gene Strand Reference Allele"], :type => :single, :namespace => organism
    dumper.init
    count = 0
    TSV.traverse pasted, :type => :array, :into => dumper do |line|
      next if line =~ /#/
      count += 1
      position, reference, *exons = line.split("\t")
      ex = exons.reject{|e| e.empty?}

      if ex.empty?
        [position, reference] 
      else
        ex_strands = ex.collect{|e| exon_position[e] and exon_position[e][0] }.compact
        if ex_strands.select{|s| s == -1}.any? and ex_strands.select{|s| s == 1}.empty?
          [position, Misc::BASE2COMPLEMENT[reference]]
        else
          [position, reference]
        end
      end
    end
  end
  export_synchronous :gene_strand_reference

  dep do |jobname,options|
    Sequence.job(:gene_strand_reference, jobname, options.merge(:positions => options[:mutations]))
  end
  input *MUTATIONS_INPUT
  input *ORGANISM_INPUT
  task :is_watson => :boolean do |mutations,organism|
    Misc.consume_stream mutations

    gene_reference = step(:gene_strand_reference)
    reference = step(:gene_strand_reference).step(:reference)

    reference.join
    gene_reference.join

    dumper = TSV.paste_streams([gene_reference, reference], :type => :list, :sort => true)

    gene = 0
    watson = 0

    TSV.traverse dumper do |mutation,values|
      generef, ref = values
      next if generef == ref
      allele = mutation.split(":")[2]
      gene += 1 if generef == allele
      watson += 1 if ref == allele
    end

    gene >= watson
  end
  export_synchronous :is_watson
end
