module Sequence

  input *POSITIONS_INPUT
  input *ORGANISM_INPUT
  input *VCF_INPUT
  input :full_reference_sequence, :boolean, "Consider the third field of the position an alternate allele and return the whole deleted sequence in indels", false 
  dep &VCF_CONVERTER
  task :reference => :tsv do |positions,organism,vcf,full_sequence|
    if dependencies.select{|d| d.task_name == :genomic_mutations}.any?
      positions.close if IO === positions
      positions = step(:genomic_mutations)
    end

    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Reference Allele"], :type => :single, :namespace => organism
    dumper.init
    chromosome_files = {}

    TSV.traverse positions, :bar => "Reference", :type => :array, :into => dumper do |position|
      begin
        chr, pos, alt = position.split(/[\s:\t]+/)
        next if pos.nil?
        chr.sub!(/^chr/i,'')
        chr = "MT" if chr == "M"
        file = chromosome_files[chr] ||= begin
                                           Sequence.chromosome_file(organism, chr)
                                         rescue Exception
                                           :missing
                                         end
        next if file == :missing

        pos = pos.to_i
        if full_sequence and alt["-"]
          file.seek pos - 2
          ref = file.read(alt.length + 1)
        else
          file.seek pos.to_i - 1
          ref = file.getc
        end
        [position, ref]
      rescue Exception
        Log.exception $!
        [position, "?"]
      end
    end
  end
  export_synchronous :reference

  dep :reference
  dep :exons
  task :gene_strand_reference => :tsv do 
    pasted = TSV.paste_streams [step(:reference), step(:exons)], :sort => false
    organism = step(:reference).inputs["organism"]
    exon_position = Sequence.exon_position(organism)
    dumper = TSV::Dumper.new :key_field => "Genomic Position", :fields => ["Gene Strand Reference Allele"], :type => :single, :namespace => organism
    dumper.init
    count = 0
    res = TSV.traverse pasted, :bar => "Gene strand reference", :type => :array, :into => dumper do |line|
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
    res
  end
  export_synchronous :gene_strand_reference

  dep :exons do |jobname,options|
    Sequence.job(:exons, jobname, options.merge(:positions => options[:mutations]))
  end
  input *MUTATIONS_INPUT
  task :to_watson => :array do |mutations|
    organism = step(:exons).inputs["organism"]
    exon_position = Sequence.exon_position(organism)
    pasted = TSV.paste_streams([mutations, step(:exons)], :sort => false)
    count = 0
    TSV.traverse pasted, :bar => "To watson", :type => :array, :into => :stream do |line|
      next if line =~ /#/
      count += 1
      mutation, *exons = line.split("\t")
      ex = exons.reject{|e| e.empty?}

      if ex.empty?
        mutation
      else
        ex_strands = ex.collect{|e| exon_position[e] and exon_position[e][0] }.compact
        if ex_strands.select{|s| s == -1}.any? and ex_strands.select{|s| s == 1}.empty?
          mutation
        else
          chr,pos,allele, *rest = mutation.split(":")
          allele = allele.split("").each{|base| Misc::BASE2COMPLEMENT[base] || base} * ""
          ([chr,pos,allele] + rest) * ":"
        end
      end
    end
  end
  export_synchronous :to_watson

  dep :gene_strand_reference do |jobname,options|
    Sequence.job(:gene_strand_reference, jobname, options.merge(:positions => options[:mutations]))
  end
  input *MUTATIONS_INPUT
  input *ORGANISM_INPUT
  task :is_watson => :boolean do |mutations,organism|
    Misc.consume_stream mutations


    gene_reference = step(:gene_strand_reference)
    Step.wait_for_jobs [gene_reference]

    reference = step(:gene_strand_reference).step(:reference)

    dumper = TSV.paste_streams([gene_reference, reference], :sort => true)

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

  dep :reference
  task :add_reference => :array do 
    TSV.traverse step(:reference), :type => :array, :into => :stream do |line|
      next if line =~ /^#/
      mut, ref = line.split "\t"
      parts = mut.split ":"
      [parts[0],parts[1],ref + ">" + parts[2]] * ":"
    end
  end
  export_asynchronous :add_reference
end
