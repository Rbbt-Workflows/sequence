module Sequence

  input :breakpoints, :array, "Breakpoints"
  input *ORGANISM_INPUT
  dep :transcript_offsets, :positions => :placeholder do |jobname,options|
    list = TSV.traverse options[:breakpoints], :type => :array, :into => [] do |breakpoint|
      res = breakpoint.split("-")
      res.extend MultipleResult
      res
    end
    options[:positions] = list
    {:jobname => jobname, :inputs => options}
  end
  task :fusions => :tsv do |breakpoints,organism|
    transcript_offsets = step(:transcript_offsets).load
    transcript_sequence ||= Sequence.transcript_sequence(organism)
    transcript_biotype ||= Sequence.transcript_biotype(organism)
    transcript_protein = Sequence.transcript_protein(organism)
    
    dumper = TSV::Dumper.new :key_field => "Breakpoint", :fields => ["Protein", "Wildtype", "Position", "Downstream"], :type => :double, :organism => organism
    dumper.init
    TSV.traverse breakpoints, :type => :array, :into => dumper do |breakpoint|
      pos5, pos3 = breakpoint.split("-")
      offs5 = transcript_offsets[pos5]
      offs3 = transcript_offsets[pos3]
      res = []
      res.extend MultipleResult
      offs5.each do |off5|
        transcript5, offset5, strand5 = off5.split(":")
        next unless transcript_biotype[transcript5].include? 'coding'

        protein5 = transcript_protein[transcript5]
        transcript_sequence5 = transcript_sequence[transcript5]
        codon5 = Sequence.codon_at_transcript_position(organism, transcript5, offset5.to_i)
        next if codon5.include? "UTR"
        offs3.each do |off3|
          transcript3, offset3, strand3 = off3.split(":")
          transcript_sequence3 = transcript_sequence[transcript3]
          codon3 = Sequence.codon_at_transcript_position(organism, transcript3, offset3.to_i)

          if strand5 == strand3
            appended = transcript_sequence3[0..offset3.to_i]
          else
            appended = transcript_sequence3[offset3.to_i..-1].reverse
          end

          downstream5, codon_num = Sequence.downstream(organism, transcript5, "+" + appended, codon5)

          res << [breakpoint, [protein5, codon_num, downstream5]]
        end
      end
      res
    end
  end
end
