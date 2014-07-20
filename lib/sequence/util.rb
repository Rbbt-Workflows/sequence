module Sequence

  input *POSITIONS_INPUT
  task :sort => :array do |positions|
    stream = TSV.traverse positions, :type => :array, :into => :stream do |line|
      line = line.strip
      next if line.empty?
      chr, pos, *rest = line.split(":")

      c = chr.to_i
      c = chr.bytes[0].to_s if c == 0
      [("%02d" % c)<<("%020d" % pos.to_i), line] * "\t"
    end

    CMD.cmd('sort -n|cut -f 2', :in => stream, :pipe => true)
  end

  input :positions, :array, "Genomic positions"
  input :source, :string, "Original organism code and build", "Hsa/may2009"
  input :target, :string, "Target organism code and build", "Hsa/jan2013"
  task :lift_over => :array do |positions, source, target|
    raise ParameterException, "No positions given" if positions.nil?
    Organism.liftOver(positions, source, target)
  end
end
