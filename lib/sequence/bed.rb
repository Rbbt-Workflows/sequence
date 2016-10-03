module Sequence

  input :bed_file, :file, "BED file"
  task :bed_index => :string do |bed_file|

    max_size = 0
    chromosome_files = {}
    TSV.traverse bed_file, :type => :array do |line|
      orig_chr, start, eend, *rest = line.split("\t")

      rest_str = rest * "\t"
      rest_str = rest_str.strip

      chr = orig_chr.sub('chr','')
      chr = "MT" if chr == "M"

      size = rest_str.length
      max_size = size if size > max_size
      chr_file = chromosome_files[chr] ||= begin
                                             FileUtils.mkdir_p file('tmp') unless File.exists? file('tmp')
                                             sout, sin = Misc.pipe 
                                             sort = CMD.cmd('sort -g -k 2 -t ":"', :in => sout, :pipe => true)
                                             t = Thread.new do
                                               Misc.consume_stream(sort, false, file('tmp/io_' + chr))
                                             end
                                             ConcurrentStream.setup(sin, :threads => [t])
                                             sin
                                           end
      chr_file.puts [start, eend, rest_str] * ":"
    end
    chromosome_files.each do |chr,file| file.close; file.join end
    file('tmp').glob('io_*').each do |io_file|
      chr = File.basename(io_file).sub('io_','')
      table = FixWidthTable.new file('index/' + chr), max_size, true
      TSV.traverse io_file, :type => :array do |line|
        start, eend, value = line.split(":")
        pos = [start.to_i, eend.to_i]
        table.add_range_point pos, value
      end
      table.read
    end

    #FileUtils.rm_rf file('tmp')

    "done"
  end

  input *POSITIONS_INPUT
  dep :bed_index
  task :match_positions_to_bed => :tsv do |mutations|

    indices = {} 
    total = 0
    step(:bed_index).file('index').glob("*").collect do |file|
      next if file == 'tmp'
      chr = File.basename(file)
      table = FixWidthTable.new file
      total += table.size
      indices[chr] = table
    end

    fields = []
    indices.values.first.value(0).split("\t").length.times do |i|
      fields << "Field #{ i }"
    end

    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => fields, :type => :double
    dumper.init
    TSV.traverse mutations, :type => :array, :into => dumper do |mutation|

      chr, pos, *rest = mutation.split(":")

      chr = chr.sub('chr','')
      chr = "MT" if chr == "M"

      index = indices[chr]
      next if index.nil?

      matches = Misc.zip_fields(index.get_range(pos.to_i).collect{|m| m.split("\t")})

      next if matches.nil? or matches.empty?
      [mutation, matches]
    end

  end

  input *POSITIONS_INPUT
  input :bed_file, :file, "BED file", nil, :stream => true
  input :sorted, :boolean, "Positions and bed file are sorted", false
  task :intersect_bed => :tsv do |positions,bed_file,sorted|
    bed_io = TSV.traverse TSV.get_stream(bed_file), :type => :array, :into => :stream do |line|
      chr, start, eend, id, *rest = line.split("\t")
      chr.sub!('chr','')
      [chr,start,eend,id] * ":"
    end
    bed_io = Misc.sort_mutation_stream(bed_io) unless sorted

    bedfile = file('bedfile')
    Open.write(bedfile, bed_io)

    position_io = sorted ? TSV.get_stream(positions) : Misc.sort_mutation_stream(TSV.get_stream(positions))

    io = Misc.open_pipe do |sin|
      Misc.intersect_streams(position_io, Open.open(bedfile), sin)
    end
  end


end
