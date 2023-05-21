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
                                             FileUtils.mkdir_p file('tmp') unless File.exist? file('tmp')
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

  input :bed_file, :file, "BED file", nil, :stream => true
  input :sorted, :boolean, "Positions and bed file are sorted", false
  task :prepare_bed_file => :text do |bed_file,sorted|
    bed_io = TSV.traverse TSV.get_stream(bed_file), :type => :array, :into => :stream do |line|
      chr, start, eend, id, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSize, blockStarts = line.split("\t")
      chr.sub!('chr','')
      [chr,start,eend,id] * ":" + "-###-" + [blockCount, blockSize, blockStarts].compact * ":"
    end
    bed_io = Misc.sort_mutation_stream(bed_io) unless sorted
    bed_io
  end

  dep :prepare_bed_file, :compute => :produce do |jobname,options|
    name = File.basename(options[:bed_file] || options["bed_file"])
    Sequence.job(:prepare_bed_file, name, options)
  end
  input *POSITIONS_INPUT
  input :sorted, :boolean, "Positions are sorted", false
  input :subset_blocks, :boolean, "Subset only matching blocks", true
  task :intersect_bed => :tsv do |positions,sorted,subset_blocks|
    position_io = sorted ? TSV.get_stream(positions) : Misc.sort_mutation_stream(TSV.get_stream(positions))

    io = Misc.open_pipe do |sin|
      sin << "#Genomic Position\tRegion" << "\n"
      Misc.intersect_streams(position_io, step(:prepare_bed_file).path.open, sin)
    end
    
    if subset_blocks
      TSV.traverse io, :type => :array, :into => :stream, :bar => "Intersecting with BED file" do |line|
        next line if line =~ /^#/
        mutation, bed_info = line.split("\t")
        mchr, mpos, *rest = mutation.split(":")

        entity_part, block_part = bed_info.split("-###-")
        bchr, bstart, beend, bid = entity_part.split(":") 
        blockCount, blockSize, blockStarts = block_part.split(":")

        if blockCount.to_i > 0
          match = false
          bstart = bstart.to_i
          mpos = mpos.to_i
          blockStarts.split(",").zip(blockSize.split(",")).each do |start, size|
            bbstart = bstart + start.to_i
            bbeend = bstart + start.to_i + size.to_i
            if mpos >= bbstart and mpos <= bbeend
              match = true 
            end
          end
          next unless match
          line.split("-###-").first
        else
          line.split("-###-").first
        end
      end
    else
      io
    end
  end


end
