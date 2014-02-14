module Sequence

  def self.liftOver_map(source, target)
    source_hg = Organism.hg_build(source)
    target_hg = Organism.hg_build(target)

    case
    when (source_hg == 'hg19' and target_hg == 'hg18')
      map_url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz" 
    when (source_hg == 'hg18' and target_hg == 'hg19')
      map_url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz" 
    else
      nil
    end
  end

  def self.liftOver_bed(source_bed, target_bed, map_url)
    TmpFile.with_file() do |map_file|

      Log.debug "Downloading map file: #{map_url}"
      
      Open.write(map_file, Open.open(map_url))

      TmpFile.with_file() do |unmapped_file|
        begin 
          FileUtils.chmod(755, Rbbt.software.opt.bin.liftOver.produce.find)
        rescue Exception
        end

        CMD.cmd("#{Rbbt.software.opt.bin.liftOver.find} '#{source_bed}' '#{map_file}' '#{target_bed}' '#{unmapped_file}'")
      end
    end
  end

  def self.swap_build(input, output, source, target, field = "Genomic Mutation")
    map_url = Sequence.liftOver_map(source, target)

    input_stream = TSV.get_stream input

    parser = TSV::Parser.new input_stream

    field_pos = parser.all_fields.index field

    line = parser.first_line

    TmpFile.with_file do |source_bed|
      Log.debug "Write to source bed file: #{ source_bed }"
      Open.write(source_bed) do |file|
        while line do
          position = line.split("\t")[field_pos]

          chr, pos = position.split(":").values_at(0,1)
          file.puts ["chr#{chr}", pos.to_i-1, pos, position] * "\t"

          line = input_stream.gets
        end
      end


      TmpFile.with_file() do |target_bed|
        Log.info "Lift over into target bed: #{ target_bed }"
        Sequence.liftOver_bed(source_bed, target_bed, map_url)

        TmpFile.with_file() do |tmpfilename|
          Log.info "Write target bed into tmpfilename as a tsv: #{ tmpfilename }"
          Open.write(tmpfilename) do |tmpfile|
            tmpfile.puts TSV.header_lines "LO_SOURCE", ["LO_TARGET"]
            Open.read(target_bed) do |line|
              chr, _pos, pos, name = line.strip.split "\t"
              chr.sub! 'chr', ''
              lifted = [chr, pos, name.split(":")[2..-1]].flatten * ":"
              tmpfile.puts [name, lifted] * "\t"
            end
          end

          Log.info "Use the tmpfilename to attach the new ids"
          num_fields = parser.all_fields.length
          add_key_fields = ([field_pos] + (0..field_pos-1).to_a + (field_pos+1..num_fields-1).to_a).flatten
          clean_fields = ((1..field_pos-1).to_a + [num_fields] + (field_pos..num_fields-1).to_a).flatten
          clean_fields = ((1..field_pos).to_a + [num_fields] + (field_pos+1..num_fields-1).to_a).flatten

          merge_source = TSV.reorder_stream(TSV.get_stream(input_stream), add_key_fields)
          merge_ids = Open.open(tmpfilename)

          Log.debug "Fix columns"
          sout = Misc.open_pipe do |sin|
            TSV.merge_different_fields(merge_source, merge_ids, sin)
          end
          Misc.sensiblewrite(output) do |file|
            io = TSV.reorder_stream(sout, clean_fields)
            file.puts "#: :type=:list#namespace=Hsa/may2009"
            while line = io.gets
              next if line =~ /#:/
              line.sub!(/LO_TARGET/, "Genomic Mutation")
              file.puts line
            end
          end
        end
      end
    end
  end

  input :positions, :array, "Position", nil
  input :source, :select, "Source build", "Hsa/may2009", :select_options => ["Hsa/may2009", "Hsa/jun2013"]
  input :target, :select, "Target build", "Hsa/jan2013", :select_options => ["Hsa/may2009", "Hsa/jun2013"]
  input :map_url, :string, "Map url", "auto"
  task :liftover => :tsv do |positions,source,target,map_url|

    map_url = Sequence.liftOver_map(source, target) if map_url.nil? or map_url.empty? or map_url == "auto"
    raise RbbtParamterException.new("No map between #{ source } and #{ target }") if map_url.nil?

    tsv = TSV.setup({}, :key_field => source, :fields => [target], :type => :single, :unnamed => true)

    TmpFile.with_file() do |source_bed|

      Log.debug "Writing source_bed for Lifover"
      Open.write(source_bed) do |file|
        positions.each{|position| 
          chr, pos = position.split(":").values_at(0,1)
          file.puts ["chr#{chr}", pos.to_i-1, pos, position] * "\t"
        }
        file.puts
      end

      TmpFile.with_file() do |target_bed|
        Sequence.liftOver_bed(source_bed, target_bed, map_url)

        positions = {}
        Open.read(target_bed) do |line|
          chr, _pos, pos, name = line.strip.split "\t"
          chr.sub! 'chr', ''
          tsv[name] = [chr, pos, name.split(":")[2..-1]].flatten * ":"
        end
      end
    end

    tsv
  end

  export_synchronous :liftover
end
