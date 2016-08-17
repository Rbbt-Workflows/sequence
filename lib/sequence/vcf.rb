module Sequence
  module VCF
    def self.header_lines(vcf)
      header_lines = []
      while line = vcf.gets
        if line =~ /^##/
          header_lines << line
        else
          return [header_lines, line]
        end
      end
      return [header_lines, line]
    end

    def self.header(vcf)
      lines, next_line = header_lines(vcf)

      header = {}
      lines.each do |line|
        if line =~ /^##([A-Z]+)=<ID=(.*),Number=(.*),Type=(.*),Description="(.*)">/
          field, id, number, type, description = $1, $2, $3, $4, $5
          subfield = {:numer => number, :type => type, :description => description}
          header[field] ||= {}
          header[field][id] = subfield
        else
        end
      end

      return [header, next_line, lines * "" ]
    end

    def self.parse_info_fields(info_fields, info)
      values = {}
      info.split(";").each{|p| k, _sep, v = p.strip.partition("="); v ||= "true"; values[k.strip] = v.strip}
      info_fields.collect{|f| values[f] }
    end

    def self.parse_format_fields(format_fields, format, samples)
      listed_format = format.split(":")
      sample_values = samples.collect do |sample|
        values = {}
        sample.split(':').each_with_index do |v,i|
          values[listed_format[i]] = v
        end
        values
      end

      res = []
      samples.each_with_index do |sample,i|
        format_fields.each{|f|
          res << sample_values[i][f]
        }
      end
      res
    end

    def self.save_stream(stream)
      parser = TSV::Parser.new stream, :type => :list
      fields = parser.fields

      raise "No fields in stream" if fields.nil? or fields.empty?
      # Sample fields
      tmp = fields[4..-1]
      sample_fields = {}
      sample_field_names = Set.new

      while field = tmp.shift
        break unless field =~ /(.+):(.+)/
        sample, name = $1,$2
        sample_fields[sample] ||= []
        sample_fields[sample] << field
        sample_field_names << name
      end

      if sample_fields.any? and 
          sample_fields.values.first.length > 1 and 
          sample_fields.values.collect{|l| l.sort}.uniq.length == 1


        sample_field_names = sample_field_names.to_a
        format = sample_field_names * ":"

        samples = sample_fields.keys
        all_sample_fields = sample_fields.values.flatten
        all_sample_field_pos = all_sample_fields.collect{|f| TSV.identify_field(nil, fields, f)+1 }

        has_samples = true
      else
        has_samples = false
      end

      # INFO fields
      info_fields = fields[4..-1]
      info_fields -= all_sample_fields if has_samples
      info_fields_pos = info_fields.collect{|f| TSV.identify_field(nil, fields, f)  }

      vcf_key = "CHROM"
      vcf_fields = %w(POS ID REF ALT QUAL FILTER)
      vcf_fields << "INFO" if info_fields.any?
      vcf_fields << "FORMAT" if has_samples
      vcf_fields += samples if has_samples

      dumper = TSV::Dumper.new :key_field => vcf_key, :fields => vcf_fields, :preamble => parser.preamble
      dumper.init
      TSV.traverse parser, :into => dumper, :type => :array do |mutation,parts|
        orig, id, qual, filter = parts

        chr, pos, ref, alt = orig.split ":"
        res = [pos, id, ref, alt, qual, filter]
        if info_fields.any?
          info_values = parts.values_at *info_fields_pos
          str = info_fields.dup.zip(info_values).reject{|p| p[1].nil? or p[1].empty?}.collect{|f,v| [f.gsub(/\s/,'_'), v] * "=" } * ";"
          res << str
        end

        if has_samples
          sample_data = {}
          sample_values = parts.values_at *all_sample_field_pos
          all_sample_fields.zip(sample_values).each do |sample_field,value|
            sample, field = sample_field.split ":"
            sample_data[sample] ||= {}
            sample_data[sample][field] = value
          end
          sample_columns = []
          samples.each do |sample|
            fields = sample_field_names
            sample_columns << sample_data[sample].values_at(*fields) * ":"
          end
          res << format
          res += sample_columns
        end

        [chr, res]
      end
      dumper.stream
    end

    def self.open_stream(vcf, no_info = false, no_format = false, no_preamble = false)
      vcf = TSV.get_stream vcf
      header, line, preamble = header vcf

      if line =~ /#/
        fields = line.sub(/^#/,'').split(/\s+/)
        line = vcf.gets
      else
        fields = %w(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample)[0..line.split(/\s+/).length-1]
      end

      stream_fields = ["Original", "RS ID", "Quality", "Filter"]

      if line.nil?
      else

        info_fields = header["INFO"].keys if header.include? "INFO"
        info_fields ||= []
        format_fields = header["FORMAT"].keys if header.include? "FORMAT"
        format_fields ||= line.split("\t")[8].split(":") if line.split("\t")[8]

        info_pos = fields.index("INFO") unless no_info
        format_pos = fields.index("FORMAT") unless no_format
        sample_fields = format_pos ? fields[format_pos+1..-1] : []
        if sample_fields.length > 1 and sample_fields.select{|s| s =~ /[_-]T$/}.any?
          sample_fields = sample_fields.collect{|s| s.sub(/[_-]T$/, '')}
        end

        stream_fields.concat sample_fields.collect{|s| format_fields.collect{|f| [s,f] * ":" }}.flatten if format_pos
        stream_fields.concat info_fields if info_pos
      end

      Misc.open_pipe false, false do |sin|
        begin
          preamble = nil if no_preamble 
          header_lines = TSV.header_lines "Genomic Mutation", stream_fields, :type => :list, :preamble => preamble
          sin.puts header_lines
          while line

            if line !~ /\w/
              line = vcf.gets
              next
            end

            line.chomp!

            line_values = []

            chr, position, id, ref, alt, qual, filter, *rest = parts = line.split(/\t/,-1)
            orig = [chr,position,ref,alt] * ":"

            chr.sub! 'chr', ''

            position, alt = Misc.correct_vcf_mutation(position.to_i, ref, alt)
            mutation = [chr, position.to_s, alt * ","] * ":"

            line_values << mutation
            line_values << orig
            line_values << id
            line_values << qual
            line_values << filter

            if format_pos
              format_values = parse_format_fields(format_fields, parts[format_pos], parts[format_pos+1..-1])
              line_values.concat format_values
            end

            if info_pos
              info_values = parse_info_fields(info_fields, parts[info_pos])
              line_values.concat info_values
            end

            sin.puts line_values * "\t"

            line = vcf.gets
          end
          vcf.join if vcf.respond_to? :join
          sin.close if sin.respond_to? :close
        rescue Exception
          sin.abort if sin.respond_to? :abort
          raise $!
        end
      end
    end

    def self.open(vcf)
      TSV.open(self.open_stream(vcf), :filename => TSV.get_filename(vcf))
    end
  end

  input :vcf_file, :text, "VCF file", nil, :stream => true
  input :info, :boolean, "Expand the INFO field", true
  input :format, :boolean, "Expand the sample FORMAT fields", true
  input :preamble, :boolean, "Retain file preamble (With field descriptions)", true
  task :expanded_vcf => :tsv do |vcf,info,format,preamble|
    Sequence::VCF.open_stream(vcf, !info, !format,!preamble)
  end

  dep do |name, options|
    options = Misc.add_defaults options, :info => false, :format => false, :preamble => false
    Sequence.job(:expanded_vcf, name, options)
  end
  input :vcf_file, :text, "VCF file", nil, :stream => true
  input :quality, :float, "Quality threshold", nil
  task :genomic_mutations => :array do |vcf_file,quality|
    expanded_vcf = step(:expanded_vcf)
    Misc.consume_stream vcf_file, true 

    TSV.traverse expanded_vcf, :bar => "Mutations from VCF", :key_field => "Genomic Mutation", :fields => ["Quality"], :cast => :to_f, :type => :single, :into => :stream do |mutation,qual|
      next if quality and qual > 0 and qual < quality
      mutation
    end
  end

  dep do |name, options|
    options = Misc.add_defaults options, :info => true, :format => true, :preamble => false
    Sequence.job(:expanded_vcf, name, options)
  end
  input :vcf_file, :text, "VCF file", nil, :stream => true
  input :quality, :float, "Quality threshold", nil
  task :cnvs => :array do |vcf_file,quality|
    expanded_vcf = step(:expanded_vcf).join
    Misc.consume_stream vcf_file, true 

    fields = TSV.parse_header(expanded_vcf.path).fields
    tumor_tcn = fields.select{|f| f =~ /TUMOUR:TCN/i }.first
    tumor_tcn ||= fields.select{|f| f =~ /_T:TCN/i }.first
    tumor_mcn = fields.select{|f| f =~ /TUMOUR:MCN/i }.first
    tumor_mcn ||= fields.select{|f| f =~ /_T:MCN/i }.first

    fields = fields - [tumor_tcn, tumor_mcn]

    normal_tcn = fields.select{|f| f =~ /:TCN/i}.first
    normal_mcn = fields.select{|f| f =~ /:MCN/i}.first

    select_fields = ["END", normal_tcn, normal_mcn, tumor_tcn, tumor_mcn]
    TSV.traverse expanded_vcf.path.tsv(:fields => select_fields, :unnamed => true), :bar => "CNVs from VCF", :type => :list, :into => :stream do |mutation,values|
      begin
        eend, ntcn, nmcn, ttcn, tmcn = values
        chr, start, _rest = mutation.split(":")
        next if [ntcn, nmcn] == [ttcn, tmcn]
        tumor_cnv = [chr, start, eend, [ttcn, tmcn]*"-",[ntcn, nmcn]*"-"] * ":"
        tumor_cnv
      rescue
        Log.exception $!
        next
      end
    end
  end
end
