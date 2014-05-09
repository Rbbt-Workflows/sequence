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

      return [header, next_line]
    end

    def self.parse_info_fields(info_fields, info)
      values = {}
      info.split(";").each{|p| k,v = p.split("="); v ||= true; values[k] = v}
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

    def self.open_stream(vcf, no_info = false, no_format = false)
      vcf = TSV.get_stream vcf
      header, line = header vcf

      if line =~ /#/
        fields = line.sub(/^#/,'').split(/\s+/)
        line = vcf.gets
      else
        fields = %w(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample)[0..line.split(/\s+/).length-1]
      end

      info_fields = header["INFO"].keys if header.include? "INFO"
      format_fields = header["FORMAT"].keys if header.include? "FORMAT"

      info_pos = fields.index("INFO") unless no_info
      format_pos = fields.index("FORMAT") unless no_format
      sample_fields = format_pos ? fields[format_pos+1..-1] : []

      stream_fields = ["RS ID", "Quality"]
      stream_fields.concat info_fields if info_pos
      stream_fields.concat sample_fields.collect{|s| format_fields.collect{|f| [s,f] * ":" }}.flatten if format_pos

      Misc.open_pipe false, false do |sin|
        begin
          sin.puts TSV.header_lines "Genomic Mutation", stream_fields, :type => :list
          while line

            if line !~ /\w/
              line = vcf.gets
              next
            end


            line_values = []

            chr, position, id, ref, alt, qual, filter, *rest = parts = line.split(/\s+/)
            chr.sub! 'chr', ''

            position, alt = Misc.correct_vcf_mutation(position.to_i, ref, alt)
            mutation = [chr, position.to_s, alt * ","] * ":"

            line_values << mutation
            line_values << id
            line_values << qual

            if info_pos
              info_values = parse_info_fields(info_fields, parts[info_pos])
              line_values.concat info_values
            end

            if format_pos
              format_values = parse_format_fields(format_fields, parts[format_pos], parts[format_pos+1..-1])
              line_values.concat format_values
            end

            sin.puts line_values * "\t"

            line = vcf.gets
          end
          vcf.join if vcf.respond_to? :join
          sin.close if sin.respond_to? :close
        rescue Exception
          Log.exception $!
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
  task :expanded_vcf => :tsv do |vcf,info,format|
    Sequence::VCF.open_stream(vcf, !info, !format)
  end

  dep do |name, options|
    Sequence.job(:expanded_vcf, name, options.merge({:info => false, :format => false}))
  end
  input :vcf_file, :text, "VCF file", nil, :stream => true
  input :quality, :float, "Quality threshold", 200
  task :genomic_mutations => :array do |vcf_file,quality|
    expanded_vcf = step(:expanded_vcf)
    Misc.consume_stream vcf_file, true 

    stream = TSV.traverse expanded_vcf, :key_field => "Genomic Mutation", :fields => ["Quality"], :cast => :to_f, :type => :single, :into => :stream do |mutation,qual|
      next if qual > 0 and qual > quality
      mutation
    end
  end
end
