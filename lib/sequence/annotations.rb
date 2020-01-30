module Sequence

  dep :mutated_isoforms_fast
  dep :splicing_mutations
  task :affected_genes => :tsv do
    organism = step(:mutated_isoforms_fast).inputs["organism"]

    protein_index = Organism.identifiers(organism).index :target => "Ensembl Gene ID", :fields => ["Ensembl Protein ID"], :persist => true
    transcript_index = Organism.gene_transcripts(organism).index :target => "Ensembl Gene ID", :fields => ["Ensembl Transcript ID"], :persist => true

    pasted = TSV.paste_streams([step(:splicing_mutations), step(:mutated_isoforms_fast)], :sort => true)
    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => organism
    dumper.init
    TSV.traverse pasted, :bar => "Affected Genes", :type => :array, :into => dumper do |line|
      next if line =~ /^#/
      m, *rest = line.split("\t")
      genes = rest.collect do |part|
        next if part.nil? or part.empty?
        e, c = part.split(":")
        g = if e =~ /ENSP/
              next unless c =~ /^([A-Z*])\d+([A-Z*]+)/i and $1 != $2
              protein_index[e]
            else
              next if c =~ /^UTR\d$/
              transcript_index[e]
            end
        g
      end.compact.uniq

      next if genes.empty?
      [m, genes]
    end
  end
  export_asynchronous :affected_genes

  def self.ablated_domains(mi, organism=Organism.default_code("Hsa"))
    require 'rbbt/sources/InterPro'
    @@ensp2uni ||= Organism.identifiers(organism).index :target => "UniProt/SwissProt Accession", :persist => true, :fields => ["Ensembl Protein ID"], :unnamed => true
    @@domain_info ||= InterPro.protein_domains.tsv :persist => true, :unnamed => true
    return [] unless mi =~ /:.*(\d+)(FrameShift|\*)$/
    pos = $1.to_i
    protein = mi.partition(":")[0]
    uni = @@ensp2uni[protein]
    if uni.nil?
      Log.warn "No UniProt/SwissProt accession for protein: #{ protein }" 
      return []
    end
    ablated_domains = []
    if uni
      domains = @@domain_info[uni]
      if domains
        Misc.zip_fields(domains).each do |domain,start,eend|
          if eend.to_i > pos
            ablated_domains << domain
          end
        end
      end
    end
    ablated_domains
  end

  def self.mi_sequence_ontology_term(mi, mut, organism)
    protein, change = mi.split(":")
    chr,pos,allele = mut.split(":")

    case 
    when ablated_domains(mi, organism).any?
      'transcript_ablation' 
    #when false
    #  'splice_acceptor_variant'
    #when false
    #  'splice_donor_variant'
    when change =~ /[A-Z]\d+\*$/ 
      'stop_gained'
    when change =~ /Frame/
      'frameshift_variant'
    when change =~ /\*\d+[A-Z]$/
      'stop_lost'
    when change =~ /M1[^M]$/
      'start_lost'
    #when false
    #  'transcript_amplification'
    when change =~ /Indel/ && ! allele.include?("-")
      'inframe_insertion'
    when change =~ /Indel/
      'inframe_deletion'
    when change =~ /([A-Z])\d+([A-Z])/ && $1 != $2
      'missense_variant'
    when false
      'protein_altering_variant'
    #when false
    #  'splice_region_variant'
    when change =~ /\d+X$/
      'incomplete_terminal_codon_variant'
    when change =~ /\*\d+\*/
      'stop_retained_variant'
    when change =~ /([A-Z*])\d+([A-Z*])/ && $1 == $2
      'synonymous_variant'
    when false
      'coding_sequence_variant'
    #when false
    #  'mature_miRNA_variant'
    when change == "UTR5"
      '5_prime_UTR_variant'
    when change == "UTR3"
      '3_prime_UTR_variant'
    #when false
    #  'non_coding_transcript_exon_variant'
    #when false
    #  'intron_variant'
    #when false
    #  'NMD_transcript_variant'
    #when false
    #  'non_coding_transcript_variant'
    #when false
    #  'upstream_gene_variant'
    #when false
    #  'downstream_gene_variant'
    #when false
    #  'TFBS_ablation'
    #when false
    #  'TFBS_amplification'
    #when false
    #  'TF_binding_site_variant'
    #when false
    #  'regulatory_region_ablation'
    #when false
    #  'regulatory_region_amplification'
    #when false
    #  'feature_elongation'
    #when false
    #  'regulatory_region_variant'
    #when false
    #  'feature_truncation'
    #when false
    #  'intergenic_variant'
    else
      raise "Not identified: #{ mi }"
    end
  end

  def self.mut_sequence_ontology_term(mut, juncs, genes, exons, up_genes, down_genes, organism)
    miRNAs = []
    chr,pos,allele = mut.split(":")
    @ense2enst ||= Organism.transcript_exons(organism).index :fields => ["Ensembl Exon ID"], :target => "Ensembl Transcript ID", :persist => true, :merge => true, :type => :double
    @enst2biotype ||= Organism.transcript_biotype(organism).tsv :persist => true, :type => :single
    transcripts = exons.collect{|e| @ense2enst[e]}.flatten.uniq

    case 
    #when (ad = ablated_domains(mi, organism)).any?
    #  'transcript_ablation:' 
    when juncs.select{|j| j =~ /acceptor\((\d+)\)/ and $1.to_i.abs <= 2}.any?
      'splice_acceptor_variant'
    when juncs.select{|j| j =~ /donor\((\d+)\)/ and $1.to_i.abs <= 2}.any?
      'splice_donor_variant'
    #when change =~ /[A-Z]\d+\*$/ 
    #  'stop_gained'
    #when change =~ /Frame/
    #  'frameshift_variant'
    #when change =~ /\*\d+[A-Z]$/
    #  'stop_lost'
    #when change =~ /M1[^M]$/
    #  'start_lost'
    #when false
    #  'transcript_amplification'
    #when change =~ /Indel/ && ! allele.include?("-")
    #  'inframe_insertion'
    #when change =~ /Indel/
    #  'inframe_deletion'
    #when change =~ /([A-Z])\d+([A-Z])/ && $1 != $2
    #  'missense_variant'
    #when false
    #  'protein_altering_variant'
    when juncs.any?
      'splice_region_variant'
    #when change =~ /\d+X$/
    #  'incomplete_terminal_codon_variant'
    #when change =~ /\*\d+\*/
    #  'stop_retained_variant'
    #when change =~ /([A-Z*])\d+([A-Z*])/ && $1 == $2
    #  'synonymous_variant'
    when transcripts.select{|t| @enst2biotype[t] == "protein_coding"}.any?
      'coding_sequence_variant'
    when (genes & miRNAs).any?
      'mature_miRNA_variant'
    #when change == "UTR5"
    #  '5_prime_UTR_variant'
    #when change == "UTR3"
    #  '3_prime_UTR_variant'
    when transcripts.select{|t| bt = @enst2biotype[t]; bt != "nonsense_mediated_decay" and bt != "protein_coding"}.any?
      'non_coding_transcript_exon_variant'
    when genes.any? && ! exons.any?
      'intron_variant'
    when transcripts.select{|t| @enst2biotype[t] == "nonsense_mediated_decay"}.any?
      'NMD_transcript_variant'
    when false
      'non_coding_transcript_variant'
    when up_genes.any?
      'upstream_gene_variant'
    when down_genes.any?
      'downstream_gene_variant'
    #when false
    #  'TFBS_ablation'
    #when false
    #  'TFBS_amplification'
    #when false
    #  'TF_binding_site_variant'
    #when false
    #  'regulatory_region_ablation'
    #when false
    #  'regulatory_region_amplification'
    #when false
    #  'feature_elongation'
    #when false
    #  'regulatory_region_variant'
    #when false
    #  'feature_truncation'
    when genes.empty? && up_genes.empty? && down_genes.empty?
      'intergenic_variant'
    else
      nil
    end
  end

  dep :mutated_isoforms_fast
  dep :exon_junctions, :positions => :mutations
  dep :genes, :positions => :mutations
  dep :exons, :positions => :mutations
  dep :TSS, :positions => :mutations
  dep :TES, :positions => :mutations
  task :sequence_ontology => :tsv do
    Workflow.require_workflow "InterPro"
    so_term_order = Rbbt.share.databases.sequence_ontology.terms.tsv :fields => ["Order"], :type => :single, :cast => :to_i
    organism = step(:mutated_isoforms_fast).inputs[:organism]
    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Mutated Isoform", "MI SO Terms", "MUT SO Terms", "SO Term"], :type => :double, :namespace => organism
    dumper.init
    pasted = TSV.paste_streams([step(:mutated_isoforms_fast), step(:exon_junctions), step(:genes), step(:exons), step(:TSS), step(:TES)], :fix_flat => true, :sort => true)
    TSV.traverse pasted, :into => dumper, :bar => "Sequence ontology" do |mut,values|
      mut = mut.first if Array === mut
      mis, juncs, genes, exons, up_genes, down_genes = values
      mi_so_terms = mis.collect{|mi| Sequence.mi_sequence_ontology_term(mi,mut,organism) }
      mut_so_terms = Sequence.mut_sequence_ontology_term(mut, juncs, genes, exons, up_genes, down_genes, organism)
      so_terms = mi_so_terms 
      so_terms = so_terms + [mut_so_terms] if mut_so_terms
      top_term = so_terms.sort_by{|t| so_term_order[t] || 1000}.first

      [mut,[mis, mi_so_terms, mut_so_terms, [top_term]]]
    end
  end
  export_asynchronous :sequence_ontology
end


