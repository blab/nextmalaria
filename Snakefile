'''
Quirks:
- had to give root metadata same columns as other metadata
Extra:
- clades
- hide root
- reorder colors
- Change exclusion criteria
- Change calling a snp criteria
'''
configfile: "config/parameters.yaml"
files = config["files"]

rule all:
    input:
        auspice_json = expand("auspice/Pf_{gene}.json", gene=config["gene"]),
        tip_frequencies_json = expand("auspice/Pf_{gene}_tip-frequencies.json", gene=config["gene"])

# ulimit -n 8000
# export AUGUR_RECURSION_LIMIT=10000

rule normalize:
    message: "Normalizing vcf file, splitting multiallelic variants into lines"
    input:
        vcf = files["input_vcf"],
        ref = files["reference"]
    params:
        chromosome = config["chromosome"],
        start = config["gene_start"],
        end = config["gene_end"]
    output:
        vcf = files["norm_vcf"]
    shell:
        """
        bcftools norm -f {input.ref} -m- -r {params.chromosome}:{params.start}-{params.end} -o {output.vcf} -O z {input.vcf}
        bcftools index --tbi {output.vcf}
        """

rule split_vcf:
    message: "Splitting VCF with all samples into individual samples"
    input:
        vcf = rules.normalize.output.vcf
    params:
        chromosome = config["chromosome"],
        start = config["gene_start"],
        end = config["gene_end"]
    output:
        vcf_dir = directory("data/vcfs")
    shell:
        """
        bcftools +split {input.vcf} -Oz -o {output.vcf_dir} -i "FORMAT/AD[0:1] > 0" -r {params.chromosome}:{params.start}-{params.end}
        """

rule normalize_samples:
    message: "Normalizing samples, joining multiallelic variants into a single line"
    input:
        vcf = "data/vcfs/{sample}.vcf.gz",
        ref = files["reference"]
    params:
        chromosome = config["chromosome"],
        start = config["gene_start"],
        end = config["gene_end"]
    output:
        vcf = "data/vcfs_norm/{sample}.vcf.gz"
    shell:
        """
        bcftools index --tbi -f {input.vcf}
        bcftools norm -f {input.ref} -m+ -r {params.chromosome}:{params.start}-{params.end} -Oz -o {output.vcf} {input.vcf}
        """

rule generate_exclude:
 # Note: Uses very conservative definition of SNV, there can be zero reads mapping to reference at that position. This could probably be updated to be less conservatve.
    message:
        """
        Generating list of samples to exclude:
            - multiallelic samples
            - samples where reference & alt present
        """
    input:
        vcfs = expand("data/vcfs_norm/{sample}.vcf.gz", sample=config["samples"])
    output:
        exclude = "config/exclude.txt"
    shell:
        """
        for f in {input.vcfs}; do
            QUERY=$(bcftools query -f '[%POS]' -i "N_ALT>1 || FORMAT/AD[0:0]>0" $f)
            if [[ $QUERY ]]; then
                FNAME=${{f##*/}}
                echo ${{FNAME%%.*}} >> {output.exclude}
            fi
        done
        """

rule consensus:
    message: "Calling consensus chromosome from vcf"
    input:
        vcfs = expand("data/vcfs/{sample}.vcf.gz", sample=config["samples"]),
        index = expand("data/vcfs/{sample}.vcf.gz.tbi", sample=config["samples"]),
        reference = files["reference"]
    params:
        chromosome = config["chromosome"],
        start = config["gene_start"],
        end = config["gene_end"]
    output:
        fasta = expand("results/{chromosome}_{{gene}}.fasta", chromosome=config["chromosome"])
    shell:
        """
        for f in {input.vcfs}; do
            FNAME=${{f##*/}}
            SAMPLE=${{FNAME%%.*}}
            samtools faidx {input.reference} {params.chromosome}:{params.start}-{params.end} | \
            bcftools consensus $f | \
            sed "0,/{params.chromosome}:{params.start}-{params.end}/s//$SAMPLE/" >> \
            {output.fasta}
        done
        """

rule reverse_complement:
    message: "Reverse complementing fasta"
    input:
        fasta = expand("results/{chromosome}_{{gene}}.fasta", chromosome=config["chromosome"])
    output:
        fasta = expand("results/RC_{chromosome}_{{gene}}.fasta", chromosome=config["chromosome"])
    shell:
        """
        seqkit seq {input.fasta} -rp > {output.fasta}
        """

rule filter:
    message:
        """
        Filtering to
          - {params.max_sequences} sequences from {params.group_by}
          - only strains where {params.query}
          - not in {input.exclude}
        """
    input:
        sequences = expand("results/RC_{chromosome}_{{gene}}.fasta", chromosome=config["chromosome"]),
        metadata = "data/Pf_6_metadata.tsv",
        exclude = "config/exclude.txt"
    output:
        sequences = "results/filtered_{gene}.fasta",
        metadata = "results/metadata_filtered_{gene}.tsv"
    params:
        group_by = "year country",
        max_sequences = 300,
        min_length = 2000,
        query =  "--query '(QC_pass==True) & (Is_returning_traveller==False)'"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --output {output.sequences} \
            --output-metadata {output.metadata} \
            --group-by {params.group_by} \
            --subsample-max-sequences {params.max_sequences} \
            {params.query} \
            --exclude {input.exclude} \
            --min-length {params.min_length}
        """

rule generate_gene_reference:
    message: "Extracting {params.gene} reference from reference {params.chromosome}"
    input:
        chromosome = files["reference"]
    params:
        chromosome = config["chromosome"],
        gene = config["gene"],
        start = config["gene_start"],
        end = config["gene_end"]
    output:
        gene = expand("config/{chromosome}_{{gene}}.fasta", chromosome=config["chromosome"])
    shell:
        """
        samtools faidx {input.chromosome} {params.chromosome}:{params.start}-{params.end} -i --mark-strand no > {output.gene}
        """

rule add_root:
    message: "Combining  filtered sequences & root"
    input:
        sequences = "results/filtered_{gene}.fasta",
        metadata = "results/metadata_filtered_{gene}.tsv",
        root = files["root_fasta"],
        root_metadata = files["root_metadata"]
    output:
        sequences = "results/prealigned_{gene}.fasta",
        metadata = "results/metadata_adjusted_{gene}.tsv"
    shell:
        """
        cat {input.sequences} > {output.sequences}
        cat {input.root} >> {output.sequences}
        cat {input.metadata} > {output.metadata}
        cat {input.root_metadata} >> {output.metadata}
        """

rule align:
    message: "Aligning genomes"
    input:
        sequences = "results/prealigned_{gene}.fasta",
        reference =  expand("config/{chromosome}_{{gene}}.fasta", chromosome=config["chromosome"])
    threads: 4
    output:
        alignment = "results/aligned_{gene}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --output {output.alignment} \
            --nthreads {threads} \
            --reference-sequence {input.reference} \
            --remove-reference
        """

rule tree:
    message: "Building tree"
    input:
        sequences = rules.align.output.alignment
    threads: 4
    output:
        tree = "results/tree_raw_{gene}.nwk"
    shell:
        """
        augur tree \
            --alignment {input.sequences} \
            --output {output.tree} \
            --nthreads {threads}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - use {params.clock_rate} substitution rate & {params.clock_std_dev} std. dev
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        metadata = "results/metadata_adjusted_{gene}.tsv"
    output:
        tree = "results/tree_{gene}.nwk",
        node_data = "results/branch_lengths_{gene}.json"
    params:
        coalescent = "skyline",
        date_inference = "marginal",
        clock_rate = 0.00000004, # pulled from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5389722/
        clock_std_dev = 0.00000005, # pulled from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5389722/
        root = config["root"]
        #clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root}
        """

        #--timetree \
        #--coalescent {params.coalescent} \
        #--date-confidence \
        #--date-inference {params.date_inference}  \
        #--clock-rate {params.clock_rate} \
        #--clock-std-dev {params.clock_std_dev} \
        #--clock-filter-iqd {params.clock_filter_iqd}

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment
    output:
        node_data = "results/nt_muts_{gene}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """#

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = expand("config/{chromosome}_{{gene}}.gb", chromosome=config["chromosome"])
    output:
        node_data = "results/aa_muts_{gene}.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """#

# rule clades: TBA

#rule traits: Don't want to include at this time.
#    message:
#        """
#        Inferring ancestral traits for {params.columns!s}
#          - increase uncertainty of reconstruction by {params.sampling_bias_correction} to partially account for sampling bias
#        """
#    input:
#        tree = rules.refine.output.tree,
#        metadata = rules.parse.output.metadata
#    output:
#        node_data = "results/traits.json",
#    params:
#        columns = "region country",
#        sampling_bias_correction = 3
#    shell:
#        """
#        augur traits \
#            --tree {input.tree} \
#            --metadata {input.metadata} \
#            --output {output.node_data} \
#            --columns {params.columns} \
#            --confidence \
#            --sampling-bias-correction {params.sampling_bias_correction}
#        """#

rule frequencies:
    message: "Calculating tip frequencies"
    input:
        tree = rules.refine.output.tree,
        metadata = "results/metadata_adjusted_{gene}.tsv"
    params:
        min_date = "2001-12-31",
        max_date = "2015-12-31",
        pivot_interval = 12,
        pivot_interval_units = "months"
    output:
        tip_frequencies = "auspice/Pf_{gene}_tip-frequencies.json"
    shell:
        """
        augur frequencies \
            --method kde \
            --metadata {input.metadata} \
            --tree {input.tree} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --pivot-interval {params.pivot_interval} \
            --pivot-interval-units {params.pivot_interval_units} \
            --output {output.tip_frequencies}
        """

rule colors:
    message: "Constructing colors file"
    input:
        ordering = files["color_ordering"],
        color_schemes = files["color_schemes"],
        metadata = "results/metadata_adjusted_{gene}.tsv"
    output:
        colors = "results/colors_{gene}.tsv"
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors} \
            --metadata {input.metadata}
        """

rule clades:
    message: "Assigning clades"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        clades = files["clades"]
    output:
        clade_data = "results/clades_{gene}.json"
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.nt_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = "results/metadata_adjusted_{gene}.tsv",
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        clades = rules.clades.output.clade_data,
        lat_longs = files["lat_longs"],
        colors = rules.colors.output.colors,
        auspice_config = files["auspice_config"],
        description = files["description"]
    params:
        title = '"Plasmodium falciparum Kelch-13 phylogenetic analysis"'
    output:
        auspice_json = "results/Pf_{gene}.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.clades} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --title {params.title} \
            --description {input.description} \
            --include-root-sequence \
            --output {output.auspice_json}
        """

rule hide_root:
	message: "Hiding root"
	input:
		auspice_json = rules.export.output.auspice_json

	output:
		auspice_json = "auspice/Pf_{gene}.json"
	shell:
		"""
		python scripts/hide_root.py \
			--input {input.auspice_json} \
			--output {output.auspice_json}
		"""
