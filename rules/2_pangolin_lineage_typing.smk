
#previous_stage defined but never used, commenting out for now
rule uk_normal_pangolin:
    input:
#        previous_stage = config["output_path"] + "/logs/1_summarize_preprocess_uk.log",
        fasta = rules.uk_add_dups_to_lineageless.output.fasta,
    params:
        outdir = config["output_path"] + "/2/normal_pangolin",
        tmpdir = config["output_path"] + "/2/normal_pangolin/tmp"
    output:
        lineages = config["output_path"] + "/2/normal_pangolin/lineage_report.csv"
    log:
        config["output_path"] + "/logs/2_uk_normal_pangolin.log"
    shell:
        """
        pangolin {input.fasta} \
        --outdir {params.outdir} \
        --tempdir {params.tmpdir}  >> {log} 2>&1
        """

#originally took rules.uk_add_previous_lineages_to_metadata.output.metadata from rule_1
rule uk_add_pangolin_lineages_to_metadata:
    input:
        metadata = rules.uk_add_del_finder_result_to_metadata.output.metadata,
        lineages = rules.uk_normal_pangolin.output.lineages
    output:
        metadata = config["output_path"] + "/2/uk.with_new_lineages.csv",
    log:
        config["output_path"] + "/logs/2_uk_add_normal_pangolin_lineages_to_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.lineages} \
          --index-column fasta_header \
          --join-on taxon \
          --new-columns lineage \
          --out-metadata {output.metadata} &>> {log}
        """


rule get_filled_analysis_instrument:
    input:
        metadata = rules.uk_add_pangolin_lineages_to_metadata.output.metadata
    output:
        metadata = config["output_path"] + "/2/filled_analysis_instrument.csv"
    run:
        import pandas as pd

        df = pd.read_csv(input.metadata)

        df.loc[:,'sequence_length'] = df.loc[:,'length']
        df.loc[:,'pango'] = df.loc[:,'lineage']

        df = df.loc[:,['central_id', 'redcap_repeat_instance', \
                        'consensus', 'ave_depth', 'sequence_length', \
                        'missing', 'gaps', 'pango', 'ph_cluster', \
                        'p323l', 'd614g', 'n439k', 'del_1605_3', \
                        'epi_week', 'analysis_complete']]

        df.to_csv(output.metadata, index=False)


rule import_analysis_instrument_to_redcap:
    input:
        metadata = rules.get_filled_analysis_instrument.output.metadata,
        redcap_db = config["redcap_access"]
    output:
        metadata = config["output_path"] + "/2/analysis_intrument_form_exported.csv"
    run:
        with open(input.redcap_db, 'r') as f:
            read_file = f.read().strip('\n')
            url = read_file.split(',')[0]
            key = read_file.split(',')[1]
        f.close()

        proj = redcap.Project(url, key)

        df = pd.read_csv(input.metadata)
        df.insert(1, 'redcap_repeat_instrument', 'analysis')
        df.loc[:,'analysis_complete'] = pd.to_numeric(df.loc[:,'analysis_complete'], downcast='integer')

        proj.import_records(df)
        df.to_csv(output.metadata)


#removed where-column line in fastafunk fetch
rule uk_output_lineage_table:
    input:
        fasta = rules.uk_filter_omitted_sequences.output.fasta,
        metadata = rules.uk_add_pangolin_lineages_to_metadata.output.metadata
    output:
        fasta = config["output_path"] + "/2/uk.matched.fasta",
        metadata = config["output_path"] + "/2/uk.matched.lineages.csv"
    log:
        config["output_path"] + "/logs/2_uk_output_full_lineage_table.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column fasta_header \
          --filter-column fasta_header adm0 adm1 adm2 \
                          sample_date epi_week \
                          lineage uk_lineage \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --low-memory \
          --restrict
        """


#commented out as it doesn't seem immediately useful
#had to include 'analysis_form' input to upload to redcap
#there's probably a better way of doing that
rule summarize_pangolin_lineage_typing:
    input:
        fasta = rules.uk_output_lineage_table.output.fasta,
        metadata = rules.uk_output_lineage_table.output.metadata,
        analysis_form = rules.import_analysis_instrument_to_redcap.output.metadata
    params:
        grapevine_webhook = config["grapevine_webhook"],
        json_path = config["json_path"],
        date = config["date"]
    log:
        config["output_path"] + "/logs/2_summarize_pangolin_lineage_typing.log"
#    shell:
#        """
#        echo '{{"text":"' > {params.json_path}/2_data.json
#        echo "*Step 2: {params.date} COG-UK pangolin typing complete*\\n" >> {params.json_path}/2_data.json
#        echo '"}}' >> {params.json_path}/2_data.json
#        echo "webhook {params.grapevine_webhook}"
#        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/2_data.json {params.grapevine_webhook}
#        """
