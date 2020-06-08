
rule treetime:
    input:
        tree = config["output_path"] + "/5/{lineage}/trees/uk_lineage_UK{i}.tree",
        metadata = config["output_path"] + "/5/cog_gisaid.lineages.with_all_traits.with_phylotype_traits.csv",
        fasta = config["output_path"] + "/3/cog_gisaid.fasta",
    params:
        lineage="{lineage}",
        i="{i}"
    output:
        temp_fasta = temp(config["output_path"] + "/6/{lineage}/trees/uk_lineage_UK{i}.temp.fasta"),
        fasta = config["output_path"] + "/6/{lineage}/trees/uk_lineage_UK{i}.fasta",
        metadata = config["output_path"] + "/6/{lineage}/trees/uk_lineage_UK{i}.timetree.csv",
        tree = config["output_path"] + "/6/{lineage}/trees/uk_lineage_UK{i}.nexus",
        directory = directory(config["output_path"] + "/6/{lineage}/trees/uk_lineage_UK{i}_timetree/")
    log:
        config["output_path"] + "/logs/6_timetree_run_lineage{lineage}_uk{i}.log"
    shell:
        """
        sed "s/'//g" {input.tree} > {output.tree}

        fastafunk extract \
          --in-fasta {input.fasta} \
          --in-tree {input.tree} \
          --out-fasta {output.fasta} &>> {log}

        fastafunk fetch \
          --in-fasta {output.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date \
          --out-fasta {output.temp_fasta} \
          --out-metadata {output.metadata} \
          --restrict &>> {log}


        set +e

        treetime \
          --aln {output.fasta} \
          --tree {output.tree} \
          --dates {output.metadata} \
          --name-column sequence_name \
          --date-column sample_date \
          --outdir {output.directory} &>> {log}

        exitcode=$?
        if [ $exitcode -eq 0 ]
        then
            exit 0
        else
            echo "timetree failed for this sample" &>> {log}
            exit 0
        fi
        """


LINEAGES = []
df = pd.read_csv(config["lineage_splits"])
for i,row in df.iterrows():
    LINEAGES.append(row['lineage'])

UK = []
for l in LINEAGES:
    UK.append(glob_wildcards(config["output_path"] + "/5/" + l + "/trees/uk_lineage_UK{i}.tree").i)

LINEAGES_REP = []
for i,x in enumerate(UK):
    LINEAGES_REP = LINEAGES_REP + [LINEAGES[i]] * len(x)

UK = [item for sublist in UK for item in sublist]

rule summarize_treetime:
    input:
        expand(config["output_path"] + "/logs/6_timetree_run_lineage{lineage}_uk{i}.log", zip, lineage = LINEAGES_REP, i = UK)
    params:
        webhook = config["webhook"],
    log:
        config["output_path"] + "/logs/6_summarize_treetime.log"
    shell:
        """
        echo '{{"text":"' > 6_data.json
        echo "*Step 6: treetime for UK lineages complete*\\n" >> 6_data.json
        echo '"}}' >> 6_data.json
        echo 'webhook {params.webhook}'
        curl -X POST -H "Content-type: application/json" -d @6_data.json {params.webhook}
        """










#