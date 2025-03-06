process COMBINE_REPORT {
    label 'process_single'
    label 'pandas'

    input:
    path(ch_statistics_out)
    path(ch_busco_out)
    path(ch_omark_out)
    path(ch_brh_out)
    path(ch_distribution_out)
    path(ch_index_stats)
    path(ch_feature_stats)
    path(ch_psauron_stats)
    path(ch_intron_stats)
    
    output:
    path("Evaluation_output.json"), emit: statistics_json
    path("Evaluation_output.txt"), emit: statistics_txt

    script:
    """
    python3 ${projectDir}/bin/combine_report.py --statistics_output ${ch_statistics_out} \\
    --busco_output ${ch_busco_out} --omark_output ${ch_omark_out} \\
    --brh_output ${ch_brh_out} --idxstats_output ${ch_index_stats} --compare_distribution_output ${ch_distribution_out} \\
    --feature_output ${ch_feature_stats} --psauron_output ${ch_psauron_stats} --intron_output ${ch_intron_stats}
    """
}
