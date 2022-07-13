import base64
from io import BytesIO
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import math, os, sys, json

def goqc( input_json ):
    """
    input_json:
    {
    'analysis_name': <ANALYSIS_NAME>
    'input_file': <FILE_NAME>
    'input_file_type': davidgo, ...
    'output_dir': <LOCAL_OUTPUT_DIR>,
    'pvalue_cutoff': cutoff for plotting (default: 0.3))
    }
    """
    print('in goqc()')
    name_label = input_json['analysis_name']
    output_dir = input_json['output_dir']
    if input_json['input_file_type'].replace('_','').lower() == 'davidgo':        
        # can input more than one GO datafile, so input_file is a list
        for fname in input_json['input_file']:
            df = pd.read_csv(fname, sep='\t')
            fdr_cutoff = -1.0*math.log(float(input_json['pvalue_cutoff']), 10)

            # gene names column
            gene_names = []
            for _gn in list(df['Genes']):
                _g = _gn.split(', ')
                gene_names.append(', '.join(_g))    
            df['gene_names'] = gene_names
                
            # log10fdr and color columns
            df['-Log10(FDR)'] = list(map(lambda x: -1.0*math.log(x,10), df['FDR']))
            df.sort_values(by=['-Log10(FDR)'], inplace=True, ascending=False)
            df['colors'] = list(map(lambda x: 'green' if ('CELLULAR_COMPONENT' in x or 'CC_DIRECT' in x) 
                                    else ('blue' if ('BIOLOGICAL_PROCESS' in x or 'BP_DIRECT' in x) 
                                          else ('red' if ('MOLECULAR_FUNCTION' in x or 'MF_DIRECT' in x)
                                                else 'purple')), df['Category']))
            
            df_CC = df[(df["colors"]=='green') & (df['-Log10(FDR)']>0)]
            df_BP = df[(df["colors"]=='blue') & (df['-Log10(FDR)']>0)]
            df_MF = df[(df["colors"]=='red') & (df['-Log10(FDR)']>0)]
            df_OT = df[(df["colors"]=='purple') & (df['-Log10(FDR)']>0)]

            # GO bar plot
            fig_go_plot = go.Figure()

            trace1 = go.Bar(
                x=df_CC['-Log10(FDR)'][0:min(5,len(df_CC))][::-1],
                y=df_CC['Term'][0:min(5,len(df_CC))][::-1],
                marker=dict(color=df_CC['colors']),
                name="cellular component", 
                hovertext=df_CC['gene_names'],
                text=df_CC['Count'],
                textposition='outside',
                orientation='h')

            trace2 = go.Bar(
                x=df_MF['-Log10(FDR)'][0:min(5,len(df_MF))][::-1],
                y=df_MF['Term'][0:min(5,len(df_MF))][::-1],
                marker=dict(color=df_MF['colors']),
                name="molecular function",
                hovertext=df_MF['gene_names'],
                text=df_MF['Count'],
                textposition='outside',    
                orientation='h')

            trace3 = go.Bar(
                x=df_BP['-Log10(FDR)'][0:min(5,len(df_BP))][::-1],
                y=df_BP['Term'][0:min(5,len(df_BP))][::-1],
                marker=dict(color=df_BP['colors']),
                name="biological process",
                hovertext=df_BP['gene_names'],
                text=df_BP['Count'],
                textposition='outside',    
                orientation='h')
            
            trace4 = go.Bar(
                x=df_OT['-Log10(FDR)'][0:min(5,len(df_OT))][::-1],
                y=df_OT['Term'][0:min(5,len(df_OT))][::-1],
                marker=dict(color=df_OT['colors']),
                name="other",
                hovertext=df_OT['gene_names'],
                text=df_OT['Count'],
                textposition='outside',    
                orientation='h')

            fig_go_plot.add_trace(trace4)
            fig_go_plot.add_trace(trace3)
            fig_go_plot.add_trace(trace2)
            fig_go_plot.add_trace(trace1)
            # significance line
            fig_go_plot.add_vline(x=-1.0*math.log(0.05,10), line_width=1) #, line_dash="dash")
            
            fig_go_plot.update_layout(title="Gene Ontology Analysis - Bar Plot of Top GO Terms",
                                      xaxis_title="-Log10(FDR)",
                                      yaxis_title="GO Term",
                                      title_x = 0.5,
                                      width=1200, height=600, autosize=False,
                                      showlegend=True,
                                      xaxis={'categoryorder':'category descending'},
                                      font=dict(size=10))

            with open(os.path.join(output_dir, '{}.goqc.barplots.html'.format(name_label)), 'a') as f:
                f.write(fig_go_plot.to_html(full_html=False, include_plotlyjs='cdn'))

    return
