"""
Inspiration from:
https://levelup.gitconnected.com/python-classes-to-standardize-plotly-figure-formatting-123fe35c8d2d
"""

import sys
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

class GenPlot:

    def __init__(self,rows=1,cols=1,shared_xaxes=True,shared_yaxes=True,subplot_titles='',legend=dict(x=0.65,y=1,font=dict(size=16))):

        # misc. figure parameters
        self.params = {'linewidth': 6,
                        'mrkrsize': 10,
                        'opacity': 1.0,
                        'width': 850,
                        'length': 700
                        }

        # font for figure labels and legend
        self.lab_dict = dict(family='Latex',
                             size=30,
                             color='black'
                             )

        # font for number labeling on axes
        self.tick_dict = dict(family='Latex',
                              size=28,
                              color='black',
                              )

        # parameters for annotations
        self.annotationFont = dict(family='Latex',
                                   size=22,
                                   color='black',
                                   )

        # line parameters
        self.line = dict(width=10
                         )

        # initialize figure as subplots
        self.fig = make_subplots(rows=rows,
                                 cols=cols,
                                 shared_xaxes=shared_xaxes,
                                 shared_yaxes=shared_yaxes,
                                 vertical_spacing=0.04,
                                 horizontal_spacing=0.04,
                                 subplot_titles=subplot_titles,
                                 )

        # set font, borders, size, background color,
        # and legend  position for figure
        self.fig.update_layout(font=self.lab_dict,
                                margin=dict(r=20,t=20,b=10),
                                autosize=False,
                                width=850,
                                height=700,
                                plot_bgcolor='white',
                                legend=legend
                                )

        # color dictionary for violin plots
        self.colors_violin = dict(MC3 = '#116C6E',
                                      MC3H = '#8F011B',
                                      # MC3H_10p = '#4A8515',
                                      CHOL = '#e28743',
                                      )

        # SLD/reflectometry colours

        # light colours: blue, green, red, brown
        self.col_light = ['#1e91ff','#6EC531','#fe0000','#fbc48d','#fc8eac']

        # dark colours: blue, green, red, brown
        self.col_dark = ['#00008a','#00693C','#8b0000','#8a4704','#c154c1']

        # size of the experimental datapoints
        self.marker_size = 20

        # marker opacity
        self.opacity = 1.0

        # symbol for the experimental data points
        self.marker_symbol =  ['circle', 'square', 'triangle-up', 'diamond', 'circle-open']

        # color dictionary for loglog plot types
        self.color_dict_loglog = dict(DLMC3 = '#32BE25', # dark green for mc3h
                                    feDLMC3 = '#32BE25',
                                    MC3 = '#32BE25',
                                    feCHL1 = '#e28743', # dark orange
                                    CHL1 = '#e28743',
                                    CHOL = '#e28743',
                                    ADE = '#cc0000', # dark red
                                    TIP3 = '#1e81b0', # dark blue
                                    )

    def show_figure(self):
        self.fig.show(renderer='png')


    def add_annotation(self,annotationText='',xshift=0,yshift=0,row=1,col=1):

        self.fig.add_annotation(x=0.36,
                               y=0.2,
                               text=annotationText,
                               showarrow=False,
                               arrowhead=5,
                               font=self.annotationFont,
                               align='right',
                               xshift=xshift,
                               yshift=yshift,
                               row=row,
                               col=col,
                               )



    def update_xaxis(self,xlabel='',showticklabels=True,axisType='linear',row=1,col=1):

        self.fig.update_xaxes(title_text='<b>'+xlabel+'</b>',
                              #range=[-1,2],
                              showline=True,
                              linecolor= 'black',
                              linewidth=5,
                              showticklabels=showticklabels,
                              ticks= 'inside',
                              mirror='allticks',
                              tickwidth=5,
                              tickcolor='black',
                              tickfont=self.tick_dict,
                              # nticks=5,
                              type=axisType,
                              exponentformat='power',
                              title_standoff=40,
                              row=row,
                              col=col
                              )


    def update_yaxis(self,ylabel='',showticklabels=True,axisType='linear',row=1,col=1):

        self.fig.update_yaxes(title_text='<b>'+ylabel+'</b>',
                              range=[0,30], # apm=[0,150]; thick=[0,30]
                              showline= True,
                              linecolor= 'black',
                              linewidth=5,
                              showticklabels=showticklabels,
                              ticks= 'inside',
                              mirror='allticks',
                              tickwidth=5,
                              tickcolor='black',
                              tickfont=self.tick_dict,
                              # nticks=10,
                              type=axisType,
                              exponentformat='power',
                              title_standoff=40,
                              row=row,
                              col=col
                              )


    def update_layout(self,axis_template,showlegend=False):

        self.fig.update_layout(
            xaxis = axis_template,
            yaxis = axis_template,
            showlegend = False,
            width = 700,
            height = 700,
            autosize = False
            )

        self.fig.update_traces(showscale=False)

        self.fig.update_xaxes(ticks='',)
        self.fig.update_yaxes(ticks='',)




    def plotViolin(self,df,lipids,cholName,row=1,col=1,showlegend=False):

        # I don't want to do this but need a quick fix.
        # re-populating lipids with all possible components
        # to plot empty dataset so they're all the same width
        lipids = ['MC3H','MC3','CHOL']

        for lipid in lipids:

            # if dataset doesn't exist, an column with out of range values is created
            try: _ = df.loc[:,lipid]
            except KeyError: df[lipid] = -1 #np.nan

            self.fig.add_trace(go.Violin(
                                    y=df.loc[:,lipid],
                                    line_color='black',
                                    fillcolor=self.colors_violin.get(lipid),
                                    showlegend=showlegend,
                                    legendgroup=lipid,
                                    scalegroup=lipid,
                                    name=lipid,
                                    opacity=0.6,
                                    ),
                                row=row,
                                col=col,
                                )

        self.fig.update_layout(violinmode='overlay') # overlay / group
        self.fig.update_traces(box_visible=True,
                               meanline_visible=True)

        self.fig.update_xaxes(type='category',
                              title_text=cholName)
        # self.fig.update_yaxes(nticks=8)





    def plotEnrich(self,enrich_list,fname,lipids,row=1,col=1,showlegend=True):

        self.fig.update_layout(
                margin=dict(l=20,r=20,t=20,b=20),
                autosize=False,
                width=700,
                height=600,
                plot_bgcolor='white',
                showlegend = True,
                legend=dict(
                            x=0.7,
                            y=1.0,
                            bgcolor='rgba(0,0,0,0)'
                            ),
                legend_font_size=20,
                )

        # initialise x and y data structures
        x, y = [], pd.DataFrame()

        # iterate through each of the enrichment dataframes
        for rep, df in enumerate(enrich_list):

            # for each df filter the needed rows via the label column into another df
            for lipid in lipids:
                _ = df.loc[df['Label']==lipid]

                # iterate through column names to get the data series
                for LIPID in lipids:

                    # store frames - frame0 as x values
                    x.append( _.loc[:,'Frame'] - _['Frame'].iloc[0] )

                    # create enrichment column name
                    column = f'{rep+1}_{lipid}-{LIPID}'

                    # convert series to list for functional data processing
                    series = list(_.loc[:,'fe'+LIPID])

                    # store in df called y
                    y[column] = series

        # initialise data structures
        Y, _, done_list = pd.DataFrame(), pd.DataFrame(), []

        # create a list of the columns in the y df
        columns = y.columns.values.tolist()

        # iterate through the columns in the list
        for column in columns:

            # get the repeat number (sample number)
            rep = column.split('_')[0]

            # get the name of the system, e.g. CHOL with respect to CHOL
            name = column.split('_')[1]

            # if the system name has not been stored in done_list, perform average
            # this prevents repeat data processing
            if name not in done_list:

                # add name to list of analyses completed 
                done_list.append(name)

                # define list of column names that have the same name (different repeat number)
                headers = [col for col in columns if col.split('_')[1] == name]

                # slice the y df via the headers list
                df = y[headers]

                # calculate the average and standard deviation
                y_av, y_std = df.mean(axis=1), df.std(axis=1)

                # add trace to figure
                self.fig.add_trace(go.Scatter(
                                    x=x[0],
                                    y=y_av,
                                    error_y=dict(
                                            type='data',
                                            array=y_std,
                                            visible=True),
                                    showlegend=True,
                                    name=name,
                                    mode='lines+markers',
                                    line=dict(width=3,
                                              color=self.colors_violin.get(name.split('-')[0]),
                                              dash='solid' if name.split('-')[0] == name.split('-')[1] else 'dot'
                                              ),
                                    marker=dict(size=10,
                                                symbol='square' if name.split('-')[0] == name.split('-')[1] else 'circle',
                                                color=self.colors_violin.get(name.split('-')[0]),
                                                line=dict(
                                                    color='black',
                                                    width=0.5,
                                                    ),
                                                ),
                                    opacity=0.8,
                                    ),
                                row=row,
                                col=1, # can't be 'col' as this is used above
                                )

        # update axis parameters
        self.fig.update_yaxes(range=[0.5,2.0],title_standoff=20)
        self.fig.update_xaxes(title_standoff=20)




    def plotScatter(self,X,Y,sysName,row=1,col=1,showlegend=True):

        # green for mc3 tail 1; red for sn2
        colour='#32BE25' if 'sn1' in sysName.split('_') else '#cc0000'

        self.fig.add_trace(go.Scatter(
                            x=X,
                            y=Y,
                            showlegend=showlegend,
                            name=str(sysName),
                            mode='lines+markers',
                            marker_symbol='square',
                            marker_color=colour,
                            line=dict(width=5,
                                      color=colour),
                            marker=dict(size=15),
                            opacity=self.params['opacity'],
                            ),
                        row=row,
                        col=col
                        )


    def plotHeatMap(self,X,Y,Z,row=1,col=1):

        max_z = np.max(Z)
        Z = Z / max_z
        # print(max_z)
        # print(Z)

        self.fig.add_trace(go.Heatmap(
                                    x=X,
                                    y=Y,
                                    z=Z,
                                    # zsmooth='best',
                                    ),
                            row=row,
                            col=col
                            )
        self.fig.update_layout(xaxis_nticks=len(X))
        self.fig.update_layout(yaxis_nticks=len(Y))

        self.fig.update_xaxes(type='category')#, tickangle=45,tickfont=dict(size=12))
        self.fig.update_yaxes(type='category')



    def plotSLD(self,SLD_x,SLD_y,contrast,nFrames,col_orig,row=1,col=1,showlegend=True):

        self.fig.update_xaxes(nticks=5)
        self.fig.update_yaxes(nticks=5)

        #self.fig.update_xaxes(range=[0,0.35])
        self.fig.update_yaxes(range=[-0.2,8])
        #self.fig.update_yaxes(range=[-8,8])

        self.fig.update_layout(
                margin=dict(l=20,r=20,t=20,b=20),
                autosize=False,
                width=700,
                height=600,
                plot_bgcolor='white',
                showlegend = True,
                legend=dict(
                            x=0.0,
                            y=1.0,
                            bgcolor='rgba(0,0,0,0)'
                            ),
                legend_font_size=20,
                )

        self.fig.add_trace(go.Scatter(
                            x=SLD_x,
                            y=SLD_y,
                            showlegend=showlegend,
                            name='%s; %d frames' %(contrast,nFrames),
                            mode='lines',
                            line=dict(width=5,
                                      color=self.col_light[col_orig],
                                      ),
                            marker=dict(size=15,
                                        ),
                            opacity=self.params['opacity'],
                            ),
                        row=row,
                        col=col
                        )



    def plotR(self,contrast,nFrames,q,R,Q,expNR,R_err_exp,labels,col_orig,row=1,col=1,showlegend=True):

        self.fig.update_xaxes(nticks=10)
        self.fig.update_yaxes(nticks=5)
        self.fig.update_xaxes(range=[0,0.35])
        self.fig.update_yaxes(range=[-6.2,0.2])

        self.fig.update_yaxes(title_standoff=10)

        self.fig.update_layout(
                margin=dict(l=20,r=20,t=20,b=20),
                autosize=False,
                width=700,
                height=600,
                plot_bgcolor='white',
                showlegend = True,
                legend=dict(
                            x=0.4,
                            y=1,
                            bgcolor='rgba(0,0,0,0)',
                            title='25 % chol. pH 3.0<br>CG: N = 338 PolyA 10nt',
                            ),
                legend_font_size=20,
                )


         # experimental data
        self.fig.add_trace(go.Scatter(
                            x=Q,
                            y=expNR,
                            error_y=dict(
                                    type='data', # value of error bar given in data coordinates
                                    array=R_err_exp,
                                    visible=True),
                            showlegend=showlegend,
                            name='%s' %(labels),
                            mode='markers',
                            marker=dict(size=10,
                                        symbol=self.marker_symbol[col_orig],
                                        color=self.col_dark[col_orig],
                                        line=dict(
                                            color='black',
                                            width=0.5,
                                            ),
                                        ),
                            opacity=self.params['opacity'],
                            ),
                        row=row,
                        col=col
                        )


        # MD data
        self.fig.add_trace(go.Scatter(
                            x=q,
                            y=R,
                            showlegend=showlegend,
                            name='%s; %d frames' %(contrast,nFrames),
                            mode='lines',
                            line=dict(width=3,
                                      color=self.col_light[col_orig],
                                      ),
                            opacity=self.params['opacity'],
                            ),
                        row=row,
                        col=col
                        )

    def reorder(self):

        self.fig.data = (self.fig.data[0],self.fig.data[2],self.fig.data[4],\
                        self.fig.data[1],self.fig.data[3],self.fig.data[5],\
                        )
