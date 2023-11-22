"""
Inspiration from:
https://levelup.gitconnected.com/python-classes-to-standardize-plotly-figure-formatting-123fe35c8d2d
"""

import numpy as np
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
        self.color_dict_violin = dict(MC3_ph7 = '#116C6E',
                                      MC3H_ph3 = '#8F011B',
                                      MC3H_10p = '#4A8515',
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
                              nticks=3,
                              type=axisType,
                              exponentformat='power',
                              title_standoff=40,
                              row=row,
                              col=col
                              )


    def update_yaxis(self,ylabel='',showticklabels=True,axisType='linear',row=1,col=1):

        self.fig.update_yaxes(title_text='<b>'+ylabel+'</b>',
                          # range=[-20,200],
                          showline= True,
                          linecolor= 'black',
                          linewidth=5,
                          showticklabels=showticklabels,
                          ticks= 'inside',
                          mirror='allticks',
                          tickwidth=5,
                          tickcolor='black',
                          tickfont=self.tick_dict,
                          nticks=3,
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
            width = 700, height = 700,
            autosize = False
            )

        self.fig.update_traces(showscale=False)

        self.fig.update_xaxes(ticks='',)
        self.fig.update_yaxes(ticks='',)




    def plotViolin_halves(self,df,sysComp,cholPerc,row=1,col=1,showlegend=False):

        # identify which pH the system is in to access correct mc3 colour
        if 'MC3' in sysComp and '10P' not in sysComp: suffix = '_ph7'
        if 'MC3H' in sysComp and '10P' not in sysComp: suffix = '_ph3'
        if 'MC3H' in sysComp and '10P' in sysComp: suffix = '_10p'

        # update cationic ionisable lipid identifier
        cil_str = sysComp[0]

        # if more than one column, plot mc3 and then cholesterol data
        if len(df.columns) > 1:

            self.fig.add_trace(go.Violin(
                                    y=df.loc[:,cil_str],
                                    side='negative',
                                    line_color=self.color_dict_violin.get(cil_str+suffix),
                                    showlegend=showlegend,
                                    legendgroup=cil_str, scalegroup=cil_str, name=cholPerc,
                                    ),
                                    row=row,
                                    col=col,
                                    )

            self.fig.add_trace(go.Violin(
                                    y=df.loc[:,'CHOL'],
                                    side='positive',
                                    line_color=self.color_dict_violin.get('CHOL'),
                                    showlegend=showlegend,
                                    legendgroup='CHOL', scalegroup='CHOL', name=cholPerc,
                                    ),
                                    row=row,
                                    col=col,
                                    )

        # else no cholesterol, plot mc3 data only
        else:
            self.fig.add_trace(go.Violin(
                                    y=df.loc[:,cil_str],
                                    side='negative',
                                    line_color=self.color_dict_violin.get(cil_str+suffix),
                                    showlegend=showlegend,
                                    legendgroup=cil_str, scalegroup=cil_str, name=cholPerc,
                                    ),
                                    row=row,
                                    col=col,
                                    )


        self.fig.update_xaxes(type='category')
        self.fig.update_traces(box_visible=True, meanline_visible=True)
        self.fig.update_layout(violinmode='overlay')
        self.fig.update_layout(yaxis_nticks=5)




    def plotViolin_full(self,df,sysComp,cholPerc,row=1,col=1,showlegend=False):

        # identify which pH the system is in to access correct mc3 colour
        suffix = '_ph3' if 'MC3H' in sysComp else '_ph7'

        self.fig.add_trace(go.Violin(
                                y=df.loc[:,'delta_av'],
                                line_color=self.color_dict_violin.get('DLMC3'+suffix),
                                showlegend=showlegend,
                                legendgroup='DLMC3', scalegroup='DLMC3', name=cholPerc,
                                ),
                                row=row,
                                col=col,
                                )

        self.fig.update_xaxes(type='category')
        self.fig.update_traces(box_visible=True, meanline_visible=True)




    def plotLogLog(self,X,Y,sysName,sysComp,speciesList,row=1,col=1,showlegend=True):

        self.fig.update_layout(
                margin=dict(l=20,r=20,t=20,b=20),
                autosize=False,
                width=700,
                height=600,
                plot_bgcolor='white',
                showlegend = True,
                legend=dict(
                            x=0.4,
                            y=1.0,
                            bgcolor='rgba(0,0,0,0)'
                            ),
                legend_font_size=20,
                )

        for idx, species in enumerate(speciesList):

            self.fig.add_trace(go.Scatter(
                                x=X[idx],
                                y=Y[idx],
                                showlegend=True,
                                name=str(sysName) + ' (' + str(species) + ')',
                                mode='lines',
                                line=dict(width=5,
                                          color=self.color_dict_loglog.get(species),
                                          dash='solid' if species != 'TIP3_lower' else 'dot'),
                                marker=dict(size=15,
                                            symbol='square',
                                            ),
                                opacity=self.params['opacity'],
                                ),
                            row=row,
                            col=col
                            )
        # self.fig.add_vrect(x0=0.5, x1=10, line_width=0, fillcolor="red", opacity=0.2)
        self.fig.update_layout(yaxis_nticks=5)
        self.fig.update_layout(xaxis_nticks=5)




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
