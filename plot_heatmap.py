import sys
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import contacts_data as contacts_data

pd.set_option('display.max_columns', None)
pio.kaleido.scope.default_format = "svg"

class Figure:
    def __init__(self):
        # initialise figure
        self.fig = go.Figure()

    def plot_heat(self, X, Y, Z, Z_str, show_legend=False):
        self.fig.add_trace(go.Heatmap(
            x=X,
            y=Y,
            z=Z,
            colorbar=dict(
                title=cbar_title,
                ticks='',
                tickwidth=2.5,
                len=1.045,
                yanchor='top',
                y=1.025,
                x=1,
            ),
            # zmin=ZMIN,
            # zmax=ZMAX,
            text=Z_str,
            texttemplate='%{text}',
            textfont=dict(size=20),
        ))
        self.fig.update_layout(xaxis_nticks=len(X))
        self.fig.update_layout(yaxis_nticks=len(Y))
        self.fig.update_xaxes(type='category')
        self.fig.update_yaxes(type='category')

    def update_fig_design(self):
        self.fig.update_layout(font=dict(family='Latex', color='black'),
                               autosize=False,
                               width=700,
                               height=600,
                               margin=dict(l=20, r=20, t=20, b=20),
                               plot_bgcolor='white',
                               showlegend=True,
                               legend=dict(
                                   x=0.5,
                                   y=legend_y,
                                   bgcolor='rgba(0,0,0,0)',
                                   orientation="h",
                                   xanchor="center",
                               ),
                               legend_font_size=fs_legend,
                               barmode='group',
                               xaxis_tickangle=x_tick_angle,
                               # title=fig_title,
                               )
        self.fig.update_xaxes(
            title_text=label_xaxis,
            showline=True,
            linecolor='black',
            linewidth=2.5,
            ticks='inside',
            mirror='allticks',
            tickwidth=2.5,
            ticklen=8,
            tickcolor='black',
            tickfont=dict(family='Latex', size=fs_xaxis, color='black'),
            side=xaxis_side,
            title_standoff=standoff_x,
        )
        self.fig.update_yaxes(
            title_text=label_yaxis,
            automargin=True,
            showline=True,
            linecolor='black',
            linewidth=2.5,
            ticks='inside',
            minor_ticks="inside",
            mirror='allticks',
            tickwidth=2.5,
            ticklen=8,
            tickcolor='black',
            tickfont=dict(family='Latex', size=fs_yaxis, color='black'),
            title_standoff=standoff_y,
        )

    def save_figure(self, title):
        output_dir = '../output/' + title + '.png'
        self.fig.write_image(output_dir, scale=5)

    def show_figure(self):
        self.fig.show(renderer='png')


# ----------------------------------------------------------------------------#
# PREAMBLE

# set axis range for y axis
yRange = True
yrange_lb = 0  # total = -2, depletion = -1.1, insertion = 0.0
yrange_ub = 90  # total = 4, depletion = 0.0, insertion = 3.0

# colorbar properties
ZMIN = 40  # apm=40, thick=10, scc=0, enrich=0,
ZMAX = 150  # apm=150, thick=25, scc=0.1, enrich=1.5

cbar_title = 'Counts'

err_width = 5.0
line_width = 2.0

# marker opacity
opacity = 1.0

# define the title of the legend
legend_title = ''

# define legend position
legend_x = 0.01  # 0.74 / 0.64
legend_y = 1.0  # 0.01 / 0.6 / 1.0

# define legend font size
fs_legend = 24  # 32
fs_xaxis = 24
fs_yaxis = 28

# axis label distance from axis
standoff_x = 25
standoff_y = 20

x_tick_angle = 0  # 0 / 45
xaxis_side = 'bottom'

# ----------------------------------------------------------------------------#

MC3H = "NP CN GLA CX C1A C1B D2A D2B D3A D3B C4A C4B"
LI5H = "OH NP CA CB GLA GLB CBX CA1 CB1 CB3 CA2 CB2 CB4"
CHOL = "ROH R1 R2 R3 R4 R5 R6 C1 C2"
RNA = "BB BB1 BB2 SC1 SC2 SC3 SC4 SC5 SC6"

name = contacts_data.name
title = f'heatmap_{name}'

label_xaxis = r'$\huge{\text{RNA}}$'
label_yaxis = r'$\huge{\text{CHOL}}$'

X = RNA.split(' ')
Y = CHOL.split(' ')

df = contacts_data.contacts
# df_str = contacts_data.std
df_str = None

fig = Figure()

fig.__init__()

fig.plot_heat(X=X, Y=Y, Z=df, Z_str=df_str, show_legend=False)

fig.update_fig_design()

fig.save_figure(title)

# ----------------------------------------------------------------------------#
