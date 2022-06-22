import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
import seaborn as sns

def plot_3d_meristem(df, gene, plot_title = "plot title", legend_title = "legend title", 
                     show_grid = True, show_legend = True, colorscale = ["yellow","red"], 
                     show_trace = False, backgrnd_col="white", font_col="black"):
    if show_grid:
        xaxis_grid=dict(showgrid=True, zeroline=False, visible=True,
                        gridcolor=font_col, showbackground=False,
                        tickfont=dict(color=font_col), title=dict(font=dict(color=font_col)))
        yaxis_grid=dict(showgrid=True, zeroline=False, visible=True,
                        gridcolor=font_col, showbackground=False,
                        tickfont=dict(color=font_col), title=dict(font=dict(color=font_col)))
        zaxis_grid=dict(showgrid=True, zeroline=False, visible=True, 
                        gridcolor=font_col, showbackground=False,
                        tickfont=dict(color=font_col), title=dict(font=dict(color=font_col)))
    else:
        xaxis_grid=dict(showgrid=False, zeroline=False, visible=False)
        yaxis_grid=dict(showgrid=False, zeroline=False, visible=False)
        zaxis_grid=dict(showgrid=False, zeroline=False, visible=False)

        

    fig = go.Figure(data=[go.Scatter3d(x=df['x'], y=df['y'], z=df['z'], 
                                       mode='markers',
                                       marker=dict(color=df[gene],
                                                   showscale=show_legend,
                                                   colorscale=colorscale,
                                                   colorbar=dict(thickness=30, tickcolor=font_col,
                                                                  tickfont=dict(color=font_col))),
                                       showlegend=show_trace)],
                    layout=go.Layout(title=dict(text=plot_title, font=dict(size=30, color=font_col)), 
                                     scene_dragmode='orbit',
                                     scene=dict(bgcolor=backgrnd_col,
                                               xaxis=xaxis_grid, yaxis=yaxis_grid, zaxis=zaxis_grid),
                                     
                                     margin=dict(l=0, r=0, t=50, b=0),
                                     paper_bgcolor=backgrnd_col,
                                     scene_camera=dict(up=dict(x=0, y=0, z=1), 
                                                      eye=dict(x=2, y=1, z=1), # zooming
                                                      center=dict(x=0, y=0, z=0))
                                    )

                    )
    # Custom Legend
    MAX = float(str(max(df[gene]))[0:4])
    MIN = float(str(min(df[gene]))[0:4])
    MIDDLE = float(str(MIN+((MAX - MIN)/2))[0:4])
    fig['data'][0]['marker'].update(colorbar=dict(thickness=50, len=0.5, tickmode='array', 
                                                  tickvals=[MIN,MIDDLE,MAX],
                                                 tickfont=dict(size=20),
                                                 title=dict(text=legend_title, font=dict(size=20, color=font_col))))
    return fig


def plot_umap(umap, plot_size_x=7, plot_size_y=7, label_size=10):
    groups = list(set(umap.loc[:,"cluster_no"].to_list()))
    customPalette = sns.color_palette("hls", len(groups))
    data = umap
    plt.figure(figsize=(plot_size_x,plot_size_y))
    #loop through labels and plot each cluster
    for cluster_no in groups:
    # for i, label in enumerate(groups):
        # print(label)
        # add data points 
        plt.scatter(x=data.loc[data['cluster_no']==cluster_no, 'UMAP_1'], 
                    y=data.loc[data['cluster_no']==cluster_no,'UMAP_2'], 
                    color=customPalette[cluster_no], 
                    alpha=1)

        #add label
        plt.annotate(cluster_no, 
                     data.loc[data['cluster_no']==cluster_no,['UMAP_1','UMAP_2']].mean(),
                     horizontalalignment='center',
                     verticalalignment='center',
                     size=label_size, weight='bold',
                     color='white',
                     backgroundcolor=customPalette[cluster_no],
                    ) 