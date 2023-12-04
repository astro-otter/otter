'''
Functions to generate different plotly express plots with summary values for the catalog
'''

# imports
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.io import to_html
import astropy.coordinates as coord
import astropy.units as u

def skymap(fig, tdes):
    '''
    Generates a skymap of the locations of the tdes
    '''

    fig.update_geos(projection_type="mollweide", showland=False,
                    showcoastlines=False, showframe=True,
                    lonaxis=dict(showgrid=True, dtick=30),
                    lataxis=dict(showgrid=True, dtick=30),
                    row=1, col=2
                    )
    
    info = {'RA [deg]':[],
            'Dec [deg]':[],
            'TDE Name':[],
            'Z': []
            }

    for t in tdes:
        coords = t['coordinate']['equitorial']
        defaultcoord = None
        for c in coords:
            if 'default' in c and c['default']:
                defaultcoord = c
                break
        if defaultcoord is None:
            defaultcoord = coords[0] # just get the first as a catch all
            
        info['RA [deg]'].append(defaultcoord['ra'])
        info['Dec [deg]'].append(defaultcoord['dec'])
        info['TDE Name'].append(t['name']['default_name'])
        if 'distance' in t and 'redshift' in t['distance']:
            info['Z'].append(float(t['distance']['redshift'][0]['value']))
        else:
            info['Z'].append(0)

    info['RA [deg]'] = coord.Angle(info['RA [deg]'], unit=u.hourangle)
    info['Dec [deg]'] = coord.Angle(info['Dec [deg]'], unit=u.deg)
    
    info['RA [deg]'] = info['RA [deg]'].deg
    info['Dec [deg]'] = info['Dec [deg]'].deg
    
    for lon in np.arange(0, 360, 30):
        fig.add_trace(go.Scattergeo(lon=[lon, lon, lon],
                                    lat=[90, 0, -90],
                                    mode='lines+text',
                                    showlegend=False,
                                    text=["", f"{lon}", ""],
                                    textposition='bottom center',
                                    line=go.scattergeo.Line(width=1, dash='dot', color='grey')))
    for lat in [-90, -60, -30, 30, 60, 90]:
        fig.add_trace(go.Scattergeo(lon=[0, 180, 360], lat=[lat, lat, lat],
                                    mode='lines+text',
                                    text=[f"{lat}"],
                                    textposition="middle right",
                                    showlegend=False,
                                    line=go.scattergeo.Line(width=1, dash='dot', color='grey')))
        
    fig.add_trace(go.Scattergeo(lon=info['RA [deg]'],
                                lat=info['Dec [deg]'],
                                hovertext=info['TDE Name'],
                                hovertemplate="<b>Name</b>: %{hovertext}<br>" +
                                              "<b>RA</b>: %{lon}<br>" +
                                              "<b>Dec</b>: %{lat}" +
                                              "<extra></extra>",
                                showlegend=False,
                                marker=dict(color=np.log10(info['Z']),
                                            colorbar=dict(thickness=20,
                                                          title='log(z)'),
                                            colorscale="Thermal"
                                            )
                                ),
                  row=1, col=2)
        
    return fig
    
def redshifts(fig, tdes):
    '''
    Generates a plot of the redshifts
    '''

    z = []
    for t in tdes:
        if 'distance' in t and 'redshift' in t['distance']:
            z.append(t['distance']['redshift'][0]['value'])
    
    z = np.array(z).astype(float)
    
    fig.add_trace(go.Histogram(x=z,
                               marker={"color": "grey"},
                               showlegend=False),
                  row=1, col=1)
    return fig    

def plotCatalogSummary(tdes):
    '''
    makes a subplot with all of them and returns the html
    '''
    fig = make_subplots(rows=1, cols=2, specs=[[{"type":"histogram"},
                                               {"type":"scattergeo"}]
                                               ]
                        )
    fig.update_xaxes(title_text="Right Ascension [deg]", row=1, col=2)
    fig.update_xaxes(title_text="Redshift", row=1, col=1)
    fig.update_yaxes(title_text="Declination [deg]", row=1, col=2)
    fig.update_yaxes(title_text="Number of TDEs", row=1, col=1)
    
    #fig.update_layout(title = {'text':f'Number of TDEs: {len(tdes)}'})
    
    fig = skymap(fig, tdes)
    fig = redshifts(fig, tdes)

    return to_html(fig, full_html=False,
                       default_width='100%',
                       default_height='500px')


def plotPhot(tde):
    return tde.plotPhotometry()

def plotSpec(tde):
    return ''
