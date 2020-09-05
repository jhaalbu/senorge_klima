import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.cbook as cbook
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.dates import DateFormatter
#from windrose import WindroseAxes
import requests
import datetime
import numpy as np
import folium
from pyproj import CRS
from pyproj import Transformer

import streamlit as st
from streamlit_folium import folium_static
import folium


def nve_api(lat, lon, startdato, sluttdato, para):
    """Henter data frå NVE api GridTimeSeries

    Args:
        lat (str): øst-vest koordinat (i UTM33)
        output er verdien i ei liste, men verdi per dag, typ ne
        lon (str): nord-sør koordinat (i UTM33)
        startdato (str): startdato for dataserien som hentes ned
        sluttdato (str): sluttdato for dataserien som hentes ned
        para (str): kva parameter som skal hentes ned f.eks rr for nedbør
        
    Returns:
        verdier (liste) : returnerer i liste med klimaverdier
        
    """

    api = 'http://h-web02.nve.no:8080/api/'
    url = api + '/GridTimeSeries/' + str(lat) + '/' + str(lon) + '/' + str(startdato) + '/' + str(sluttdato) + '/' + para + '.json'
    r = requests.get(url)

    verdier = r.json()
    return verdier

def transformer(lat, lon):
    transformer = Transformer.from_crs(5973, 4326)
    trans_x, trans_y =  transformer.transform(lat, lon)
    return trans_x, trans_y


def klima_dataframe(lat, lon, startdato, sluttdato):
    """Funksjonen tar inn ei liste  med lister med klimaparameter og lager dataframe

    Args:
        lat (str): øst-vest koordinat (i UTM33)
        output er verdien i ei liste, men verdi per dag, typ ne
        lon (str): nord-sør koordinat (i UTM33)
        startdato (str): startdato for dataserien som hentes ned
        sluttdato (str): sluttdato for dataserien som hentes ned
        para (str): kva parameter som skal hentes ned f.eks rr for nedbør

    Returns
        df (dataframe): Returnerer ei pandas dataframe
    """
    bar = st.progress(0)
    rr = nve_api(lat, lon, startdato, sluttdato, 'rr') #Nedbør døgn
    bar.progress(20)
    fsw = nve_api(lat, lon, startdato, sluttdato, 'fsw') #Nysnø døgn
    bar.progress(40)
    sdfsw3d = nve_api(lat, lon, startdato, sluttdato, 'sdfsw3d') #Nynsø 3 døgn
    bar.progress(60)
    sd = nve_api(lat, lon, startdato, sluttdato, 'sd') #Snødybde
    bar.progress(80)
    tm = nve_api(lat, lon, startdato, sluttdato, 'tm') #Døgntemperatur
    bar.progress(100)
    start = datetime.datetime(int(startdato[0:4]), int(startdato[5:7]), int(startdato[8:10]))#'1960-01-01' 
    end = datetime.datetime(int(sluttdato[0:4]), int(sluttdato[5:7]), int(sluttdato[8:10]))
    #Etablerer pandas dataframe frå rr liste
    df = pd.DataFrame(rr['Data'])
    #Lager kolonne med datoer
    df['dato'] = pd.date_range(start, end)
    #Gir nytt navn til rr kolonna
    df.rename({0 : rr['Theme']}, axis=1, inplace=True)
    #Setter datokolonna til indexkolonna, slik at det går an å bruke grouperfunksjon for å sortere på tidsserie
    df.set_index('dato', inplace=True)
    #Etablerer kolonner i dataframefor andre parameter
    df[fsw['Theme']] = fsw['Data']
    df[sd['Theme']] = sd['Data']
    df[sdfsw3d['Theme']] = sdfsw3d['Data']
    df[tm['Theme']] = tm['Data']
    df['rr3'] = df.rr.rolling(3).sum() #Summerer siste 3 døgn 
    df[df > 60000] = 0
    return df



#OBS!! Datoer for vinddata oppgis lenger ned i arket (eksister ikkje i like lang periode)
startdato = '1958-01-01' 
sluttdato = '2019-12-31'

def plot_normaler(df, ax1=None):
    #Lager dataframe med månedsvise mengder av temperatur og nedbør GRAF1
    #df['month'] = pd.DatetimeIndex(df.index).month #Lager ei kolonne med månedsnummer

    mon_rr = df['rr'].groupby(pd.Grouper(freq='M')).sum() #Grupperer nedbør etter måneder per år og summerer
    mon_tm = df['tm'].groupby(pd.Grouper(freq='M')).mean() #Grupperer temperatur etter måneder per år og tar snitt
    month_rr = mon_rr.to_frame() #Lager dataframe, unødvendig? Eklere å plotte?
    month_tm = mon_tm.to_frame() #Lager dataframe, unødvendig? Eklere å plotte?
    month_rr['m'] = pd.DatetimeIndex(month_rr.index).month #lager kolonne for månedsnummer
    month_tm['m'] = pd.DatetimeIndex(month_tm.index).month #Lager kolonne for månedsnummer
    month_mean_tm = month_tm.groupby(['m']).mean()
    if ax1 is None:
        ax1 = plt.gca()
    ax1.set_title('Gjennomsnittlig månedsnedbør og temperatur ' + startdato[0:4] + ' til ' + sluttdato[0:4])
    ax1.bar(month_rr['m'], month_rr['rr'], width=0.5, snap=False)
    ax1.set_xlabel('Måned')
    ax1.set_ylabel('Nedbør (mm)')
    ax1.set_ylim(0, month_rr['rr'].max()+50)
    #ax1.text('1960', aar_df['rr'].max()+20, "Gjennomsnittlig månedsnedbør:  " + str(int(snitt)) + ' mm')

    ax2 = ax1.twinx()#Setter ny akse på høgre side 
    ax2.plot(month_mean_tm.index, month_mean_tm['tm'], 'r', label='Gjennomsnittstemperatur', linewidth=3.5)
    ax2.set_ylim(month_mean_tm['tm'].min()-2, month_mean_tm['tm'].max()+5)
    ax2.set_ylabel(u'Temperatur (\u00B0C)')
    ax2.yaxis.set_tick_params(length=0)
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax2.get_yaxis().set_visible(True)
    ax2.legend()

    return ax1, ax2

st.title('AV-Klima')
st.write("Test av klimadata streamlit")

#Gi in kordinater for posisjon og start og sluttdato for dataserien.
lon = st.text_input("Gi NORD koordinat")
#lon = 6822565  #Y
lat = st.text_input("Gi ØST koordinat")
#lat = 67070      #X

if lat and lon:
    transformer = Transformer.from_crs(5973, 4326)
    trans_x, trans_y =  transformer.transform(lat, lon)

    m = folium.Map(location=[trans_x, trans_y], zoom_start=10) #Bruker transformerte koordinater (ikkej funne ut korleis bruke UTM med folium)
    folium.Circle(
        radius=1000,
        location=[trans_x, trans_y],
        fill=False,
    ).add_to(m)
    #Legger til Norgeskart WMS
    folium.raster_layers.WmsTileLayer(
        url='https://opencache.statkart.no/gatekeeper/gk/gk.open_gmaps?layers=topo4&zoom={z}&x={x}&y={y}',
        name='Norgeskart',
        fmt='image/png',
        layers='topo4',
        attr=u'<a href="http://www.kartverket.no/">Kartverket</a>',
        transparent=True,
        overlay=True,
        control=True,
        
    ).add_to(m)
    folium_static(m)

    df = klima_dataframe(lat, lon, startdato, sluttdato) 
    st.write(df)


    fig, ax = plt.subplots(1)
    ax1, ax2 = plot_normaler(df)

    st.write(fig)
#plot_something(data1, ax1, color='blue')
#plot_something(data2, ax2, color='red')
