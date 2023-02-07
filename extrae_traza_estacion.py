


from obspy import UTCDateTime
from obspy.clients import fdsn, seedlink
from obspy.clients.fdsn.client import Client
#from obspy.clients.seedlink.basic_client import Client #llama Client no por fdsn sino por seedlink.basic_client
from obspy.core.stream import Stream
from obspy import read
import os
from io import open
from datetime import datetime, timedelta
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, DateFormatter, date2num, num2date

def extrae_traza(ancho_linea,time_i,estacion,ventana,net,loc,ch,nacional,obt_ppsd,log_linear,escalar,guardar_traza):

		
	t1=UTCDateTime(time_i)
	t2 = t1+(ventana*60) #la ventana se da en minutos

	if nacional==True:
		client_fdsn = fdsn.Client("http://10.100.100.232:8091")
		st = client_fdsn.get_waveforms("CM", estacion, loc, ch, t1, t2)
		
		inv = client_fdsn.get_stations(network="CM", station=estacion,starttime=t1,endtime=t2, channel=ch, level='response')
		
	if nacional==False:
		#PARA TRAER estaciones desde iris
		client = Client("IRIS")
		
		st = client.get_waveforms(net,estacion,loc, ch, t1, t2)
		inv = client.get_stations(network=net, station=estacion,starttime=t1, endtime=t2,channel=ch, level='response')
	


	fig1 = plt.figure(figsize=(21., 11.))
	
	print(st)

	gaps_overlaps=st.get_gaps()

	#--- codigo angel
	gaps_a = [gap_a[6] for gap_a in st.get_gaps()]
	gaps_b = []
	gap_sum_a = np.sum(gaps_a)
	


	print(gap_sum_a)
	#
	#
	
	time_i_gaps_HHZ=[]
	time_i_gaps_HHN=[]
	time_i_gaps_HHE=[]
	time_i_over_HHZ=[]
	time_i_over_HHN=[]
	time_i_over_HHE=[]

	time_f_gaps_HHZ=[]
	time_f_gaps_HHN=[]
	time_f_gaps_HHE=[]
	time_f_over_HHZ=[]
	time_f_over_HHN=[]
	time_f_over_HHE=[]

	num_gaps_z=0
	num_gaps_n=0
	num_gaps_e=0
	num_over_z=0
	num_over_n=0
	num_over_e=0
	num_overlaps=0

	for go in gaps_overlaps:

		print (go)
		
		ch_go=go[3]
					
		if go[-2] > 0: #mayor a cero gap
			
			

			time_i_gap=go[4].datetime
			time_f_gap=go[5].datetime
			
			if ch_go=="HHZ":
				time_i_gaps_HHZ.append(time_i_gap)
				time_f_gaps_HHZ.append(time_f_gap)
				num_gaps_z=num_gaps_z+1
				num_over_z=num_over_z+1
			if ch_go=="HHN":
				time_i_gaps_HHN.append(time_i_gap)
				time_f_gaps_HHN.append(time_f_gap)
				num_gaps_n=num_gaps_n+1
				num_over_n=num_over_n+1	
			if ch_go=="HHE":
				time_i_gaps_HHE.append(time_i_gap)
				time_f_gaps_HHE.append(time_f_gap)
				num_gaps_e=num_gaps_e+1
				num_over_e=num_over_e+1
			
		if go[-2] < 0: #menor a cero overlap
			
			num_overlaps=num_overlaps+1

			time_i_over=go[4].datetime
			time_f_over=go[5].datetime
			
			if ch_go=="HHZ":
				time_i_over_HHZ.append(time_i_over)
				time_f_over_HHZ.append(time_f_over)
			if ch_go=="HHN":
				time_i_over_HHN.append(time_i_over)
				time_f_over_HHN.append(time_f_over)	
			if ch_go=="HHE":
				time_i_over_HHE.append(time_i_over)
				time_f_over_HHE.append(time_f_over)


	print("Numero de Gaps Z: "+str(num_gaps_z))
	print("Numero de Gaps N: "+str(num_gaps_n))
	print("Numero de Gaps E: "+str(num_gaps_e))
	print("Numero de Overlaps Z: "+str(num_over_z))
	print("Numero de Overlaps  N: "+str(num_over_n))
	print("Numero de Overlaps  E: "+str(num_over_e))
	st.merge(method=0, fill_value='interpolate', interpolation_samples=0)

	
	c=1

	for tr in st:		
		
		#tr.filter("bandpass", freqmin=0.0015, freqmax=0.035, corners=4, zerophase=False) #aplica filtro but pasabanda
		data=tr.data
		vec_t=[]
		ti=tr.stats.starttime
		delta=tr.stats.delta
		#print("delta",1/delta)
		for i in range(0,len(data)):
			#print(ti+(i*delta))
			vec_t.append((ti+(i*delta)).datetime)

		ax1 = fig1.add_subplot(3,1,c) 		
		ax1.plot(vec_t,data,linewidth=ancho_linea,color="k",alpha=0.9)
		#ax1.axhline(y=0,linewidth=0.3,color="k",linestyle="--")

		ax1.set_yscale(log_linear)

		if ventana <= 80:
			interval_min=5
			interval_hor=1			
		if 80 < ventana <= 180:
			interval_min=15
			interval_hor=1			
		if 180 < ventana <= 360:
			interval_min=20
			interval_hor=2
		if ventana > 360:
			interval_min=40
			interval_hor=3
			

		ax1.xaxis.set_minor_locator( MinuteLocator(interval = interval_min))
		ax1.xaxis.set_minor_formatter( DateFormatter('%H:%M') )
		ax1.xaxis.set_major_locator( HourLocator(interval = interval_hor) )
		ax1.xaxis.set_major_formatter( DateFormatter('\n%b%d-%H:%M UT'))
		
		d=max(tr.data)-min(tr.data)
		max_ylim=max(tr.data)+(0.1*d)
		min_ylim=min(tr.data)-(0.1*d)
		offs=round(np.mean(tr.data),1)
		ax1.axhline(y=offs,linewidth=0.3,color="k",linestyle="--")
		if escalar==True:
			ax1.set_ylim(min_ylim,max_ylim)
		ax1.text(t1.datetime+timedelta(minutes=2),max(tr.data),"\nOffset "+str(offs)+", Cuentas: Max "+str(max(tr.data)) + ", Min "+str(min(tr.data)),fontsize=8)
		#ax1.text(t1.datetime+timedelta(minutes=2),offs,"Offset "+str(offs),fontsize=8,color="b")
		ax1.set_ylabel(estacion+"_"+tr.stats.channel)

		if tr.stats.channel=="HHZ":
			ax1.vlines(time_i_gaps_HHZ, min_ylim, max_ylim, colors='r',linewidth=0.5,alpha=0.6)
			ax1.vlines(time_f_gaps_HHZ, min_ylim, max_ylim, colors='r',linewidth=0.5,alpha=0.6)
			ax1.vlines(time_i_over_HHZ, min_ylim, max_ylim, colors='g',linewidth=0.5,alpha=0.6)
			ax1.vlines(time_f_over_HHZ, min_ylim, max_ylim, colors='g',linewidth=0.5,alpha=0.6)
		if tr.stats.channel=="HHN":
			ax1.vlines(time_i_gaps_HHN, min_ylim, max_ylim, colors='r',linewidth=0.5,alpha=0.6)
			ax1.vlines(time_f_gaps_HHN, min_ylim, max_ylim, colors='r',linewidth=0.5,alpha=0.6)
			ax1.vlines(time_i_over_HHN, min_ylim, max_ylim, colors='g',linewidth=0.5,alpha=0.6)
			ax1.vlines(time_f_over_HHN, min_ylim, max_ylim, colors='g',linewidth=0.5,alpha=0.6)
		if tr.stats.channel=="HHE":
			ax1.vlines(time_i_gaps_HHE, min_ylim, max_ylim, colors='r',linewidth=0.5,alpha=0.6)
			ax1.vlines(time_f_gaps_HHE, min_ylim, max_ylim, colors='r',linewidth=0.5,alpha=0.6)
			ax1.vlines(time_i_over_HHE, min_ylim, max_ylim, colors='g',linewidth=0.5,alpha=0.6)
			ax1.vlines(time_f_over_HHE, min_ylim, max_ylim, colors='g',linewidth=0.5,alpha=0.6)
			

		print("Offset "+tr.stats.channel+" ",offs) 	


		if obt_ppsd==True:
			ppsd = PPSD(tr.stats, metadata=inv)
			ppsd.add(st)			
			try:		
				
				ppsd.plot(cmap=pqlx, xaxis_frequency=True, show_mean=True,show=False)
			except:
					
				print("El tiempo de extraccion es muy corto para calcular PPSD")			

		c=c+1
	
		
	plt.subplots_adjust(hspace=0.35)
	plt.show()
		

	if guardar_traza==True:
		extrae_traza=os.getcwd()+"/traza_extraida" #crea carpeta donde extrae continuo cortado cada media hora
		dir_est=extrae_traza+"/"+estacion

		if not os.path.exists(extrae_traza):
			os.makedirs(extrae_traza)
			print(f"\nSe creo el directorio {extrae_traza}, donde se almacenarán las formas de onda extraidas en la estacion {estacion}")

		if not os.path.exists(dir_est):
			os.makedirs(dir_est)
			print(f"Se creo el directorio {dir_est}, donde se almacenarán las formas de onda extraidas\n")

				#---------------------------

		if len(st)>0:
			st_Stream = Stream(traces=st) #con este comando se unen trazas en un Stream
			name_time=t1.datetime.strftime("%Y%m%d%H%M%S")
			#st_Stream.write(dir_est+"/"+estacion+"_"+name_time+"_"+str(ventana)+"min.mseed", format="MSEED") #activar para guardar traza
			#print("Se guardo la traza "+name_time+"_"+str(ventana)+"min.mseed"+" en el directorio "+dir_est)
		

time_i="2022-08-09 03:00:00" #UT tiempo de inicio de extraccion
ventana=60*6 # (numero de horas despues del 60*)

net="CM" #sirve para internacionales, cuando es nacional por defecto se pone CM
estacion="BLMC" #BLMC 
loc="00"
ch="HH*" #Internacionales suelen ser BH

nacional=True
obt_ppsd=True #si se quiere calcular el ppsd

log_linear="linear" #elegir la escala Y: "log" or "linear"
escalar=True #escalar los 3 ejes a los mismos limites en Y
guardar_traza=False
ancho_linea=0.4

extrae_traza(ancho_linea,time_i,estacion,ventana,net,loc,ch,nacional,obt_ppsd,log_linear,escalar,guardar_traza) #True si es nacional