
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.animation import FFMpegWriter
# import geopandas as gpd
# from mpl_toolkits.basemap import Basemap, shiftgrid
import sys
import random

import datetime


fig = plt.figure()
# fig, axs = plt.subplots(2, 1) # For Map


# world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

# m = Basemap(width=12000000,height=9000000,projection='lcc',resolution='c',lat_1=45.,lat_2=55,lat_0=50,lon_0=-107.)
# m.drawcoastlines()
# m.drawmapboundary(fill_color='aqua')
# m.fillcontinents(color='coral',lake_color='aqua')

# m = Basemap(width=12000000,height=9000000,projection='lcc',resolution='l',lat_1=45.,lat_2=55,lat_0=50,lon_0=-107.)
# m = Basemap(resolution='l', ax=axs[1])
# m.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
# m.bluemarble()

# m.shadedrelief()
# m.etopo()
# m.drawcoastlines()



    
def greatCircleDistance_km(pnt1, pnt2):
    """
    This uses the special case of the vincenty formula
    https://en.wikipedia.org/wiki/Great-circle_distance
    
    This can be the heuristic for A*
    """
    # print("GREAT",pnt1)
    # pnt1, pnt2 are Latitude-Longitude Pairs in units of decimal degrees
    # phi1 = 2*np.pi * pnt1[0] / 360
    # phi2 = 2*np.pi * pnt2[0] / 360
    toRadians = lambda x:x*np.pi/180
    phi1, phi2, lam1, lam2 = map(toRadians, [pnt1[0], pnt2[0], pnt1[1], pnt2[1]])
    # dLam = 2*np.pi * np.abs(pnt2[1]-pnt1[1]) / 360
    dLam = np.abs(lam2-lam1)
    earthRadius_km = 6371.009
    
    x = np.sqrt((np.cos(phi2)*np.sin(dLam))**2 + 
        (np.cos(phi1)*np.sin(phi2) - np.sin(phi1)*np.cos(phi2)*np.cos(dLam))**2)
    y = np.sin(phi1)*np.sin(phi2) + np.cos(phi1)*np.cos(phi2)*np.cos(dLam)
    deltaSigma_rad = np.arctan2(x,y)
    length_km = earthRadius_km * deltaSigma_rad
    
    return length_km
    
def greatCircleVector(pnt, dir, dist):
    # pnt is a Latitude-Longitude Pair in decimal degrees, dir is direction in degrees, dist is length in km. Returns point at distance from pnt along dir
    pass
    
def greatCircleWaypoint(pnt1, pnt2, dist):
    # pnt1, pnt2 are Latitude-Longitude Pairs in units of decimal degrees, dist is length in km. Returns point along the line
    pass
    
def ecefVector(pnt):
    # pnt is a Latitude-Longitude Pair
    lat = np.pi * (90-pnt[0]) / 180
    lon = np.pi * pnt[1] / 180

    nx = np.cos(lat)*np.cos(lon)
    ny = np.cos(lat)*np.sin(lon)
    nz = np.sin(lat)
    return np.array([nx,ny,nz])

def latLonPair(x,y,z):
    # https://gis.stackexchange.com/a/304028
    r = np.sqrt(x**2 + y**2 + z**2)
    clat = np.arccos(z/r)/np.pi*180
    lat = 90.0 - clat
    lon = np.arctan2(y,x)/np.pi*180
    lon = (lon+360)%360
    return np.array([lat,lon]) #,np.ones(lat.shape)])

def intermediatePoint(pnt1, pnt2, fraction):
    if(fraction == 0): return np.array([pnt1[0],pnt1[1]])
    
    toRadians = lambda x:x*np.pi/180
    phi1, phi2, lam1, lam2 = map(toRadians, [pnt1[0], pnt2[0], pnt1[1], pnt2[1]])

    dLam = lam2-lam1
        
    x = np.sqrt((np.cos(phi2)*np.sin(dLam))**2 + (np.cos(phi1)*np.sin(phi2) - np.sin(phi1)*np.cos(phi2)*np.cos(dLam))**2)
    y = np.sin(phi1)*np.sin(phi2) + np.cos(phi1)*np.cos(phi2)*np.cos(dLam)
    deltaSigma_rad = np.arctan2(x,y)
    
    A =np.sin((1-fraction)*deltaSigma_rad)/np.sin(deltaSigma_rad)
    B =np.sin((fraction)*deltaSigma_rad)/np.sin(deltaSigma_rad)
    
    x = A * np.cos(phi1) * np.cos(lam1) + B * np.cos(phi2) * np.cos(lam2)
    y = A * np.cos(phi1) * np.sin(lam1) + B * np.cos(phi2) * np.sin(lam2)
    z = A * np.sin(phi1)                + B * np.sin(phi2)
    
    phi3 = np.arctan2(z, np.sqrt(x**2 + y**2))
    lam3 = np.arctan2(y,x)
    lat = phi3*180/np.pi
    lon = lam3*180/np.pi
    return np.array([lat,lon]) # ,np.ones(lat.shape)])


class BaseAircraftClass():
    def __init__(self):
        self.aircraft_name         = None
        self.aircraft_mass_kg      = None
        self.batt_capacity_kJ      = None
        self.batt_threshold_full   = None
        self.batt_threshold_cutoff = None
        self.solar_MPP_W           = None
        self.power_consup_W        = None
        self.cruise_airspeed_m_s   = None
        self.climb_rate_m_s        = None
        self.sink_rate_m_s         = None
        self.cruise_altitude_m     = None
        self.sprint_airspeed_m_s   = None
    
class GannetSUWAVE(BaseAircraftClass):
    # Gannet SUWAVE TVBS:
    def __init__(self):
        # super().__init__()
        self.aircraft_name         = "Gannet"
        self.aircraft_mass_kg      = 2.5
        self.batt_capacity_kJ      = 550 # 4s10Ah = 550
        self.batt_threshold_full   = 0.65
        self.batt_threshold_cutoff = 0.15
        self.solar_MPP_W           = 71.4 # Ideal-Conditions solar array power production, including debris and attenuation
        self.power_consup_W        = 85.7 # Typical Aircraft power consumption
        self.cruise_airspeed_m_s   = 14.5 # Accurate for Gannet Prime
        self.climb_rate_m_s        = 0.6  # Total Guess
        self.sink_rate_m_s         = -1.1 # Total Guess
        self.cruise_altitude_m     = 2000 # 
        self.sprint_airspeed_m_s   = lambda x : np.sqrt(x*0.4)+self.cruise_airspeed_m_s

class LaysianQuadplane(BaseAircraftClass):
    # Laysian Solar Quadplane VTOL:
    def __init__(self):
        # super().__init__()
        self.aircraft_name         = "Laysian24 Quadplane"
        self.aircraft_mass_kg      = 4.0
        self.batt_capacity_kJ      = 550 # 4s10Ah = 550   225 for half
        self.batt_threshold_full   = 0.95 # 0.65 at 550
        self.batt_threshold_cutoff = 0.15  # 
        self.solar_MPP_W           = 71.4 # 95.2 for -32, 71.4 for -24
        self.power_consup_W        = 55.7
        self.cruise_airspeed_m_s   = 12.5
        self.climb_rate_m_s        = 0.6
        self.sink_rate_m_s         = -1.1
        self.cruise_altitude_m     = 2000
        self.sprint_airspeed_m_s   = lambda x : np.sqrt(x*0.4)+self.cruise_airspeed_m_s

class LaysianPlus28(BaseAircraftClass):
    # Laysian-Plus VTOL:
    def __init__(self):
        # super().__init__()
        self.aircraft_name         = "LaysianVTOL-28 TTRSLT"
        self.aircraft_mass_kg      = 3.7
        self.batt_capacity_kJ      = 225 # 4s10Ah = 550   225 for half
        self.batt_threshold_full   = 0.95 # 0.65 at 550
        self.batt_threshold_cutoff = 0.15  # 
        self.solar_MPP_W           = 83.3 # 3.5W per cell * 0.85 = 2.975
        self.power_consup_W        = 52.2 # 31.3W / 0.6
        self.cruise_airspeed_m_s   = 13.3
        self.climb_rate_m_s        = 0.6
        self.sink_rate_m_s         = -1.1
        self.cruise_altitude_m     = 2000
        self.sprint_airspeed_m_s   = lambda x : np.sqrt(x*0.4)+self.cruise_airspeed_m_s

class GannetPlus28(BaseAircraftClass):
    # Gannnet-Plus VTOL:
    def __init__(self):
        # super().__init__()
        self.aircraft_name         = "GannetVTOL-28 TTRSLT"
        self.aircraft_mass_kg      = 3.5
        self.batt_capacity_kJ      = 225 # 4s10Ah = 550
        self.batt_threshold_full   = 0.95
        self.batt_threshold_cutoff = 0.15
        self.solar_MPP_W           = 83.3
        self.power_consup_W        = 55.8
        self.cruise_airspeed_m_s   = 13.5
        self.climb_rate_m_s        = 0.6
        self.sink_rate_m_s         = -1.1 
        self.cruise_altitude_m     = 2000
        self.sprint_airspeed_m_s   = lambda x : np.sqrt(x*0.4)+self.cruise_airspeed_m_s

class MiniHawk(BaseAircraftClass):
    # MiniHawk-VTOL
    def __init__(self):
        # super().__init__()
        self.aircraft_name         = "MiniHawk-VTOL"
        self.aircraft_mass_kg      = 1.4
        self.batt_capacity_kJ      = 60
        self.batt_threshold_full   = 0.95
        self.batt_threshold_cutoff = 0.35  # Hover landing needs 0.35, 8.3km hops
        self.solar_MPP_W           = 8.0
        self.power_consup_W        = 120.7
        self.cruise_airspeed_m_s   = 18.0
        self.climb_rate_m_s        = 0.6
        self.sink_rate_m_s         = -1.1
        self.cruise_altitude_m     = 2000
        self.sprint_airspeed_m_s   = lambda x : np.sqrt(x * 0.4)+self.cruise_airspeed_m_s

class MAXiHawk16(BaseAircraftClass):
    # MAXiHawk-VTOL (Large MiniHawk, 16 Cells)
    def __init__(self):
        # super().__init__()
        self.aircraft_name         = "MAXiHawk16-NosePuller" # "MAXiHawk16-NosePuller" "MAXiHawk16-Stock"
        self.aircraft_mass_kg      = 2.8
        self.batt_capacity_kJ      = 180   # 4s3700mAh with reserve, 13 minutes at 140W
        self.batt_threshold_full   = 0.95
        self.batt_threshold_cutoff = 0.35  # 19km hops (13km more likely)
        self.solar_MPP_W           = 48.0  # 3.5*16*0.85
        self.power_consup_W        = 95.0 # 45.0/0.6=75, but who am I kidding, twin pullers optimized for hover ~= 140W, use 95.0 for Nose-Puller
        self.cruise_airspeed_m_s   = 14.5
        self.climb_rate_m_s        = 0.6
        self.sink_rate_m_s         = -1.1
        self.cruise_altitude_m     = 2000
        self.sprint_airspeed_m_s   = lambda x : np.sqrt(x * 0.4)+self.cruise_airspeed_m_s

    # MaxiHawk16_Elevons_-3deg T2-VLM2+Drag3 (Wetted @ 0.010 + 1.5dm^2 @ 1.2), CL=0.4 alpha=6, pitch moment is zero
    # Power Curve points:
    # Vx    Fx*Vx
    # 12    45
    # 13    55
    # 14    70
    # 15    82
    # 16    98
    # 18    144
    # 20    210

# Needs to be turned into an object
def flyDay(dayNumber, initialPosition, batt_energy_kJ, aircraft, graphDay=False, tailwind_m_s=2.0, ocean_drift_m_s=0.1, cloud_cover_density=0.6):
    # Essentials for simulating the day of flying
    lat_radians = initialPosition[0]*np.pi/180
    # dayNumber %= 365
    declination_radians = -1.0*np.arcsin(0.39779*np.cos((0.98565*np.pi/180)*(dayNumber+10)+(1.914*np.pi/180)*np.sin((0.98565*np.pi/180)*(dayNumber-2))))
    
    # Parameters
    startHour = 4   # 0400 in the morning
    endHour   = 28  # 0400 / 4AM next day
    step_hr   = 0.01 # 0.1 # Step every 6 minutes 
    
    # tailwind_m_s = 2.0 # 4.0
    # ocean_drift_m_s = 0.1

    # cloud_cover_density = 0.6
    cloud_cover_opacity = 0.05
    cloud_cover_op_dev  = 0.2

    cloud_attenuation = 1.0 # ersistant Variable, not a parameter
    
    
    # Recording Elements
    arr_times           = np.arange(startHour, endHour, step_hr)
    arr_power_flux_W    = np.zeros(arr_times.shape)
    arr_solar_power_W   = np.zeros(arr_times.shape)
    arr_batt_energy_kJ  = np.zeros(arr_times.shape)
    arr_altitude_m      = np.zeros(arr_times.shape)
    arr_airspeed_m_s    = np.zeros(arr_times.shape)
    # arr_distance_km     = np.zeros(arr_times.shape)
    
    # Runtime Variables for the ForLoop
    motorActiveState = False
    distanceFlown_km      = 0
    # batt_energy_kJ        = batt_capacity_kJ * batt_threshold_cutoff # Assume that the aircraft cut and glided into the night and entered sleep mode previously
    altitude_m            = 0
    prevAltitude_m        = 0
    airspeed_m_s          = 0
    hopsCount             = 0
    landedTime            = 0
    
    
    
    # for currentHour in arr_times: # range(startHour, stopHour, step_hr):
    # for i, hour in np.ndenumerate(arr_times): # Probably breaking tradition and culture by doing it this way
    for i, hour in enumerate(arr_times): # https://stackoverflow.com/a/49384823/6228436
        hourAngle_radians = (hour-12.0)/12*np.pi
        solarFactor = max(np.sin(lat_radians)*np.sin(declination_radians)+np.cos(lat_radians)*np.cos(declination_radians)*np.cos(hourAngle_radians),0)
        # solarAngle_deg = 90.0-np.arccos(solarFactor)*180/np.pi

        if(i%20 == 0):
            rv = random.random()
            cloud_attenuation = max(min((0.5-rv)*cloud_cover_op_dev + cloud_cover_opacity, 1), 0) if(rv < cloud_cover_density) else 1.0

        solarFactor *= cloud_attenuation
            
        
        # Regardless of if the aircraft is flying or not, panel is contributing something
        solar_power_W    = aircraft.solar_MPP_W * solarFactor
        arr_solar_power_W[i] = solar_power_W
        # solar_contrib_kJ = solar_power_W * 3.6 * step_hr
        
        # Find the power budget
        if(motorActiveState):
            if(altitude_m < aircraft.cruise_altitude_m):
                power_flux_W = (solar_power_W - aircraft.power_consup_W - aircraft.aircraft_mass_kg * 9.8 * aircraft.climb_rate_m_s)
            else:
                power_flux_W = (solar_power_W - aircraft.power_consup_W)
        else: 
            power_flux_W = solar_power_W
        arr_power_flux_W[i] = power_flux_W
        
        # Credit or Subtract however much altitude has been gained/lost (Needs to be driven by power_flux_W, assumed values for now)
        vertical_rate_m_s = aircraft.climb_rate_m_s if motorActiveState else aircraft.sink_rate_m_s
        altitude_m = min(max(altitude_m + vertical_rate_m_s*step_hr*3600,0),aircraft.cruise_altitude_m)
        arr_altitude_m[i] = altitude_m
        
        if(altitude_m<=0 and prevAltitude_m > 0): landedTime = hour
        prevAltitude_m = altitude_m
        
        # Calculate what our airspeed has been, determined after altitude so that airspeed forced to zero when altitude was/is zero
        if(altitude_m<=0):
            airspeed_m_s = ocean_drift_m_s
        elif(altitude_m < aircraft.cruise_altitude_m):
            airspeed_m_s = aircraft.cruise_airspeed_m_s + tailwind_m_s
        elif(power_flux_W > 0 and motorActiveState):
            if(batt_energy_kJ >= aircraft.batt_capacity_kJ):
                airspeed_m_s = aircraft.sprint_airspeed_m_s(power_flux_W) + tailwind_m_s
                power_flux_W *= 0.02 # We are consuming the extra power to cruise faster, the above lambda yields equilibrium for power budget
            else:
                airspeed_m_s = aircraft.cruise_airspeed_m_s + tailwind_m_s
        else:
            airspeed_m_s = aircraft.cruise_airspeed_m_s + tailwind_m_s # Default for gliding and cruise
        arr_airspeed_m_s[i] = airspeed_m_s
            
        # Distance Accumulated from the past interval up to this point
        distanceFlown_km += (airspeed_m_s) * 3.6 * step_hr
        # arr_distance_km(i) = distanceFlown_km
        
        
        # Credit or Subtract from Battery Energy State whatever has been gained or lost up to this point
        batt_energy_kJ += power_flux_W * step_hr * 3.6 # J/s * 0.1hr * 3600sec/hr *0.001 kJ/J
        batt_energy_kJ = aircraft.batt_capacity_kJ if(batt_energy_kJ > aircraft.batt_capacity_kJ) else 0 if(batt_energy_kJ < 0) else batt_energy_kJ
        arr_batt_energy_kJ[i] = batt_energy_kJ
        
        
        # Determine what the next motor state should be
        # If battery is above threshold for launch, we turn on motor. If we are already flying and above battery cut, we keep the motor on. Power flux not necessarily greater than power_consup_W at launch time.
        # if(solar_MPP_W >= power_consup_W):
        if(False):
            aboveLaunchThresholdState = (batt_energy_kJ >= (batt_capacity_kJ*batt_threshold_full)) and (power_flux_W >= power_consup_W)
        else:
            aboveLaunchThresholdState = (batt_energy_kJ >= (aircraft.batt_capacity_kJ*aircraft.batt_threshold_full))
            
        aboveCutThresholdState = (batt_energy_kJ >= (aircraft.batt_capacity_kJ*aircraft.batt_threshold_cutoff)) and motorActiveState # isFlyingState
        nextMotorActiveState = (aboveLaunchThresholdState or aboveCutThresholdState)
        if(not motorActiveState and nextMotorActiveState): hopsCount +=1
        motorActiveState = nextMotorActiveState
        # fluxPositiveState = solar_power_W > power_consup_W # Aircraft that have an MPP ratio below 1.0 will never meed this condition, can only do hops, no continuous daytime flight
        
        
        # print('{0:2d}:{1:4.1f},{2:0.2f}'.format(dayNumber, hour, batt_energy_kJ) )
        if(batt_energy_kJ <= 0): print('# Battery Emptied on {0:2d}:{1:4.1f},{2:0.2f}'.format(dayNumber, hour, batt_energy_kJ) )
        # print('{0:2d},{1:2f}'.format(int(i), distanceFlown_km) )
        
    # if(distanceFlown_km<=0):
    #     # print('# No-Fly on day {0:2d}, Lat: {1:0.4f}, Dec: {2:0.4f}'.format(dayNumber, lat_radians, declination_radians))
    #     bestSolarFactor = max(np.sin(lat_radians)*np.sin(declination_radians)+np.cos(lat_radians)*np.cos(declination_radians)*1.0,0)
    #     # print('# No-Fly on day {0:2d}, Solar: {1:0.4f}'.format(dayNumber, bestSolarFactor))
    
    # if(hopsCount > 0):
    if(graphDay):
        plt.plot(arr_times, arr_power_flux_W,   '.r',     label="Sys.InOut (W)") #axs[0]. instead of plt.
        plt.plot(arr_times, arr_solar_power_W,  color='orange', marker=',', linestyle='none',   label="Solar.In (W)")
        plt.plot(arr_times, arr_batt_energy_kJ*0.1, '-c', label="Battery (*10kJ)")
        plt.plot(arr_times, arr_altitude_m*0.01, '-k',    label="Altitude (*100m)")
        plt.plot(arr_times, arr_airspeed_m_s,   '-b',     label="Airspeed (m/s)")
        plt.legend()
        plt.suptitle(aircraft.aircraft_name + " - " + ( datetime.datetime(2023, 1, 1) + datetime.timedelta(dayNumber - 1) ).strftime("%d %b") )
        plt.title("{lat:0.4f},{lon:0.4f}".format(lat=initialPosition[0],lon=initialPosition[1]))
        plt.ylim(bottom= -10, top=100)
        plt.gca().set_xlabel('Local Solar Time')
        # plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        plt.gca().set_xticks([4, 8, 12, 16, 20, 24, 28])
        current_values = plt.gca().get_xticks()
        # plt.gca().set_xticklabels(['{:%H}'.format(x, x%60) for x in current_values])
        plt.gca().set_xticklabels( [(datetime.datetime(2023, 1, 1) + datetime.timedelta(dayNumber - 1 + (x/24.0) ) ).strftime("%H:%M") for x in current_values] )
        # plt.gca().set_xticklabels(['{:f}:'])
        # plt.gca().xaxis_date()
        plt.text(3, 80, "Distance Flown: {dist:.1f}km\nTailwind: {tail:0.1f}m/s\nSea Drift: {drift:0.1f}m/s\nCloud Cover: {clouds:0.0f}%".format(
            dist=distanceFlown_km, tail=tailwind_m_s, drift=ocean_drift_m_s, clouds=(cloud_cover_density*100)) )
        # plt.draw()

        # plt.show()

    if(False):
        axs[0].plot(arr_times, arr_power_flux_W,   '.r',     label="Sys.InOut (W)") #axs[0]. instead of plt.
        axs[0].plot(arr_times, arr_solar_power_W,  '--y',    label="Solar.In (W)")
        axs[0].plot(arr_times, arr_batt_energy_kJ*0.1, '-c', label="Battery (*10kJ)")
        axs[0].plot(arr_times, arr_altitude_m*0.01, '-k',    label="Altitude (*100m)")
        axs[0].plot(arr_times, arr_airspeed_m_s,   '-b',     label="Airspeed (m/s)")
        axs[0].legend()
        # plt.suptitle(aircraft.aircraft_name + " - Day-of-Year #" + str(dayNumber) )
        plt.suptitle(aircraft.aircraft_name + " - " + ( datetime.datetime(2023, 1, 1) + datetime.timedelta(dayNumber - 1) ).strftime("%d %b") )
        plt.title("{lat:0.4f},{lon:0.4f}".format(lat=initialPosition[0],lon=initialPosition[1]))
        # axs[0].set_title("Test")
        # axs[0].ylim(bottom= -5, top=100)
        # plt.yscale('symlog', linthreshy=100)
        axs[0].text(4, -3, "Distance Flown: {dist:.1f}km".format(dist=distanceFlown_km) )
        plt.draw()

        # if(graphDay): plt.show()
        # plt.show()

        # plt.show(block=False)
        # plt.savefig

    # if(True):
    #     m.drawgreatcircle(-119.3598, 37.155, -155.9159, 29.696,del_s=50,color='red', lw=1.)
    #     plt.show()
    
    
    return distanceFlown_km, hopsCount, batt_energy_kJ, landedTime
    

class Equator():
    name = "Equator"
    waypoints = [
        [ 00.000000000,    0.0000000000, 20.0, 0.0, 0.0, 0.0],
        [ 00.000000000,   10.0000000000, 20.0, 0.0, 0.0, 0.1],
        [ 00.000000000,   20.0000000000, 20.0, 0.0, 0.0, 0.2],
        [ 00.000000000,   30.0000000000, 20.0, 0.0, 0.0, 0.3],
        [ 00.000000000,   40.0000000000, 20.0, 0.0, 0.0, 0.4],
        [ 00.000000000,   50.0000000000, 20.0, 0.0, 0.0, 0.5],
        [ 00.000000000,   60.0000000000, 20.0, 0.0, 0.0, 0.6],
        [ 00.000000000,   70.0000000000, 20.0, 0.0, 0.0, 0.7],
        [ 00.000000000,   80.0000000000, 20.0, 0.0, 0.0, 0.8],
        [ 00.000000000,   90.0000000000, 20.0, 0.0, 0.0, 0.9],
        [ 00.000000000,  100.0000000000, 20.0, 0.0, 0.0, 1.0]]

class TropicCancer():
    name = "TropicCancer"
    waypoints = [
        [ 23.436806000, -120.0000000000, 20.0],
        [ 23.436806000, -130.0000000000, 20.0]]

class Parallel39():
    name = "Parallel39"
    waypoints = [
        [ 39.000000000,    0.0000000000, 20.0, 0.0, 0.0, 0.0],
        [ 39.000000000,   10.0000000000, 20.0, 0.0, 0.0, 0.1],
        [ 39.000000000,   20.0000000000, 20.0, 0.0, 0.0, 0.2],
        [ 39.000000000,   30.0000000000, 20.0, 0.0, 0.0, 0.3],
        [ 39.000000000,   40.0000000000, 20.0, 0.0, 0.0, 0.4],
        [ 39.000000000,   50.0000000000, 20.0, 0.0, 0.0, 0.5],
        [ 39.000000000,   60.0000000000, 20.0, 0.0, 0.0, 0.6],
        [ 39.000000000,   70.0000000000, 20.0, 0.0, 0.0, 0.7],
        [ 39.000000000,   80.0000000000, 20.0, 0.0, 0.0, 0.8],
        [ 39.000000000,   90.0000000000, 20.0, 0.0, 0.0, 0.9],
        [ 39.000000000,  100.0000000000, 20.0, 0.0, 0.0, 1.0]]

class CaliforniaHawaii():
    name = "CA-to-HI"
    waypoints = [
        [ 37.155236000, -122.3598450000, 20.0, 3.0, 0.2, 0.6], # CA
        [ 20.696066103, -155.9159486517, 20.0, 6.0, 0.2, 0.4]] # HI

class RenoDubai():
    name = "Reno-to-Dubai"
    waypoints = [
        [ 39.5, -119.8, 20.0],
        [ 25.2, 55.2, 20.0]]

class Circumnavigation():
    name = "Circumnavigation"
    # Lat, Lon, WpRadius, tailwind, oceanCurrent, cloudCover
    waypoints = [
        [ 20.696066103, -155.9159486517, 20.0, 6.0, 0.2, 0.5],
        [-10.792915986,  156.7953091010, 20.0, 2.0, 0.2, 0.7],
        [ -9.947399142,  132.7279547619, 20.0, 2.0, 0.2, 0.7],
        [-15.564947689,   77.2859058055, 20.0, 6.0, 0.2, 0.5],
        [-21.499054523,   55.6506165861, 20.0, 4.0, 0.2, 0.5],
        [-35.527079059,   20.5090348490, 20.0, 0.0, 0.2, 0.4],
        [ 10.086734865,  -30.7684334790, 20.0, 4.0, 0.2, 0.5],
        [ 12.752370940,  -74.1045676950, 20.0, 4.0, 0.2, 0.4],
        [  6.314611730,  -79.8605277598, 20.0, 4.0, 0.2, 0.4],
        [ 20.696066103, -155.9159486517, 20.0, 6.0, 0.2, 0.6]]



# Lat, Long, Radius_km
# waypoints = [
    
    # Seattle to Boston
    # [ 47.6, -122.3, 20.0],
    # [ 42.4, -71.1, 20.0]]
    
    # Latitude ZigZag
    # [ 39.5, -119.8, 20.0],
    # [ 39.5, -129.8, 20.0]]
    # [ 39.5, -119.8, 20.0],
    # [ 39.5, -129.8, 20.0],
    # [ 39.5, -119.8, 20.0],
    # [ 39.5, -129.8, 20.0]]

    # Lake Tahoe Laps
    # [ 39.230, -119.970, 2.0],
    # [ 39.150, -120.030, 2.0],
    # [ 39.000, -120.030, 2.0],
    # [ 39.230, -119.970, 2.0],
    # [ 39.150, -120.030, 2.0],
    # [ 39.000, -120.030, 2.0],
    # [ 39.230, -119.970, 2.0],
    # [ 39.150, -120.030, 2.0],
    # [ 39.000, -120.030, 2.0],
    # [ 39.230, -119.970, 2.0],
    # [ 39.150, -120.030, 2.0],
    # [ 39.000, -120.030, 2.0]]
    

# mission = Equator
# mission = Parallel39
mission = CaliforniaHawaii
# mission = RenoDubai
# mission = Circumnavigation

# aircraftList = [GannetSUWAVE(), LaysianQuadplane(), MAXiHawk16()]
# aircraftList = [MAXiHawk16(), GannetSUWAVE(), LaysianQuadplane()]
# aircraftList = [LaysianQuadplane(), MAXiHawk16()]
aircraftList = [LaysianPlus28(), GannetPlus28(), MAXiHawk16()]
# aircraftList = [MAXiHawk16()]
# aircraftList = [MiniHawk()]

startingDayList = range(160, 161)  # Mid-June, typically for Hawaii
# startingDayList = range(172, 173)  # Summer Solstice
# startingDayList = range(240, 241)  # Late-August
# startingDayList = range(315, 316)  # Mid-November, typically for circumnav
# startingDayList = range(15, 16)  # Mid-January
# startingDayList = range(79, 80)  # Spring Equinox, typically for Equator Baseline test
# startingDayList = range(0,365,15)
# startingDayList = range(0,30,1)



graphDay = True
# saveVideo = False
staticDate = False


bestDayOffset = 0
bestDaysCount = 1000

waypoints = mission.waypoints

for aircraft in aircraftList:

    # aircraft()
    # print(type(aircraft))
    if(len(waypoints)<2): sys.exit("Need at least two waypoints")

    print("# " + aircraft.aircraft_name)

    start_day_str = ( datetime.datetime(2023, 1, 1) + datetime.timedelta(startingDayList[0]) ).strftime("%d %b")
    video_filename = str(mission.name)+"_"+start_day_str+"_"+str(aircraft.aircraft_name)+".mp4" # "_map.mp4"
    video_filename = "-".join( video_filename.split())

    metadata = dict(title=str(mission.name), artist='Matplotlib', comment='Migratory Mission Path Planner')
    writer = FFMpegWriter(fps=2, metadata=metadata)
    fig = plt.figure()
    # fig, axs = plt.subplots(2, 1)
    # fig.subplots_adjust(hspace=0.5)

    # world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    # world.plot(ax=axs[1])

    llcorner = [ min([x[0] for x in mission.waypoints]), min([x[1] for x in mission.waypoints]) ]
    urcorner = [ max([x[0] for x in mission.waypoints]), max([x[1] for x in mission.waypoints]) ]
    center = [ (llcorner[0]+urcorner[0])/2, (llcorner[1]+urcorner[1])/2 ]
    marginFactor = 1.2
    span = [ (urcorner[0]-llcorner[0])*marginFactor/2, (urcorner[1]-llcorner[1])*marginFactor/2 ]
    # print(llcorner, urcorner)
    # print(span)
    llcorner = [ (center[0]-span[0]), (center[1]-span[1]) ]
    urcorner = [ (center[0]+span[0]), (center[1]+span[1]) ]
    # print(llcorner, urcorner)

    # sys.exit("Stuff that needs to be solved but won't: Find Great Circle max/min latitude passed, ")

    # m = Basemap(llcrnrlat=llcorner[0], llcrnrlon=llcorner[1], urcrnrlat=urcorner[0], urcrnrlon=urcorner[1], resolution='l', ax=axs[1])
    # m.drawcoastlines()

    random.seed(12345)
    

    with writer.saving(fig, video_filename, 300):
    # if True:
        for dayOffset in startingDayList:
            dayNumber = dayOffset # Day of year
            waypointIndex = 1

            prevWaypoint = waypoints[0]
            nextWaypoint = waypoints[waypointIndex]
            
            missionDaysCount = 1
            distanceTraveled_km = 0
            totalHopsCount = 0
            segmentFlown_km = 0
            segmentTotal_km = greatCircleDistance_km(prevWaypoint, nextWaypoint)

            routeArcs = []
            # m.drawcoastlines()
            
            batt_energy_kJ = 1
            bestDayDistance_km = 0

            if(graphDay): print('date, waypointIndex, dayFlown_km, hopsCount, distanceTraveled_km, landedTime')

            # fig, axs = plt.subplots(2, 1)


            while(waypointIndex < len(waypoints)):
                
                # Solve for the current lat and long given previously accumulated segment distance
                fraction = segmentFlown_km / segmentTotal_km
                startingPosition = intermediatePoint(prevWaypoint, nextWaypoint, fraction)
                routeArcs.append(startingPosition)

                # print(routeArcs[-2][1],routeArcs[-2][0])

                # if(len(routeArcs) > 1): m.drawgreatcircle(routeArcs[-2][1],routeArcs[-2][0],routeArcs[-1][1],routeArcs[-1][0],del_s=5,color='red', lw=1.)
                
                tailwind_m_s        = prevWaypoint[3] * (1.0-fraction) + nextWaypoint[3] * fraction
                ocean_drift_m_s     = prevWaypoint[4] * (1.0-fraction) + nextWaypoint[4] * fraction
                cloud_cover_density = prevWaypoint[5] * (1.0-fraction) + nextWaypoint[5] * fraction
                
                # Roll the current day and get the distance traveled
                dayFlown_km, hopsCount, batt_energy_kJ, landedTime = flyDay(dayNumber, startingPosition, batt_energy_kJ, aircraft, graphDay, 
                    tailwind_m_s, ocean_drift_m_s, cloud_cover_density) # , axs=axs)
                # ani = animation.FuncAnimation(fig, animate, len(y), interval=dt*1000, blit=True)

                if(graphDay):                
                    writer.grab_frame()
                    plt.clf() # For Graphs, disable for map

                totalHopsCount += hopsCount
                
                bestDayDistance_km = max(bestDayDistance_km,dayFlown_km)
                
                # Update distance accumulations
                segmentFlown_km += dayFlown_km
                distanceTraveled_km += dayFlown_km
                
                # Check if we have passed the current waypoint, and increment if within range
                if((segmentFlown_km >= segmentTotal_km) or (np.abs(segmentTotal_km-segmentFlown_km)<nextWaypoint[2])):
                    # Should probably use greatCircleDistance_km, but extra computing cost?
                    waypointIndex +=1
                    segmentFlown_km = max((segmentFlown_km-segmentTotal_km),0)
                    
                    if(waypointIndex >=len(waypoints)): break
                    prevWaypoint = nextWaypoint # waypoints[waypointIndex-1]
                    nextWaypoint = waypoints[waypointIndex]

                    # Find the length of this segment for calculating when we have completed it
                    segmentTotal_km = greatCircleDistance_km(prevWaypoint, nextWaypoint)
                
                # At this point, we have completed a day of flying, and have potentially passed the last waypoint and are on the next one.
                # print('{0:5d},{1:5d},{2:10.1f},{3:10.1f},{4:12.8f},{5:10.8f}'.format(dayNumber, waypointIndex, dayFlown_km, distanceTraveled_km, startingPosition[0], startingPosition[1]) )
                # print('{0:5d},{1:5d},{2:10.1f},{3:3d},{4:10.1f}'.format(dayNumber, waypointIndex, dayFlown_km, hopsCount, distanceTraveled_km) )

                date_str = ( datetime.datetime(2023, 1, 1) + datetime.timedelta(dayNumber - 1) ).strftime("%d %b")
                if(graphDay): print('{0},{1:5d},{2:10.1f},{3:3d},{4:10.1f},{5:5.1f}'.format(date_str, waypointIndex, dayFlown_km, hopsCount, distanceTraveled_km, landedTime) )
                
                # Increment the Day and repeat
                dayNumber = (dayNumber +1) if(staticDate == False) else (dayNumber)
                missionDaysCount +=1
                
            # bestDayOffset = dayOffset if(missionDaysCount<bestDaysCount)
            # daysList.append([dayOffset, missionDaysCount, totalHopsCount, bestDayDistance_km])
            print('Day Offset: {0:4d}, {1:3d} days, {2:3d} hops, Best Day: {3:0.2f}km'.format(dayOffset,missionDaysCount,totalHopsCount,bestDayDistance_km))
            # print("")

        print("")