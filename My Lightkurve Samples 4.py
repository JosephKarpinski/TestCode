

##################################  001  #########################################
# %%
print("\nSection 001\n")
import matplotlib.pyplot as plt
%matplotlib inline

import astropy.units as u

# %%


##################################  001B  #########################################
# %%
print("\nSection 001B\n")
from importlib.metadata import version

import lightkurve as lk
print(version('lightkurve'))

from lightkurve import search_targetpixelfile
from lightkurve import search_lightcurve

# %%


##################################  002  #########################################
# %%
print("\nSection 002\n")

def Periodogram(target, targetquarter):
    lc = search_lightcurve(target, author="Kepler", quarter=targetquarter, cadence="long").download().remove_nans()
    
    pg = lc.to_periodogram(oversample_factor=10)
    pg.plot();
    plt.title("Periodogram Oversample Factor [10]");

    pg = lc.to_periodogram(oversample_factor=1)
    pg.plot();
    plt.title("Periodogram Oversample Factor [1]");

    pg.plot(view='period', scale='log');
    plt.title("Log Scale Period Space");

    lc.fold(period=pg.period_at_max_power).scatter();
    plt.title("Period at Max Power");

    lc.fold(period=2*pg.period_at_max_power, wrap_phase=0.2).scatter();
    plt.title("Plot With Double the Period");

    lc.fold(period=4*pg.period_at_max_power, wrap_phase=0.2).scatter();
    plt.title("Plot With Quadruple the Period");

    print("\nPeriod at Max Power\n" + str(pg.period_at_max_power) + "\n")
    


def FlatTarget(target, targetquarter, targetexp):
    pixelfile = search_targetpixelfile(target, quarter=targetquarter, exptime=targetexp).download();
    
    lc = pixelfile.to_lightcurve(aperture_mask=pixelfile.pipeline_mask);
    lc.scatter();

    flat_lc = lc.flatten(window_length=401);
    flat_lc.plot();
    plt.title("Flat")
    

def FoldTarget(target, targetquarter, targetexp, targetperiod):
    pixelfile = search_targetpixelfile(target, quarter=targetquarter, exptime=targetexp).download();
    #print("\nDownloaded")

    pixelfile.plot(frame=1);

    lc = pixelfile.to_lightcurve(aperture_mask=pixelfile.pipeline_mask);
    #plt.title("Scatter Plot");
    lc.scatter();
    

    lc.plot();
    plt.title("Light Curve")
    #print("\nMask Done!")

    """lc = pixelfile.to_lightcurve(aperture_mask='all');
    lc.plot();
    print("\nAll Done!")"""

    flat_lc = lc.flatten(window_length=401);
    flat_lc.plot();
    plt.title("Flat")
    #print("\nFlat Done!")

    folded_lc = flat_lc.fold(period=targetperiod);
    folded_lc.plot();
    plt.title("Folded")
    #print("\nFold Done!")

    binned_lc = folded_lc.bin(time_bin_size=0.01);
    binned_lc.plot();
    plt.title("Binned")
    #print("\nBinn Done!")

    

def ProcessTarget(target, targetquarter, targetexp):
    pixelfile = search_targetpixelfile(target, quarter=targetquarter, exptime=targetexp).download();
    #print("\nDownloaded")

    pixelfile.plot(frame=1);

    lc = pixelfile.to_lightcurve(aperture_mask=pixelfile.pipeline_mask);
    lc.scatter()
    
    lc.plot();
    #print("\nMask Done!")

    """lc = pixelfile.to_lightcurve(aperture_mask='all');
    lc.plot();
    print("\nAll Done!")"""

    flat_lc = lc.flatten(window_length=401);
    flat_lc.plot();
    #print("\nFlat Done!")
    
# %%


###################################  003  #########################################
# %%
print("\nSection 003\n")

t = """
Tabby's Star (also known as Boyajian's Star and WTF Star, and designated KIC 8462852 
in the Kepler Input Catalog) is an F-type main-sequence star in the constellation Cygnus 
approximately 1,470 light-years (450 pc) from Earth. 
Unusual light fluctuations of the star, including up to a 22% dimming in brightness.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 8462852', author='Kepler')
print(search_result)
print("\n")

ProcessTarget("KIC 8462852", 1, 1800, )

Periodogram("KIC 8462852", 1)

print("\nThe plot reveals a short-lived 20% dip in the brightness of the star. It looks like we re-discovered one of the intriguing dips in Tabby's star.\n")


pixelfile = lk.search_targetpixelfile("KIC 8462852")[1].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power
lc.fold(period).scatter();
plt.title("KIC 8462852 (1)");

for i in range(0,18):
    pixelfile = lk.search_targetpixelfile("KIC 8462852")[i].download();
    lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
    period = lc.to_periodogram("bls").period_at_max_power
    lc.fold(period).scatter();
    plt.title("KIC 8462852 (" + str(i) +")");

# %%



##################################  004  #########################################
# %%
print("\nSection 004\n")

# https://docs.lightkurve.org/tutorials/1-getting-started/what-are-lightcurve-objects.html

# http://www.openexoplanetcatalogue.com/planet/Kepler-8%20b/

search_result = lk.search_lightcurve('KIC 6922244', author='Kepler')
print(search_result)
print("\n")

FoldTarget("KIC 6922244", 4, 1800, 3.5225)
plt.title("KIC 6922244 (Kepler-8) Period (3.5225)");

Periodogram("KIC 6922244", 10)

lc = search_lightcurve('KIC 6922244', author="Kepler", quarter=10, cadence="long").download().remove_nans()
pg = lc.to_periodogram(minimum_period=3.3*u.day, maximum_period=3.7*u.day, oversample_factor=10)
print("\nRevised Period at Max Power\n" + str(pg.period_at_max_power) + "\n")

lc.fold(period=pg.period_at_max_power,  wrap_phase=0.2).scatter();
plt.title("Revised Period at Max Power");

pixelfile = lk.search_targetpixelfile("KIC 6922244")[5].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power
lc.fold(period).scatter();
plt.title("KIC 6922244 (5) (Kepler-8) Period (3.5225)");

pixelfile = lk.search_targetpixelfile("KIC 6922244")[10].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power
lc.fold(period).scatter();
plt.title("KIC 6922244 (10)");
lc.fold(period).truncate(-0.2, 0.5).scatter();
plt.title("KIC 6922244 (10) (Kepler-8) Period (3.5225)");

pixelfile = lk.search_targetpixelfile("KIC 6922244")[13].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power
lc.fold(period).scatter();
plt.title("KIC 6922244 (13) (Kepler-8) Period (3.5225)");

pixelfile = lk.search_targetpixelfile("KIC 6922244")[17].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power
lc.fold(period).scatter();
plt.title("KIC 6922244 (17) (Kepler-8) Period (3.5225)");

pixelfile = lk.search_targetpixelfile("KIC 6922244")[21].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power
lc.fold(period).scatter();
plt.title("KIC 6922244 (21) (Kepler-8) Period (3.5225)");

pixelfile = lk.search_targetpixelfile("KIC 6922244")[23].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power
lc.fold(period).scatter();
plt.title("KIC 6922244 (23) (Kepler-8) Period (3.5225)");

pixelfile = lk.search_targetpixelfile("KIC 6922244")[24].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power
lc.fold(period).scatter();
plt.title("KIC 6922244 (24) (Kepler-8) Period (3.5225)");

pixelfile = lk.search_targetpixelfile("KIC 6922244")[28].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power
lc.fold(period).scatter();
plt.title("KIC 6922244 (28) (Kepler-8) Period (3.5225)");

pixelfile = lk.search_targetpixelfile("KIC 6922244")[32].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power
lc.fold(period).scatter();
plt.title("KIC 6922244 (32) (Kepler-8) Period (3.5225)");

pixelfile = lk.search_targetpixelfile("KIC 6922244")[36].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power
lc.fold(period).scatter();
plt.title("KIC 6922244 (36) (Kepler-8) Period (3.5225)");

pixelfile = lk.search_targetpixelfile("KIC 6922244")[40].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power
lc.fold(period).scatter();
plt.title("KIC 6922244 (40) (Kepler-8) Period (3.5225)");

for i in range(44,49):
    pixelfile = lk.search_targetpixelfile("KIC 6922244")[i].download();
    lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
    period = lc.to_periodogram("bls").period_at_max_power
    lc.fold(period).scatter();
    plt.title("KIC 6922244 (" + str(i) +")" + " (Kepler-8) Period (3.5225)");


# %%


##################################  004B  #########################################
# %%
print("\nSection 004B\n")
# https://docs.lightkurve.org/tutorials/1-getting-started/what-are-periodogram-objects.html

search_result = lk.search_lightcurve('KIC 10264202', author='Kepler')
print(search_result)
print("\n")

FoldTarget("KIC 10264202", 4, 1800, 1.035)
plt.title("KIC 10264202 Period (1.035)");

Periodogram("KIC 10264202", 10)

lc = search_lightcurve('KIC 10264202', author="Kepler", quarter=10, cadence="long").download().remove_nans()
pg = lc.to_periodogram(minimum_period=0.9*u.day, maximum_period=1.2*u.day, oversample_factor=10)
print("\nRevised Period at Max Power\n" + str(pg.period_at_max_power) + "\n")

lc.fold(period=pg.period_at_max_power,  wrap_phase=0.2).scatter();
plt.title("Revised Period at Max Power");


pixelfile = lk.search_targetpixelfile("KIC 10264202")[10].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power

lc.fold(period).scatter();
plt.title("KIC 10264202 (10)");

lc.fold(period).truncate(-0.3, 0.1).scatter();
plt.title("KIC 10264202 (10)");

# %%


# https://arxiv.org/pdf/1712.05044.pdf


##################################  005  #########################################
# %%
print("\nSection 005\n")

# http://www.openexoplanetcatalogue.com/planet/Kepler-90%20b/

t = """
KIC 11442793 (Kepler-90). This system hosts seven transit- ing planets, ranging from slightly larger than Earth
to Jupiter-sized planets, all of which have been ei- ther validated or confirmed via transit-timing variations 
(Cabrera et al. 2014; Schmitt et al. 2014; Lissauer et al. 2014). 
We detect an additional TCE with a period of 14.4 days with S/N of 8.7."""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 11442793', author='Kepler')
print(search_result)
print("\n")

FoldTarget("KIC 11442793", 15, 1800, 14.4)
plt.title("KIC 11442793 (Kepler-90) Period (14.4)");

Periodogram("KIC 11442793", 10)

lc = search_lightcurve('KIC 11442793', author="Kepler", quarter=10, cadence="long").download().remove_nans()
pg = lc.to_periodogram(minimum_period=14.2*u.day, maximum_period=14.6*u.day, oversample_factor=10)
print("\nRevised Period at Max Power\n" + str(pg.period_at_max_power) + "\n")

lc.fold(period=pg.period_at_max_power,  wrap_phase=0.2).scatter();
plt.title("Revised Period at Max Power");


pixelfile = lk.search_targetpixelfile("KIC 11442793")[1].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power

lc.fold(period).scatter();
plt.title("KIC 11442793 (1)");

# %%


##################################  006  #########################################
# %%
print("\nSection 006\n")

# http://www.openexoplanetcatalogue.com/planet/Kepler-647%20b/

t = """
KIC 8480285 (Kepler-647). This system hosts two planets - one previously validated Earth-sized planet 
(Morton et al. 2016) and one candidate super-Earth in periods of 16.2 days and 29.7 days respectively. 
We detect an additional TCE at S/N = 8.4 with a period of 10.9 days."""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 8480285', author='Kepler')
print(search_result)
print("\n")

FoldTarget("KIC 8480285", 15, 1800, 16.2)
plt.title("KIC 8480285 Period (16.2)");

FoldTarget("KIC 8480285", 15, 1800, 29.7)
plt.title("KIC 8480285 Period (29.7)");

FoldTarget("KIC 8480285", 15, 1800, 10.9)
plt.title("KIC 8480285 (Kepler-647) Period (10.9)");

Periodogram("KIC 8480285", 10)

lc = search_lightcurve('KIC 8480285', author="Kepler", quarter=10, cadence="long").download().remove_nans()
pg = lc.to_periodogram(minimum_period=10.7*u.day, maximum_period=11.1*u.day, oversample_factor=10)
print("\nRevised Period at Max Power\n" + str(pg.period_at_max_power) + "\n")

lc.fold(period=pg.period_at_max_power,  wrap_phase=0.2).scatter();
plt.title("Revised Period at Max Power");


pixelfile = lk.search_targetpixelfile("KIC 8480285")[10].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power

lc.fold(period).scatter();
plt.title("KIC 8480285 (10)");

lc.fold(period).truncate(-6.5, -4.5).scatter();
plt.title("KIC 8480285 (10)");

# %%


##################################  007  #########################################
# %%
print("\nSection 007\n")

# http://www.openexoplanetcatalogue.com/planet/Kepler-534%20b/

t = """
KIC 11568987 (Kepler-534). This system hosts one validated super-Earth (Morton et al. 2016) and one candidate Earth-sized planet in 16 and 7.3 day periods, respectively. 
We detect a new TCE with a period of 27.1 days, S/N = 9.8 and a depth corresponding to a roughly Earth-sized planet.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 11568987', author='Kepler')
print(search_result)
print("\n")

FoldTarget("KIC 11568987", 15, 1800, 16)
plt.title("KIC 11568987 Period (16)");

FoldTarget("KIC 11568987", 15, 1800, 7.3)
plt.title("KIC 11568987 Period (7.3)");

FoldTarget("KIC 11568987", 15, 1800, 27.1)
plt.title("KIC 11568987 (Kepler-534) Period (27.1)");

Periodogram("KIC 11568987", 10)


pixelfile = lk.search_targetpixelfile("KIC 11568987")[10].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power

lc.fold(period).scatter();
plt.title("KIC 11568987 (10)");

"""lc = search_lightcurve('KIC 11568987', author="Kepler", quarter=15, cadence="long").download().remove_nans()
pg = lc.to_periodogram(minimum_period=26.9*u.day, maximum_period=27.3*u.day, oversample_factor=10)
print("\nRevised Period at Max Power\n" + str(pg.period_at_max_power) + "\n")

lc.fold(period=pg.period_at_max_power,  wrap_phase=0.2).scatter();
plt.title("Revised Period at Max Power");

"""
# %%


##################################  008  #########################################
# %%
print("\nSection 008\n")

# http://www.openexoplanetcatalogue.com/planet/Kepler-80%20b/

t = """
KIC 4852528 (Kepler-80). This system hosts five confirmed small planets, the outer four of which are in a rare dynamical configuration with the middle three planets and the outer three planets each locked in a three-body resonance (Lissauer et al. 2014; Mac- Donald et al. 2016). We find a new TCE with a period of 14.64 days at S/N = 8.6 exterior to all five previously known planets.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 4852528', author='Kepler')
print(search_result)
print("\n")

FoldTarget("KIC 4852528", 15, 1800, 14.64)
plt.title("KIC 4852528 (Kepler-80) Period (14.64)");

Periodogram("KIC 4852528", 8)

lc = search_lightcurve('KIC 4852528', author="Kepler", quarter=8, cadence="long").download().remove_nans()
pg = lc.to_periodogram(minimum_period=14.44*u.day, maximum_period=14.84*u.day, oversample_factor=10)
print("\nRevised Period at Max Power\n" + str(pg.period_at_max_power) + "\n")

lc.fold(period=pg.period_at_max_power,  wrap_phase=0.2).scatter();
plt.title("Revised Period at Max Power");


pixelfile = lk.search_targetpixelfile("KIC 4852528")[8].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power

lc.fold(period).scatter();
plt.title("KIC 4852528 (8)");

lc.fold(period).truncate(-2.8, -1.8).scatter();
plt.title("KIC 4852528 (8)");

# %%


##################################  009  #########################################
# %%
print("\nSection 009\n")

# http://www.openexoplanetcatalogue.com/planet/Kepler-487%20b/

t = """
KIC 5972334 (Kepler-487). This system has two validated planets (Morton et al. 2016) and two candidates, including a warm Jupiter, 
an ultra-short-period planet (Sanchis-Ojeda et al. 2014), and a super- Earth in between, 
giving it a distinct architecture reminiscent of the 55 Cnc 
(Nelson et al. 2014) and WASP-47 (Becker et al. 2015) planetary systems. 
We detect a strong (S/N = 10.2) TCE with a period of 6.02 days which, if real, 
would have a radius a bit larger than that of the Earth. 
We see extra scatter in the light curve during the transit of this new TCE, however, 
which needs to be better understood before this TCE can be considered a good candidate.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 5972334', author='Kepler')
print(search_result)
print("\n")

FoldTarget("KIC 5972334", 15, 1800, 6.02)
plt.title("KIC 5972334 (Kepler-487) Period (6.02)");

Periodogram("KIC 5972334", 12)

lc = search_lightcurve('KIC 5972334', author="Kepler", quarter=12, cadence="long").download().remove_nans()
pg = lc.to_periodogram(minimum_period=5.8*u.day, maximum_period=6.2*u.day, oversample_factor=10)
print("\nRevised Period at Max Power\n" + str(pg.period_at_max_power) + "\n")

lc.fold(period=pg.period_at_max_power,  wrap_phase=0.2).scatter();
plt.title("Revised Period at Max Power");


pixelfile = lk.search_targetpixelfile("KIC 5972334")[1].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power

lc.fold(period).scatter();
plt.title("KIC 5972334 (1)");

lc.fold(period).truncate(-8, -6).scatter();
plt.title("KIC 5972334 (1)");

# %%


##################################  010  #########################################
# %%
print("\nSection 010\n")

# http://www.openexoplanetcatalogue.com/planet/Kepler-783%20b/

t = """
KIC 10337517 (Kepler-783). This system hosts a sub-Earth- sized validated planet in a 4.3 day orbit (Morton et al. 2016) and a super-Earth-sized planet can- didate in a 7 day orbit. We detect an additional TCE in an 11.1 day orbit, which also would likely be sub-Earth-sized.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 10337517', author='Kepler')
print(search_result)
print("\n")

FoldTarget("KIC 10337517", 15, 1800, 4.3)
plt.title("KIC 10337517 Period (4.3)");

FoldTarget("KIC 10337517", 15, 1800, 7)
plt.title("KIC 10337517 Period (7)");

FoldTarget("KIC 10337517", 15, 1800, 11.1)
plt.title("KIC 10337517 (Kepler-783) Period (11.1)");

Periodogram("KIC 10337517", 10)

lc = search_lightcurve('KIC 10337517', author="Kepler", quarter=10, cadence="long").download().remove_nans()
pg = lc.to_periodogram(minimum_period=10.9*u.day, maximum_period=11.3*u.day, oversample_factor=10)
print("\nRevised Period at Max Power\n" + str(pg.period_at_max_power) + "\n")

lc.fold(period=pg.period_at_max_power,  wrap_phase=0.2).scatter();
plt.title("Revised Period at Max Power");


pixelfile = lk.search_targetpixelfile("KIC 10337517")[1].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power

lc.fold(period).scatter();
plt.title("KIC 10337517 (1)");

# %%


##################################  011  #########################################
# %%
print("\nSection 011\n")

t = """
KIC 11030475. This system has 4 planet candidates, none of which have been validated. We detect a new TCE at an orbital period of 4.75 days and S/N = 8.7, with the same time of transit, depth, duration, and exactly half the period of a previously detected candidate, KOI 2248.02. This “new” TCE, therefore, appears to be the true pe- riod of the candidate KOI 2248.02.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 11030475', author='Kepler')
print(search_result)
print("\n")

FoldTarget("KIC 11030475", 15, 1800, 4.75)
plt.title("KIC 11030475 Period (4.75)");

Periodogram("KIC 11030475", 10)


pixelfile = lk.search_targetpixelfile("KIC 11030475")[1].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power
lc.fold(period).scatter();
plt.title("KIC 11030475 (1)");

"""lc = search_lightcurve('KIC 111030475', author="Kepler", quarter=15, cadence="long").download().remove_nans()
pg = lc.to_periodogram(minimum_period=4.55*u.day, maximum_period=4.95*u.day, oversample_factor=10)
print("\nRevised Period at Max Power\n" + str(pg.period_at_max_power) + "\n")

lc.fold(period=pg.period_at_max_power,  wrap_phase=0.2).scatter();
plt.title("Revised Period at Max Power");"""

# %%


##################################  012  #########################################
# %%
print("\nSection 012\n")

# http://www.openexoplanetcatalogue.com/planet/Kepler-1581%20b/

t = """
KIC 4548011 (Kepler-1581). This system hosts one validated planet (Morton et al. 2016) in a 6.3 day period and one additional candidate in a 9.1 day period. We detect a new TCE with a period of 13.95 days with a relatively strong S/N = 9.7, which if real, would likely be a sub-Earth-sized planet. Subse- quently, this signal was also detected by the Kepler pipeline and included in the DR25 release of the Kepler planet candidate catalog as KOI 4288.04.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 4548011', author='Kepler')
print(search_result)
print("\n")

FoldTarget("KIC 4548011", 15, 1800, 6.3)
plt.title("KIC 4548011 Period (6.3)");

FoldTarget("KIC 4548011", 15, 1800, 9.1)
plt.title("KIC 4548011 Period (9.1)");

FoldTarget("KIC 4548011", 15, 1800, 13.95)
plt.title("KIC 4548011 (Kepler-1581) Period (13.95)");

Periodogram("KIC 4548011", 10)

lc = search_lightcurve('KIC 4548011', author="Kepler", quarter=10, cadence="long").download().remove_nans()
pg = lc.to_periodogram(minimum_period=13.75*u.day, maximum_period=14.15*u.day, oversample_factor=10)
print("\nRevised Period at Max Power\n" + str(pg.period_at_max_power) + "\n")

lc.fold(period=pg.period_at_max_power,  wrap_phase=0.2).scatter();
plt.title("Revised Period at Max Power");

# %%


##################################  013  #########################################
# %%
print("\nSection 013\n")

# 

t = """
KIC 3832474 (Kepler-30). This system hosts two confirmed giant planets, and a smaller confirmed super-Earth- sized planet called Kepler-30 b, which orbits with a period of 29.16 days (Fabrycky et al. 2012). All of these planets show transit timing variations (TTVs), and Kepler-30 d shows extremely large TTVs with an amplitude of more than a day (Fab- rycky et al. 2012; Panichi et al. 2017). We detect a new TCE with a period of 29.45 days - very close to the mean period of Kepler-30 d. We suspect this TCE is due to residual transits of Kepler-30 d that we were unable to completely remove from the light curve due to the large TTVs.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 3832474', author='Kepler')
print(search_result)
print("\n")

FoldTarget("KIC 3832474", 15, 1800, 29.16)
plt.title("KIC 3832474 Period (29.16)");

FoldTarget("KIC 3832474", 15, 1800, 29.45)
plt.title("KIC 3832474 (Kepler-30) Period (29.45)");

Periodogram("KIC 3832474", 10)

"""lc = search_lightcurve('KIC 3832474', author="Kepler", quarter=7, cadence="long").download().remove_nans()
pg = lc.to_periodogram(minimum_period=29.25*u.day, maximum_period=29.65*u.day, oversample_factor=10)
print("\nRevised Period at Max Power\n" + str(pg.period_at_max_power) + "\n")

lc.fold(period=pg.period_at_max_power,  wrap_phase=0.2).scatter();
plt.title("Revised Period at Max Power");"""


# %%


##################################  014  #########################################
# %%
print("\nSection 014\n")

search_result = lk.search_lightcurve('Tau Ceti')
print(search_result)
print("\n")

lk.search_lightcurve('Tau Ceti')[0].download().plot();
plt.title("Tau Ceti (0)");

lk.search_lightcurve('Tau Ceti')[1].download().plot();
plt.title("Tau Ceti (1)");

lk.search_lightcurve('Tau Ceti')[2].download().plot();
plt.title("Tau Ceti (2)");

lk.search_lightcurve('Tau Ceti')[3].download().plot();
plt.title("Tau Ceti (3)");

lk.search_lightcurve('Tau Ceti')[4].download().plot();
plt.title("Tau Ceti (4)");

lk.search_lightcurve('Tau Ceti')[5].download().plot();
plt.title("Tau Ceti (5)");

# %%


##################################  015  #########################################
 # %%
print("\nSection 015\n")

search_result = lk.search_lightcurve('Proxima Cen')
print(search_result)
print("\n")

lk.search_lightcurve('Proxima Cen')[0].download().plot();
plt.title("Proxima Cen (0)");

lk.search_lightcurve('Proxima Cen')[1].download().plot();
plt.title("Proxima Cen (1)");

lk.search_lightcurve('Proxima Cen')[2].download().plot();
plt.title("Proxima Cen (2)");

lk.search_lightcurve('Proxima Cen')[3].download().plot();
plt.title("Proxima Cen (3)");


pixelfile = lk.search_targetpixelfile("Proxima Cen")[1].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power

lc.fold(period).scatter();
plt.title("Proxima Cen");


# %%


##################################  016  #########################################
# %%
print("\nSection 016\n")

search_result = lk.search_lightcurve('Wolf 359')
print(search_result)
print("\n")

lk.search_lightcurve('Wolf 359')[0].download().plot();
plt.title("Wolf 359 (0)");
lc = lk.search_lightcurve('Wolf 359')[0].download();
flat_lc = lc.flatten(window_length=401);
flat_lc.plot();
plt.title("Flat");

lk.search_lightcurve('Wolf 359')[1].download().plot();
plt.title("Wolf 359 (1)");

#FlatTarget("Wolf 359", 1, 1800)

lc = lk.search_lightcurve('Wolf 359')[1].download();
flat_lc = lc.flatten(window_length=401);
flat_lc.plot();
plt.title("Flat");


pixelfile = lk.search_targetpixelfile("Wolf 359")[1].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power

lc.fold(period).scatter();
plt.title("Wolf 359");


# %%


##################################  017  #########################################
# %%
print("\nSection 017\n")

# http://www.openexoplanetcatalogue.com/planet/TRAPPIST-1%20b/

search_result = lk.search_targetpixelfile("Trappist-1");
print(search_result)
print("\n")

pixelfile = lk.search_targetpixelfile("Trappist-1")[1].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power

lc.fold(period).scatter();
plt.title("Trappist-1");

lc.fold(period).truncate(-0.5, -0.25).scatter();
plt.title("Trappist-1");

# %%


##################################  017B  #########################################
# %%
print("\nSection 017B\n")

# http://www.openexoplanetcatalogue.com/planet/TRAPPIST-1%20b/

search_result = lk.search_targetpixelfile("Trappist-1b");
print(search_result)
print("\n")

pixelfile = lk.search_targetpixelfile("Trappist-1b")[1].download();
lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
period = lc.to_periodogram("bls").period_at_max_power

lc.fold(period).scatter();
plt.title("Trappist-1b");

lc.fold(period).truncate(-0.5, -0.25).scatter();
plt.title("Trappist-1b");

# %%


##################################  018  #########################################
# %%
print("\nSection 018\n")

# http://www.openexoplanetcatalogue.com/planet/TRAPPIST-1%20b/

t = """
Trappist-1
"""

print(t + "\n")

search_result = lk.search_lightcurve('Trappist-1')
print(search_result)
print("\n")

lk.search_lightcurve('Trappist-1')[0].download().plot();
plt.title("Trappist-1 (0)");
lc = lk.search_lightcurve('Trappist-1')[0].download();
flat_lc = lc.flatten(window_length=401);
flat_lc.plot();
plt.title("Flat");


FoldTarget("Trappist-1", 0, 1800, 1.51)
plt.title("Trappist-1 Period (1.51)");

#Periodogram("Trappist-1", 0)

#lc = search_lightcurve(target, author="Kepler", quarter=targetquarter, cadence="long").download().remove_nans()
    
pg = lc.to_periodogram(oversample_factor=10)
pg.plot();
plt.title("Periodogram Oversample Factor [10]");

pg = lc.to_periodogram(oversample_factor=1)
pg.plot();
plt.title("Periodogram Oversample Factor [1]");

pg.plot(view='period', scale='log');
plt.title("Log Scale Period Space");

lc.fold(period=pg.period_at_max_power).scatter();
plt.title("Period at Max Power");

lc.fold(period=2*pg.period_at_max_power, wrap_phase=0.2).scatter();
plt.title("Plot With Double the Period");

lc.fold(period=4*pg.period_at_max_power, wrap_phase=0.2).scatter();
plt.title("Plot With Quadruple the Period");

print("\nPeriod at Max Power\n" + str(pg.period_at_max_power) + "\n")




"""
search_result = lk.search_lightcurve('Trappist-1', author='K2')
print(search_result)
print("\n")

FoldTarget("Trappist-1", 1, 1800, 1.51)
plt.title("Trappist-1 Period (1.51)");

Periodogram("KIC 5972334", 1)

lc = search_lightcurve('Trappist-1', author="K2", quarter=1, cadence="long").download().remove_nans()
pg = lc.to_periodogram(minimum_period=1.3*u.day, maximum_period=1.7*u.day, oversample_factor=10)
print("\nRevised Period at Max Power\n" + str(pg.period_at_max_power) + "\n")

lc.fold(period=pg.period_at_max_power,  wrap_phase=0.2).scatter();
plt.title("Revised Period at Max Power");
"""




# %%

