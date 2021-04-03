

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

def Pldtarget(target, i, targetstring):
    pixelfile = lk.search_targetpixelfile(target)[i].download();
    lc = pixelfile.to_lightcurve(method="pld").remove_outliers().flatten()
    period = lc.to_periodogram("bls").period_at_max_power
    lc.fold(period).scatter();
    plt.title(str(target) + " (" + str(i) + ")" + str(targetstring));

# %%


##################################  003  #########################################
# %%
print("\nSection 003\n")

import ipywidgets as widgets
from IPython.display import display, clear_output

r1_select_variable = widgets.Dropdown(
     options=[

'Kepler-1','Kepler-2','Kepler-3','Kepler-4','Kepler-5','Kepler-6','Kepler-7','Kepler-8','Kepler-9',
'Kepler-10','Kepler-11','Kepler-12','Kepler-13','Kepler-14','Kepler-15','Kepler-16','Kepler-17','Kepler-18','Kepler-19',
'Kepler-20','Kepler-21','Kepler-22','Kepler-23','Kepler-24','Kepler-25','Kepler-26','Kepler-27','Kepler-28','Kepler-29',
'Kepler-30','Kepler-31','Kepler-32','Kepler-33','Kepler-34','Kepler-35','Kepler-36','Kepler-37','Kepler-38','Kepler-39',
'Kepler-40','Kepler-41','Kepler-42','Kepler-43','Kepler-44','Kepler-45','Kepler-46','Kepler-47','Kepler-48','Kepler-49',
'Kepler-50','Kepler-51','Kepler-52','Kepler-53','Kepler-54','Kepler-55','Kepler-56','Kepler-57','Kepler-58','Kepler-59',
'Kepler-60','Kepler-61','Kepler-62','Kepler-63','Kepler-64','Kepler-65','Kepler-66','Kepler-67','Kepler-68','Kepler-69',
'Kepler-70','Kepler-71','Kepler-72','Kepler-73','Kepler-74','Kepler-75','Kepler-76','Kepler-77','Kepler-78','Kepler-79',
'Kepler-80','Kepler-81','Kepler-82','Kepler-83','Kepler-84','Kepler-85','Kepler-86','Kepler-87','Kepler-88','Kepler-89',
'Kepler-90','Kepler-91','Kepler-92','Kepler-93','Kepler-94','Kepler-95','Kepler-96','Kepler-97','Kepler-98','Kepler-99',

"KIC 3832474", "Kepler-30",
"KIC 4548011", "Kepler-1581",
"KIC 4852528", "Kepler-80",
"KIC 5972334", "Kepler-487",
'KIC 6922244', 'Kepler-8',
'KIC 8462852', 
"KIC 8480285", "Kepler-647",
"KIC 10264202",
"KIC 10337517", "Kepler-783",
"KIC 11030475",
"KIC 11442793", "Kepler-90",
"KIC 11568987", "Kepler-534",

"Proxima Cen",
"Tau Ceti",
"Trappist-1",
"Wolf 359",

         
    'Choose a target'],
    value='Choose a target',
    description='Target:',
    disabled=False,
)
def get_variable(b):
    clear_output
    print(r1_select_variable.value)
    
display(r1_select_variable)


# %%


##################################  004  #########################################
# %%
print("\nSection 004\n")

# Print a list of all target entries

target1 = r1_select_variable.value
print("Target: " + str(target1))

search_result = lk.search_lightcurve(target1)
print(search_result)
print("\n")

Target_Name = " "

# %%


##################################  004B  #########################################
# %%
print("\nSection 004B\n")

# Print the first entry
print("\nThe first entry\n")

Pldtarget(target1, 0, Target_Name)

# %%


##################################  004C  #########################################
# %%
print("\nSection 004C\n")

# Print the graphs of the first 5 entries
print("\nGraphs of the first 5 entries\n")

for i in range(0, 5):
    Pldtarget(target1, i, Target_Name)

# %%


#################################  004D  #########################################
# %%
print("\nSection 004D\n")

# Print all the entries
print("\nGraphs of all the entries\n")

for i in range(0, int(len(search_result))):
    Pldtarget(target1, i, Target_Name)

# %%



##################################  005  #########################################
# %%
print("\nSection 005\n")

# Print the graphs of the first half of the entries
print("\nGraphs of the first half of the entries\n")

for i in range(0, int(len(search_result)/2)):
    Pldtarget(target1, i, Target_Name)

# %%


##################################  006  #########################################
# %%
print("\nSection 006\n")

# Print the graphs of the second half of the entries
print("\nGraphs of the second half of the entries\n")

for i in range(int(len(search_result)/2) -1, len(search_result)):
    Pldtarget(target1, i, Target_Name)

# %%


##################################  007  #########################################
# %%
print("\nSection 007\n")

# Print the graphs of the first 10 entries
print("\nGraphs of the first 10 entries\n")

for i in range(0, 10):
    Pldtarget(target1, i, Target_Name)

# %%


##################################  008  #########################################
# %%
print("\nSection 008\n")

# Print the 5th entry
print("\nGraphs of the f5th entry\n")

Pldtarget(target1, 5, Target_Name)

# %%


##################################  009  #########################################
# %%
print("\nSection 009\n")

# https://github.com/openexoplanetcatalogue/open_exoplanet_catalogue/

import xml.etree.ElementTree as ET, urllib.request, gzip, io
url = "https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems.xml.gz"
oec = ET.parse(gzip.GzipFile(fileobj=io.BytesIO(urllib.request.urlopen(url).read())))

# %%


##################################  010  #########################################
# %%
print("\nSection 010\n")
# Output mass and radius of all planets 
for planet in oec.findall(".//planet"):
    print([planet.findtext("mass"), planet.findtext("radius")])

# %%


##################################  011  #########################################
# %%
print("\nSection 011\n")

# Find all circumbinary planets 
for planet in oec.findall(".//binary/planet"):
    print(planet.findtext("name"))

# %%


##################################  012  #########################################
# %%
print("\nSection 012\n")

# Output distance to planetary system (in pc, if known) and number of planets in system
for system in oec.findall(".//system"):
    print(system.findtext("distance"), len(system.findall(".//planet")))

# %%


#################################  013 #########################################
# %%
print("\nSection 013\n")

# Find all circumbinary planets 
for star in oec.findall(".//star"):
    print(star.findtext("star"))

# %%