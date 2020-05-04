import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from country_aliases import country_aliases

# If you have saved a local copy of the CSV file as LOCAL_CSV_FILE,
# set READ_FROM_URL to True
READ_FROM_URL = True

# Start the plot on the day when the number of confirmed cases reaches MIN_CASES
MIN_CASES = 1

# Plot for MAX_DAYS days after the day on which each country reaches MIN_CASES.
MAX_DAYS = 40

#PLOT_TYPE = 'deaths'
PLOT_TYPE = 'confirmed cases'

# These are the GitHub URLs for the Johns Hopkins data in CSV format.
if PLOT_TYPE == 'confirmed cases':
    data_loc = ('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/'
                'master/csse_covid_19_data/csse_covid_19_time_series/'
                'time_series_covid19_confirmed_global.csv')
    LOCAL_CSV_FILE = 'covid-19-cases.csv'
elif PLOT_TYPE == 'deaths':
    data_loc = ('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/'
                'master/csse_covid_19_data/csse_covid_19_time_series/'
                'time_series_covid19_deaths_global.csv')
    LOCAL_CSV_FILE = 'covid-19-deaths.csv'

# Read in the data to a pandas DataFrame.
if not READ_FROM_URL:
    data_loc = LOCAL_CSV_FILE


df = pd.read_csv(data_loc)
df.rename(columns={'Country/Region': 'Country'}, inplace=True)

# Read in the populations file as a Series (squeeze=True) indexed by country.
populations = pd.read_csv('country_populations.csv', index_col='Country',
                          squeeze=True)

# Group by country and sum over the different states/regions of each country.
grouped = df.groupby('Country')
df2 = grouped.sum()

#print("df2: Group sum: " + str(df2))

df2.rename(index=country_aliases, inplace=True)

pd.set_option('display.max_rows', 2000)
#print("df2: Rename(index=country_aliases: " + str(df2))

def make_bar_plot(country):
    """Make the bar plot of case numbers and change in numbers line plot."""

    # Extract the Series corresponding to the case numbers for country.
    c_df = df2.loc[country, df2.columns[3:]].astype(int)
    # Convert index to a proper datetime object
    c_df.index = pd.to_datetime(c_df.index)
    # Discard rows before the number reaches the threshold MIN_CASES.
    c_df = c_df[c_df >= MIN_CASES]
    n = len(c_df)
    if n == 0:
        print('Too few data to plot: minimum number of {}s is {}'
                .format(PLOT_TYPE, MIN_CASES))
        sys.exit(1)

    #fig = plt.Figure()
    fig = plt.figure(figsize=(10, 8))

    # Arrange the subplots on a grid: the top plot (case number change) is
    # one quarter the height of the bar chart (total confirmed case numbers).
    ax2 = plt.subplot2grid((4,1), (0,0))
    ax1 = plt.subplot2grid((4,1), (1,0), rowspan=3)
    ax1.bar(range(n), c_df.values)
    # Force the x-axis to be in integers (whole number of days) in case
    # Matplotlib chooses some non-integral number of days to label).
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

    c_df_change = c_df.diff()
    ax2.bar(range(n), c_df_change.values, fc='tab:green')
    ax2.set_xticks([])

    ax1.set_xlabel('Days since {} {}'.format(MIN_CASES, PLOT_TYPE))
    ax1.set_ylabel(f'Number of {PLOT_TYPE}, $N$')
    ax2.set_ylabel('$\Delta N$')

    # Add a title reporting the latest number of cases available.
    title = '{}\n{} {} on {}'.format(country, c_df[-1], PLOT_TYPE,
                c_df.index[-1].strftime('%d %B %Y'))
    plt.suptitle(title)


def make_comparison_plot(countries, normalize=False):
    """Make a plot comparing case numbers in different countries."""
     
    #print("countries: " + str(countries))
    
    pd.set_option('display.max_rows', 2000)
    
    #print("df2: " + str(df2.sort_values(by=['5/2/20'], ascending=False)))
    
    # Extract the Series corresponding to the case numbers for countries.
    #df2.fillna(0, inplace=True)
    #df2.replace(np.nan, 0, inplace=True)
    #ddf2 = df2.fillna(0)
    
    c_df = df2.loc[countries, df2.columns[3:]].astype(int)
    
    ##### print("c_df: Extract the Series " + str(c_df))

    # Discard any columns with fewer than MIN_CASES.
    c_df = c_df[c_df >= MIN_CASES]
    
    ##### print("c_df: Discard any columns " + str(c_df))

    if normalize:
        # Calculate confirmed case numbers per 1,000,000 population.
        
        c_df  = c_df.div(populations.loc[countries], axis='index') * 1000000
        
        
        ##### print("c_df: normalize " + str(c_df))

    # Rearrange DataFrame to give countries in columns and number of days since
    # MIN_CASES in rows.
    c_df = c_df.T.apply(lambda e: pd.Series(e.dropna().values))
    ##### print("c_df:  Rearrange DataFrame " + str(c_df))

    # Truncate the DataFrame after the maximum number of days to be considered.
    c_df = c_df.truncate(after=MAX_DAYS-1)
    ##### print("c_df:  Truncate the DataFrame " + str(c_df))
    
    # Plot the data.
    #fig = plt.figure()
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot()
    for country, ser in c_df.iteritems():
        ax.plot(range(len(ser)), np.log10(ser.values), label=country)

    if not normalize:
        # Set the tick marks and labels for the absolute data.
        ymin = int(np.log10(MIN_CASES))
        ymax = int(np.log10(np.nanmax(c_df))) + 1
        yticks = np.linspace(ymin, ymax, ymax-ymin+1, dtype=int)
        yticklabels = [str(10**y) for y in yticks]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        ax.set_ylim(ymin, ymax)
        ax.set_ylabel(f'Number of {PLOT_TYPE}')
    else:
        # Set the tick marks and labels for the per 1,000,000 population data.
        ax.set_ylim(np.log10(np.nanmin(c_df)), np.log10(np.nanmax(c_df)))
        ax.set_ylabel(f'Number of {PLOT_TYPE} per 1,000,000 population')

    # Label the x-axis
    ax.set_xlim(0, MAX_DAYS)
    ax.set_xlabel(f'Number of days since first {MIN_CASES} {PLOT_TYPE}')
    ax.set_xlabel(f'Number of days since first {MIN_CASES} {PLOT_TYPE}')
    
    



    def plot_threshold_lines(doubling_lifetime):
        """Add a line for the growth in numbers at a given doubling lifetime."""

        # Find the limits of the line for the current plot region.
        x = np.array([0, MAX_DAYS])
        y = np.log10(MIN_CASES) + x/doubling_lifetime * np.log10(2)
        ymin, ymax = ax.get_ylim()
        if y[1] > ymax:
            y[1] = ymax
            x[1] = doubling_lifetime/np.log10(2) * (y[1] - np.log10(MIN_CASES))
        ax.plot(x, y, ls='--', color='#aaaaaa')

        # The reason this matters is that we want to label the line at its
        # centre, rotated appropriately.
        s = f'Doubling every {doubling_lifetime} days'
        p1 = ax.transData.transform_point((x[0], y[0]))
        p2 = ax.transData.transform_point((x[1], y[1]))
        xylabel = ((x[0]+x[1])/2, (y[0]+y[1])/2)

        dy = (p2[1] - p1[1])
        dx = (p2[0] - p1[0])
        angle = np.degrees(np.arctan2(dy, dx))
        ax.annotate(s, xy=xylabel, ha='center', va='center', rotation=angle)

    if not normalize:
        # If we're plotting absolute numbers, indicate the doubling time.
        plot_threshold_lines(2)
        plot_threshold_lines(3)
        plot_threshold_lines(5)

    ax.legend()


# The country to plot the data for.
#country = 'Austria'
#country = 'United States'
#country = 'Iran'
#country = 'United Kingdom'
#country = 'China'
#make_bar_plot(country)
#plt.show()



countries4 = ['Vatican City']


countries10 = ['United States', 'Spain', 'Italy', 'United Kingdom', 
              'France', 'Germany', 'Turkey', 'Russia']

countries11 = ['Brazil', 'China', 'Canada', 
              'Belgium', 'Peru', 'Netherlands' , 'India']

countries12 = ['Switzerland', 'Ecuador', 'Saudi Arabia', 
              'Portugal', 'Mexico', 'Sweden' , 'Ireland', 'Pakistan']

countries13 = ['Chile', 'Singapore', 'Israel', 'Belarus', 'Austria',
              'Qatar', 'Japan', 'United Arab Emirates' , 'Poland', 'Romania']

countries14 = ['Ukraine', 'Indonesia', 'South Korea', 'Denmark', 
              'Philippines', 'Bangladesh', 'Norway' , 'Czech Republic', 'Dominican Republic']

countries15 = ['Colombia', 'Panama', 'Australia', 'South Africa', 
              'Egypt', 'Malaysia', 'Finland' , 'Morocco', 'Argentina']


make_comparison_plot(countries4,  normalize=False)
plt.show()

make_comparison_plot(countries10,  normalize=False)
plt.show()

make_comparison_plot(countries11,  normalize=False)
plt.show()

make_comparison_plot(countries12,  normalize=False)
plt.show()

make_comparison_plot(countries13,  normalize=False)
plt.show()

make_comparison_plot(countries14,  normalize=False)
plt.show()



