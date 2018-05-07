# README
Last updated by Stephanie Labou on November 17, 2017

---

# V1 - Kara Woo

This code is the beginnings of preparing the long-term Lake Baikal plankton data
for analysis with
[MARSS](https://cran.r-project.org/web/packages/MARSS/index.html). This analysis
will ultimately incorporate data on *Epischura baikalensis*, temperature,
chlorophyll, and phytoplankton. The data is not in this repo. The data lives in
a Subversion repository that you'll need access to in order to run the code. In
the future it would probably be better to symlink the data into the repository,
but I haven't done that here.

## What has been done so far?

Mostly a lot of exploring various aspects of the data to help guide decisions
about what data to use. See the `epi_data_for_marss.R` script for the code that
does everything described below.

### Epischura data

The data for zooplankton has codes to identify species. Certain codes have
previously been identified as "double counts" -- often these are codes that
represent the sum of various other codes in the data (as a fictional example,
species code D might be defined as the sum of codes A, B, and C), and so if
these species codes are included in the analysis we will be double counting some
individuals. I have removed the double count data.

In the Epischura data, counts are given as **individuals x 1000 / m^2**. To
convert to individuals per liter, I use the following function:

``` R
m2_to_l <- function(x, interval) {
  stopifnot(is.numeric(x))                           # count should be numeric
  ## takes units in 1000 individuals/m2 and converts to individuals per liter
  individuals <- x * 1000                            # convert to individuals/m2
  count_per_liter <- individuals / (interval * 1000) # convert to indiv./liter
  return(count_per_liter)
}
```

### Temperature and chlorophyll data

I have aggregated the temperature data into monthly averages.

### Understanding depth

An irritatingly tricky part of this process is that the different types of data
(zooplankton, temperature, chlorophyll, etc.) do not get sampled at the same
depths. I have made some histograms that show the distributions of maximum
sampling depth for Epischura, temperature, and chlorophyll.

I also created a bar chart showing average Epischura counts and average
chlorophyll by depth.

Finally, I've written some code to look at which depths are most frequently
sampled, as this might help us choose which depth intervals to work with.

### Phytoplankton data

I have attempted to identify some of the most abundant phytoplankton taxa in
winter and summer.

Due to a change in algae preservation technique, only phytoplankton data from
1975 and later should be included.

I separated unidentified phytoplankton taxa into picoplankton, nanoplankton, and
unknown based on information in the original phytoplankton key
(`Phytoplankton_Key_20120725.xls`, which is located in the Subversion repository
containing the data under `Longterm_data/phyto/Phytoplankton_key/`).

For each date and depth, I found the genera that constitute >=10% of
phytoplankton abundance. Then, I created a list of all genera that at some point
constitute >=10% of abundance. I plotted the most common genera in summer (July,
August, September) and winter (I first defined winter as January, February, and
March, but then tried defining it as Feb/Mar/Apr/May, which might be more
appropriate per discussions with Steph H).

Finally, I did a little work to see which unidentified taxa are most abundant.
The vast majority are bacteria.

## Known issues

* The unit conversion for Epischura counts introduced a few errors because
  either the upper and lower layers in the interval were swapped (leading to
  negative counts) or the same (leading to infinite counts).

## Next steps

* Decide which depth layers of data to use. These should be the deepest depths
  that still typically have all relevant data types (temperature, chlorophyll,
  etc.).
* Decide on which phytoplankton taxa to include in the model.

# V2 - Stephanie Labou
	
Depth resolution: Data has been limited to 0-50m (inclusive). For temperature, I have included the average temperature 0-50m. For chla, I have both the average and sum total for 0-50m. For plankton, values are sums within 0-50m. Zooplankton are provided from depth ranges, so reported values are total zoop summed over 0-10, 10-25, and 25-50m. 
Temporal resolution: Values were aggregated to month within years. Monthly data was then subset to only periods when all data is available: 1979-1999 (primarily limited by phytos, which only go through 1999).

NAs replacement: 
* Single NA: retained
* Two or more sequential NA: replaced by monthly average for that variable across whole time period

Replacing zeros:
* For each variable, I replaced zeros with half the minimum non-zero value for that value, plus a random number between 0 and 1. 
* Note: each variable received a different random number and random number remained consistent within variable. Therefore, within a variable, all zero values became the same non-zero value.
* Note: this also means that some cases, an NA value was appropriately replaced with a zero (monthly avg) which was then replaced appropriately with a non-zero. Should be fine, but somewhat squishy, so making aware.

