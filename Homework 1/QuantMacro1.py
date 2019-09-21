import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Loading the data from harddrive
url1 = "C:/Users/Tino/Desktop/Table 1.1.5.xls" # table including GDP data
url2 = 'C:/Users/Tino/Desktop/Table 1.12.xls' # table including T-S, NMI, CE
url3 = 'C:/Users/Tino/Desktop/Table 5.2.5.xls' # table including IPP
url4 = 'C:/Users/Tino/Desktop/Table 1.13.xls' # table including figures for Corporate Business Income, Corporate CE, Corporate TS, Corporate NMI
url5 = 'C:/Users/Tino/Desktop/Table 1.10.xls' # table including Net Operational Surplus
url6 = 'C:/Users/Tino/Desktop/Table FA 1.1.xls' # table including Net Fixed Assets Constant Costs
# check 6.2 and 7.2 for net fixed assets

data1 = pd.read_excel(url1)
data2 = pd.read_excel(url2)
data3 = pd.read_excel(url3)
data4 = pd.read_excel(url4)
data5 = pd.read_excel(url5)
data6 = pd.read_excel(url6)

# Picking the appropriate series' from the datasets
# Variables for part 1 and 2
GDP_annual = data1.iloc[0:1, 2:92]
tax_subsidies_annual = data2.iloc[18:20, 2:92]
NMI_annual = data2.iloc[8:9, 2:92]
CE_annual = data2.iloc[1:2, 2:92]
IPP_annual_private = data3.iloc[18:19, 2:92]
IPP_annual_public = data3.iloc[46:47, 2:92]
IPP_annual = pd.concat([IPP_annual_private, IPP_annual_public], axis=0) # merging the rows for the two IPP measures

# Variables for part 3
Corp_inc = data4.iloc[2:3,2:73]
Corp_CE = data4.iloc[3:4,2:73]
Corp_NMI = data4.iloc[7:8, 2:73]
Corp_TS = data4.iloc[8:9, 2:73]

# Variables for part 4
Prop_inc = data2.iloc[8:9, 32:92] # shortened series because observations for Net Operating surplus start in 1959
NOS = data5.iloc[8:9, 32:92]
Fixed_assets = data6.iloc[1:2, 36:96]


# vectorizing the data frames
t_s_annual_vec = pd.DataFrame(tax_subsidies_annual).to_numpy()
t_s_annual_vec = -1*np.diff(t_s_annual_vec, axis=0) # subtracting both rows to get the taxes minus subsidies measure
t_s_annual_vec = t_s_annual_vec[0]
IPP_annual_vec = pd.DataFrame(IPP_annual).to_numpy()
IPP_annual_vec = np.sum(IPP_annual_vec, axis=0) # summing the rows to get the measure for public+private IPP
CE_annual_vec = pd.DataFrame(CE_annual).to_numpy()
NMI_annual_vec = pd.DataFrame(NMI_annual).to_numpy()
GDP_annual_vec =pd.DataFrame(GDP_annual).to_numpy()

Corp_CE_vec = pd.DataFrame(Corp_CE).to_numpy()
Corp_inc_vec = pd.DataFrame(Corp_inc).to_numpy()
Corp_TS_vec =pd.DataFrame(Corp_TS).to_numpy()
Corp_NMI_vec = pd.DataFrame(Corp_NMI).to_numpy()
Corp_IPP = data3.iloc[18:19, 21:92] # shorten series since other data only available from 1948
Corp_IPP_vec = pd.DataFrame(Corp_IPP).to_numpy()

Prop_inc_vec = pd.DataFrame(Prop_inc).to_numpy()
NOS_vec = pd.DataFrame(NOS).to_numpy()
Fixed_assets_vec = pd.DataFrame(Fixed_assets).to_numpy()

### Question 1
## Exercise 1
# Calculating ratios
TS_GDP = np.divide(t_s_annual_vec, GDP_annual_vec)[0]
NMI_GDP =np.divide(NMI_annual_vec, GDP_annual_vec)[0]
IPP_GDP =np.divide(IPP_annual_vec, GDP_annual_vec)[0]

# Plotting ratios
x = [1929, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2018]

TS_GDP, = plt.plot(TS_GDP, label ='TS/GDP')
NMI_GDP, = plt.plot(NMI_GDP, label ='NMI/GDP')
IPP_GDP, = plt.plot(IPP_GDP, label ='IPP/GDP')
plt.xlabel('Year')
plt.ylabel('Percentage')
plt.title('Exercise 1.1 ratio time series\'')
plt.xticks(np.arange(0, 92, 10), x)
plt.legend(handles = [TS_GDP, NMI_GDP, IPP_GDP], loc='upper right')
plt.show()

## Exercise 2
# Calculating the ratios given in a), b), c)
LS_0 = np.divide(CE_annual_vec, GDP_annual_vec)[0]
LS_1 = np.divide(CE_annual_vec, GDP_annual_vec-t_s_annual_vec)[0]
LS_2 = np.divide(CE_annual_vec, GDP_annual_vec-t_s_annual_vec-NMI_annual_vec)[0]

# Plotting the ratios
CE_GDP, = plt.plot(LS_0, label ='CE/GDP')
CE_GDP_TS, = plt.plot(LS_1, label ='CE/(GDP-TS)')
CE_GDP_TS_NMI, = plt.plot(LS_2, label ='CE/(GDP-TS-NMI)')
plt.xlabel('Year')
plt.ylabel('Percentage')
plt.title('Exercise 1.2 ratio time series\'')
plt.xticks(np.arange(0, 92, 10), x)
plt.legend(handles = [CE_GDP, CE_GDP_TS, CE_GDP_TS_NMI], loc='upper right')
plt.show()

### Question 2
LS_0_2 = np.divide(CE_annual_vec, GDP_annual_vec-IPP_annual_vec)[0]
LS_1_2 = np.divide(CE_annual_vec, GDP_annual_vec-t_s_annual_vec-IPP_annual_vec)[0]
LS_2_2 = np.divide(CE_annual_vec, GDP_annual_vec-t_s_annual_vec-NMI_annual_vec-IPP_annual_vec)[0]

# Plotting the ratios
CE_GDP_2, = plt.plot(LS_0_2, label ='CE/(GDP-IPP)')
CE_GDP_TS_2, = plt.plot(LS_1_2, label ='CE/(GDP-TS-IPP)')
CE_GDP_TS_NMI_2, = plt.plot(LS_2_2, label ='CE/(GDP-TS-NMI-IPP)')
plt.xlabel('Year')
plt.ylabel('Percentage')
plt.title('Question 2 ratio time series\'')
plt.xticks(np.arange(0, 92, 10), x)
plt.legend(handles = [CE_GDP_2, CE_GDP_TS_2, CE_GDP_TS_NMI_2], loc='upper right')
plt.show()

### Question 3
# Calculating the ratios, this time using only variables on the corporate level
LS_0_3 = np.divide(Corp_CE_vec, Corp_inc_vec)[0]
LS_1_3 = np.divide(Corp_CE_vec, Corp_inc_vec-Corp_TS_vec)[0]
LS_2_3 = np.divide(Corp_CE_vec, Corp_inc_vec-Corp_TS_vec-Corp_NMI_vec)[0]

LS_0_4 = np.divide(Corp_CE_vec, Corp_inc_vec-Corp_IPP_vec)[0]
LS_1_4 = np.divide(Corp_CE_vec, Corp_inc_vec-Corp_TS_vec-Corp_IPP_vec)[0]
LS_2_4 = np.divide(Corp_CE_vec, Corp_inc_vec-Corp_TS_vec-Corp_NMI_vec-Corp_IPP_vec)[0]

# Plotting the ratis
xa =[1948, 1960, 1970, 1980, 1990, 2000, 2010, 2018]
Corp_CE_Inc, = plt.plot(LS_0_3, label ='Corp CE/Inc')
Corp_CE_Inc_TS, = plt.plot(LS_1_3, label ='Corp CE/(Inc-TS)')
Corp_CE_Inc_TS_NMI, = plt.plot(LS_2_3, label ='Corp CE/(Inc-TS-NMI)')
plt.xlabel('Year')
plt.ylabel('Percentage')
plt.title('Question 3 ratio time series\'')
plt.xticks(np.arange(0, 73, 10), xa)
plt.legend(handles = [Corp_CE_Inc, Corp_CE_Inc_TS, Corp_CE_Inc_TS_NMI], loc='upper left')
plt.show()

Corp_CE_Inc_2, = plt.plot(LS_0_4, label ='Corp CE/(Inc-IPP)')
Corp_CE_Inc_TS_2, = plt.plot(LS_1_4, label ='Corp CE/(Inc-TS-IPP)')
Corp_CE_Inc_TS_NMI_2, = plt.plot(LS_2_4, label ='Corp CE/(Inc-TS-NMI-IPP)')
plt.xlabel('Year')
plt.ylabel('Percentage')
plt.title('Question 3 ratio time series\'')
plt.xticks(np.arange(0, 73, 10), xa)
plt.legend(handles = [Corp_CE_Inc_2, Corp_CE_Inc_TS_2, Corp_CE_Inc_TS_NMI_2], loc='upper left')
plt.show()

## Question 4
# Calculating the capital share for each labor share and SNA between 1959 and 2018
alpha1 = LS_0[30:90] # capital share according to naive measure SNA2008
alpha2 = LS_1[30:90] # capital share adjusted for taxes and subsidies according to SNA2008
alpha3 = LS_2[30:90] # capital share adjusted for taxes, subsidies and NMI according to SNA2008
alpha4 = LS_0_2[30:90] # capital share according to naive measure SNA1993 proxy
alpha5 = LS_1_2[30:90] # capital share adjusted for taxes and subsidies according to SNA1993 proxy
alpha6 = LS_2_2[30:90] # capital share adjusted for taxes, subsidies and NMI according to SNA1993 proxy

ROC1 = np.divide(NOS_vec - alpha1*Prop_inc_vec,Fixed_assets_vec)[0]
ROC2 = np.divide(NOS_vec - alpha2*Prop_inc_vec,Fixed_assets_vec)[0]
ROC3 = np.divide(NOS_vec - alpha3*Prop_inc_vec,Fixed_assets_vec)[0]

ROC4 = np.divide(NOS_vec - alpha4*Prop_inc_vec,Fixed_assets_vec)[0]
ROC5 = np.divide(NOS_vec - alpha5*Prop_inc_vec,Fixed_assets_vec)[0]
ROC6 = np.divide(NOS_vec - alpha6*Prop_inc_vec,Fixed_assets_vec)[0]

# Plotting the ratios
xb =[1959, 1970, 1980, 1990, 2000, 2010, 2018]
ROC1_plt, = plt.plot(ROC1, label ='ROC LS naiv NSA2008')
ROC2_plt, = plt.plot(ROC2, label ='ROC LS-TS NSA2008')
ROC3_plt, = plt.plot(ROC3, label ='ROC LS-TS-NMI NSA2008')
plt.xlabel('Year')
plt.ylabel('Percentage')
plt.title('Question 4 ratio time series\' NSA2008')
plt.xticks(np.arange(0, 60, 10), xb)
plt.legend(handles = [ROC1_plt, ROC2_plt, ROC3_plt], loc='upper right')
plt.show()

xb =[1959, 1970, 1980, 1990, 2000, 2010, 2018]
ROC4_plt, = plt.plot(ROC4, label ='ROC LS naiv NSA1993')
ROC5_plt, = plt.plot(ROC5, label ='ROC LS-TS NSA1993')
ROC6_plt, = plt.plot(ROC6, label ='ROC LS-TS-NMI NSA1993')
plt.xlabel('Year')
plt.ylabel('Percentage')
plt.title('Question 4 ratio time series\' NSA1993 proxy')
plt.xticks(np.arange(0, 60, 10), xb)
plt.legend(handles = [ROC4_plt, ROC5_plt, ROC6_plt], loc='upper right')
plt.show()
