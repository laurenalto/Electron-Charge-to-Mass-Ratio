import numpy as np
import math
import matplotlib.pyplot as plt
import plotly.express as px
import statsmodels.api as sm
from scipy.stats import chisquare

# #first plot - constant voltage (i.e. current)
# def b_value(array): # we don't need this
#     uo  = 4*np.pi *10**(-7) 
#     R = 0.15
#     b = [0,0,0,0,0,0,0,0,0,0,0,0]
#     # WHAT IS Z? assuming R/2 because bulb in center of two coils
#     # Q: do we need to multiply field value by 2 for both coils? 
#     # 
#     for i in range(len(array)):
#         #print(array[i], radius [i])
#         z = R/2
#         b[i] = (uo *array[i] *(R)**2)/(2*((R)**2 + z**2))**1.5
#     return b

# def b_uncertainity(current):
#     R = 0.15
#     z = R/2
#     uo  = 4*np.pi *10**(-7) 
#     d_I = 0.0001
#     d_R = 0.0005
#     d_z = d_R
#     d_B = [0,0,0,0,0,0,0,0,0,0,0,0]
#     for i in range(len(current)):
#         d_B[i] = math.sqrt((uo*R**2/(2*(R**2 + z**2)**1.5)*d_I)**2 + (uo*current[i]*(-R**2 + 2*z**2)/(2*(R**2 + z**2)**2.5)*d_R)**2 + (-3*uo*current[i]*R**2/(2*(R**2 + z**2)**2.5)*d_z)**2)
#     return d_B

# def mass(radius, b): #  WILL NEED TO CONSIDER CHANGING VOLTAGE FOR OTHER CALCULATIONS - > might not need mass value if done with ratio?
#     voltage = 213.9
#     e = 1.60217663 * 10**(-19)
#     m = [0,0,0,0,0,0,0,0,0,0,0,0]
    
#     for i in range(len(radius)):
#         m[i] = (b[i]**2 * (radius[i])**2 * e)/(2 * voltage)
#     return m

# def m_uncertainty(radius,b, d_b): # WILL NEED TO CONSIDER CHANGING VOLTAGE FOR OTHER CALCULATIONS
#     voltage = 213.9
#     d_v = 0.01
#     e = 1.60217663 * 10**(-19)
#     d_m = [0,0,0,0,0,0,0,0,0,0,0,0]
#     d_r = 0.0005
#     for i in range(len(radius)):
#         d_m[i] = math.sqrt((b[i]*radius[i]**2*e/voltage * d_b[i])**2 + (b[i]**2*radius[i]*e/voltage * d_r)**2 + (-b[i]**2*radius[i]**2*e/(2*voltage**2) * d_v)**2)
#     return d_m

def Bc_value(current):
    uo  = 4*np.pi *10**(-7) 
    R = 0.16
    n  = 130
    Bc = [0,0,0,0,0,0,0,0,0,0,0,0]
    a = 0.8**1.5
    c = uo *n  /R
    for i in range(len(current)):
        Bc[i] = a* c* current[i]
    return Bc

def Bc_uncertainty(current):
    uo  = 4*np.pi *10**(-7) 
    R = 0.16
    d_R = 0.0005
    d_I = 0.0001
    n  = 130
    d_Bc = [0,0,0,0,0,0,0,0,0,0,0,0]
    for i in range(len(current)):
        d_Bc[i] = math.sqrt( (0.8**1.5 * uo * n / R *d_I )**2 + (-0.8**1.5 *uo*n*current[i]/R**2 *d_R)**2)
    return d_Bc
    
def correction_factor(radius, b):
    R = 0.16
    b_corrected = [0,0,0,0,0,0,0,0,0,0,0,0]
    correction_factor = [0,0,0,0,0,0,0,0,0]

    for i in range(len(correction_factor)):
       x = (radius[i])**4/(R**4*(0.6583 + 0.29*((radius[i])**2/R**2))**2)
       # print(x)
       correction_factor[i] = 1-x
       b_corrected[i] = correction_factor[i] * b[i]
       
    for i in range(9,12):
        b_corrected[i] = b[i]*1.00075
    return correction_factor, b_corrected

def uncertainty_inv_radius(radius, uncertainty):
    uncertaintyr =  [0,0,0,0,0,0,0,0,0,0,0,0]
    for i in range(len(uncertainty)):
        uncertaintyr[i] = uncertainty[i]/(radius[i]**2)
    return uncertaintyr

# CHECK HERE
def ratio(voltage, radius_volt, current, radius_cur):
    n = 130
    R = 0.16
    uo  = 4*np.pi *10**(-7) 
    k = 1/math.sqrt(2) * 0.8**1.5 * uo*n/R
    ratio_cur = [0,0,0,0,0,0,0,0,0,0,0,0]
    ratio_volt= [0,0,0,0,0,0,0,0,0,0]
    #io = 0.000036/k
    io = 0.08710826686488
    volt = 124.7
    cur = 1.87
    for i in range(len(current)):
        a = math.sqrt(volt)/ radius_cur[i]
        b = (current[i] + io/math.sqrt(2))
        ratio_cur[i] = a/b * 1/k
        ratio_cur[i] = ratio_cur[i]**2
        # ratio_cur[i] = volt/(radius_cur[i]**2 * k**2 * (current[i] + io/math.sqrt(2))**2)
        print(k,io)
    
    for i in range(len(voltage)):
        ratio_volt[i] = voltage[i]/(radius_volt[i]**2 * k**2 * (cur + io/math.sqrt(2))**2)
    return ratio_cur, ratio_volt

def un_k():
    d_R = 0.0005
    n = 130
    R = 0.16
    uo  = 4*np.pi *10**(-7) 
    d_k = 1/math.sqrt(2) * 0.8**1.5 * uo*n/R**2 * d_R
    return d_k


def un_io(d_k):
    d_Be = 0.000004
    Be = 0.000036
    k = 0.0005165985898213973
    a = (d_Be/k)**2
    b = (Be/k**2 * d_k)**2
    d_io = math.sqrt(a + b)
    return d_io

def un_em(d_io, d_k, radius):
    d_r = 0.0005
    d_I = 0.0001
    I = 1.87
    V = 213.9
    d_V = 0.01
    k = 0.0005165985898213973
    io = 0.06968660137544358
    a = (io/math.sqrt(2) + I)**2
    b = (d_V/(radius[4]**2 * k**2 *a))**2
    c = (d_r * 2*V/(radius[4]**3 * k**2 * a))**2
    f = (d_k * 2*V/(radius[4]**2 * k**3 * a))**2
    d = (d_io * 4*V/(radius[4]**2 * k**2 * a**1.5))**2
    e = (d_I * 2*V/(radius[4]**2 * k**2 * a))**2
    d_em = math.sqrt(b + c + d + e + f)
    return d_em

# data ( radius(y) in cm)
x = np.array([0.899, 0.917, 0.961, 1.033, 1.073, 1.136, 1.204, 1.265, 1.388, 1.536, 1.625, 1.870]) #A
y = np.array([0.048, 0.0455, 0.044, 0.041, 0.039, 0.037, 0.0355, 0.0335, 0.0305, 0.027, 0.0265, 0.023]) # radius [m]
dx = np.array([0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001])
dy = np.array([0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005])

temp_er = [0,0,0,0,0,0,0,0,0,0,0,0] # change this with appropriate errors


#taking inverse of radius
rad_inv = [0,0,0,0,0,0,0,0,0,0,0,0] 

for i in range(len(y)):
    rad_inv[i] = 1/y[i]

# b = b_value(x)
# print(b, "^b\n")


# d_b = b_uncertainity(x)
# print(d_b, "^d_b\n")
Bc = Bc_value(x)
print(Bc, "^Bc\n")
d_Bc = Bc_uncertainty(x)
print(d_Bc, "^d_Bc\n")

correction_factorvalues, bc_corrected = correction_factor(y, Bc)
print(bc_corrected, "^bc corrected\n")
d_Bc_corrected = Bc_uncertainty(x)
print(d_Bc, "^d_Bc\n")
print(correction_factorvalues, "^correction factor\n")
print(rad_inv, "^inverse radius\n")

rad_un_inv = uncertainty_inv_radius(y, dy)
print(rad_un_inv, "^ inverse radius uncertainity\n")

x2 = np.array([124.70, 150.20, 175.40, 201.50, 225.50, 250.50, 275.80, 300.20, 325.60, 342.70]) #Voltage
y2 = np.array([0.0165, 0.019, 0.0205, 0.0215, 0.023, 0.0245, 0.026, 0.0275, 0.0285, 0.03])
# 
test_current = np.array([1.485, 1.349, 1.22, 1.125, 1.058, 0.956, 0.896, 0.85, 0.8, 0.765])
test_radius = np.array([0.0285, 0.0315, 0.033, 0.036, 0.039, 0.0425, 0.0445, 0.0475, 0.05, 0.0525])
#cur, volt = ratio(x2,y2, test_current, test_radius)
cur, volt = ratio(x2,y2, x, y)
print(cur, " ^ e/m current\n", volt, "^e/m voltage\n")

d_k = un_k()
d_io = un_io(d_k)
d_em = un_em(d_io, d_k, y)
print(d_k, d_io, d_em)



# print(y) # printing radius
# m = mass(y, b)
# print(m, "^mass\n")
# d_m = m_uncertainty(y,b,d_b)
# print(d_m, "^d_mass\n")

# bc_corrected = np.array([0.0006454343119309968, 0.000660516035004621, 0.0006934102918462907, 0.0007475972666039858, 0.0007778491553566384, 0.0008247121287476549, 0.0008749067692954954, 0.0009202400876536487, 0.0010110528594557144, 0.0011230136333986666, 0.0011880840848130423, 0.0013672106083694702])
# rad_inv = np.array([20.833333333333332, 21.978021978021978, 22.72727272727273, 24.390243902439025, 25.641025641025642, 27.027027027027028, 28.169014084507044, 29.850746268656714, 32.78688524590164, 37.03703703703704, 37.735849056603776, 43.47826086956522])
# # PLOTTING FOR 3
# # Assuming you have already defined rad_inv, b_corrected, dy, temp_er
# fig = px.scatter(x=rad_inv, y=bc_corrected, error_x=rad_un_inv, error_y=d_Bc, trendline="ols", trendline_color_override='darkblue')
# fig.update_traces(marker=dict(size=12, line=dict(width=2, color='DarkSlateGrey')), selector=dict(mode='markers'))
# fig.update_layout(
#     title=dict(text="Coil Magnetic Field (Bc) as a Function of Electron Path Radius Inverse (1/r)", font=dict(size=30), x=0.5),
#     xaxis_title=dict(text="Inverse Radius [1/m]", font=dict(size=30)),
#     yaxis_title=dict(text="Magnetic Field from Coils Bc [T]", font=dict(size=30))
# )
# fig.update_xaxes(tickfont=dict(size=30))
# fig.update_yaxes(tickfont=dict(size=30))
# fig.show()
# results = px.get_trendline_results(fig)
# results = results.iloc[0]["px_fit_results"].summary()
# # print(results)

# x_w_c = sm.add_constant(rad_inv)
# model = sm.OLS(bc_corrected,x_w_c)
# results = model.fit()
# print("Intercept (b):", results.params[0])
# print("Slope (m_2):", results.params[1])
# #find error in slope and intercept
# print(results.summary())

# #Bc RESIDUALS
# residuals = bc_corrected - (results.params[1]*rad_inv + results.params[0])
# chisq = chisquare(bc_corrected, f_exp=(results.params[1]*rad_inv + results.params[0]))
# print(chisq)
# print(chisq)
# fig = px.scatter(x=rad_inv, y=residuals, error_x=rad_un_inv, error_y=d_Bc, trendline="ols", trendline_color_override='darkblue')
# fig.update_layout(title=dict(text="Residuals for Coils Magnetic Field Plot", font=dict(size=30), x=0.5), xaxis_title=dict(text="Inverse Radius [1/m]", font =dict(size=30)), yaxis_title=dict(text="Magnetic Field from Coils Bc [T]", font=dict(size=30)))
# fig.update_traces(marker=dict(size=12,
#                               line=dict(width=2,
#                                         color='DarkSlateGrey')),
#                   selector=dict(mode='markers'))
# fig.update_xaxes(tickfont=dict(size=30))
# fig.update_yaxes(tickfont=dict(size=30))
# fig.show()
# results
# results = px.get_trendline_results(fig)
# results = results.iloc[0]["px_fit_results"].summary()
# print(results)

# for i in range(len(bc_corrected)):
#     print(bc_corrected[i])
#     print(bc_corrected[i] + -0.000036, "\n")
    
# ## FIRST PLOT - Changing Current
# fig = px.scatter(x=x, y=y, error_x=dx, error_y=dy ,trendline="ols", trendline_color_override='darkblue')
# fig.update_traces(marker=dict(size=12,
#                               line=dict(width=2,
#                                         color='DarkSlateGrey')),
#                   selector=dict(mode='markers'))
# fig.update_layout(title=dict(text="Coil Current in Relation to the Electron Path Radius Given a Constant Voltage",font=dict(size=30), x=0.5), xaxis_title=dict(text="Current [A]", font=dict(size=30)), yaxis_title=dict(text="Electron Path Radius [cm]", font=dict(size=30)))
# fig.update_xaxes(tickfont=dict(size=30))
# fig.update_yaxes(tickfont=dict(size=30))
# fig.show()

# #chi squared test
# #DO results.bse

# # #display trendline results on plot
# # results = px.get_trendline_results(fig)
# # print(results)

# x_w_c = sm.add_constant(x)
# model = sm.OLS(y,x_w_c)
# results = model.fit()
# print("Intercept (b):", results.params[0])
# print("Slope (m_2):", results.params[1])
# #find error in slope and intercept
# print(results.summary())

# #CURRENT RESIDUALS
# residuals = y - (results.params[1]*x + results.params[0])
# chisq = chisquare(y, f_exp=(results.params[1]*x + results.params[0]))
# print(chisq)
# print("hi")
# print(chisq)
# fig = px.scatter(x=x, y=residuals, error_x=dx, error_y=dy, trendline="ols", trendline_color_override='darkblue')
# fig.update_layout(title=dict(text="Residuals Plot for Coil Current in Relation to Electron Path Radius", font=dict(size=30), x=0.5), xaxis_title=dict(text="Current [A]", font =dict(size=30)), yaxis_title=dict(text="Electron Path Radius [cm]", font=dict(size=30)))
# fig.update_traces(marker=dict(size=12,
#                               line=dict(width=2,
#                                         color='DarkSlateGrey')),
#                   selector=dict(mode='markers'))
# fig.update_xaxes(tickfont=dict(size=30))
# fig.update_yaxes(tickfont=dict(size=30))
# fig.show()


# #SECOND PLOT - Changing Voltage
# x2 = np.array([124.70, 150.20, 175.40, 201.50, 225.50, 250.50, 275.80, 300.20, 325.60, 342.70]) #Voltage
# y2 = np.array([1.65, 1.9, 2.05, 2.15, 2.3, 2.45, 2.6, 2.75, 2.85, 3])
# dx2 = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01,0.01, 0.01, 0.01, 0.01])
# dy2 = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])

# fig = px.scatter(x=x2, y=y2, error_x=dx2, error_y=dy2, trendline="ols", trendline_color_override='darkblue')
# fig.update_layout(title = dict(text="Coil Voltage in Relation to the Electron Path Radius Given a Constant Current", font=dict(size=30), x=0.5), xaxis_title=dict(text="Voltage [V]", font=dict(size=30)), yaxis_title=dict(text="Electron Path Radius [cm]", font=dict(size=30)))
# fig.update_traces(marker=dict(size=12,
#                               line=dict(width=2,
#                                         color='DarkSlateGrey')),
#                   selector=dict(mode='markers'))
# fig.update_xaxes(tickfont=dict(size=30))
# fig.update_yaxes(tickfont=dict(size=30))
# fig.show()

# #display trendline results on plot
# #results = px.get_trendline_results(fig)
# #print(results)

# # VOLTAGE RESIDUALS
# x2_w_c = sm.add_constant(x2)
# model = sm.OLS(y2,x2_w_c)
# results = model.fit()
# print("Intercept (b):", results.params[0])
# print("Slope (m_2):", results.params[1])
# print(results.summary())
# chisq = chisquare(y2, f_exp=(results.params[1]*x2 + results.params[0]))
# print(chisq)
# residuals2 = y2 - (results.params[1]*x2 + results.params[0])
# fig = px.scatter(x=x2, y=residuals2, error_x=dx2, error_y=dy2, trendline="ols", trendline_color_override='darkblue')
# fig.update_layout(title = dict(text="Residuals Plot for Coil Voltage in Relation to Electron Path Radius", font=dict(size=30), x=0.5), xaxis_title=dict(text="Voltage [V]", font =dict(size=30)), yaxis_title=dict(text="Electron Path Radius [cm]", font=dict(size=30)))
# fig.update_traces(marker=dict(size=12,
#                               line=dict(width=2,
#                                         color='DarkSlateGrey')),
#                   selector=dict(mode='markers'))
# fig.update_xaxes(tickfont=dict(size=30))
# fig.update_yaxes(tickfont=dict(size=30))
# fig.show()