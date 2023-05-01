# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 20:48:39 2023

@author: Simanta Limbu
"""

import numpy as np
import vpython as vp

description = "Triple Pendulum Simulation: The program simulates a three pendulum system.\n" \
              "The first three sliders are for changing mass and the last three are for length\n" \
              "The initial release angle is 90 degrees. The checkboxes for different planet sets"\
              "the gravity at that planet. The trace toggle will turn on trace per user choice."

# Initialization for mass and length
m1 = 1
m2 = 1
m3 = 1
l1 = 1
l2 = 1
l3 = 1

# Initial conditions
theta1 = np.pi/2 #We have released the pendulum form 90 degrees from the axis.
theta2 = np.pi/2
theta3 = np.pi/2
omega1 = 0
omega2 = 0
omega3 = 0
trail_curve = None

# Constant//Initialization at Earth(m/s^2)
g = 9.81

# Planetary gravity values
moon_g = 1.62
earth_g = 9.81
mars_g = 3.71
jupiter_g = 24.79

scene = vp.canvas(width = 1000, height = 500, title = description, caption = 'Mass1(1-20)                 Mass2(1-20)                    Mass3(1-20)                Length1(1-10)                 Length2(1-10)             Length3(1-10)\n')

# Creating sliders for mass and length input
slider_m1 = vp.slider(min= 1, max=20, value=1, length=200, bind=None, text = 'Mass1 (1-20)' )
slider_m2 = vp.slider(min= 1, max=20, value=1, length=200, bind=None, text = 'Mass2(1-20)')
slider_m3 = vp.slider(min= 1, max=20, value=1, length=200, bind=None, text = 'Mass3(1-20')
slider_l1 = vp.slider(min= 1, max=10, value=1, length=200, bind=None, text = 'Length1(1-10)')
slider_l2 = vp.slider(min= 1, max=10, value=1, length=200, bind=None, text = 'Length2(1-10)')
slider_l3 = vp.slider(min= 1, max=10, value=1, length=200, bind=None, text = 'Length3(1-10)')

#Creating labels for the mass
m1_label = vp.label(pos=vp.vec(4, 3, 0), box=False, opacity=0, color = vp.color.red)
m2_label = vp.label(pos=vp.vec(4, 2.5, 0), box=False, opacity=0, color=vp.color.green)
m3_label = vp.label(pos=vp.vec(4, 2, 0), box=False, opacity=0, color=vp.color.blue)
l1_label = vp.label(pos=vp.vec(4, 1.5, 0), box=False, opacity=0,color = vp.color.red )
l2_label = vp.label(pos=vp.vec(4, 1, 0), box=False, opacity=0, color=vp.color.green)
l3_label = vp.label(pos=vp.vec(4, 0.5, 0), box=False, opacity=0, color=vp.color.blue)


# Creating objects
pivot = vp.vector(0, 0, 0)
bob1 = vp.sphere(pos=pivot+vp.vector(l1*np.sin(theta1), -l1*np.cos(theta1), 0), radius=l1/10, color=vp.color.red)
bob2 = vp.sphere(pos=bob1.pos+vp.vector(l2*np.sin(theta2), -l2*np.cos(theta2), 0), radius=l2/10, color=vp.color.green)
bob3 = vp.sphere(pos=bob2.pos+vp.vector(l3*np.sin(theta3), -l3*np.cos(theta3), 0), radius=l3/10, color=vp.color.blue)
rod1 = vp.cylinder(pos=pivot, axis=bob1.pos-pivot, radius=l1/50, color=vp.color.red)
rod2 = vp.cylinder(pos=bob1.pos, axis=bob2.pos-bob1.pos, radius=l2/50, color=vp.color.green)
rod3 = vp.cylinder(pos=bob2.pos, axis=bob3.pos-bob2.pos, radius=l3/50, color=vp.color.blue)

bob3_trail = vp.points(color=vp.color.red, radius=0.05)

# Creating checkboxes for planetary gravity
moon_checkbox = vp.checkbox(text='Moon', bind=lambda checkbox: (update_gravity(moon_g), uncheck_other_checkboxes([moon_checkbox, earth_checkbox, mars_checkbox, jupiter_checkbox], checkbox)))
earth_checkbox = vp.checkbox(text='Earth', bind=lambda checkbox: (update_gravity(earth_g), uncheck_other_checkboxes([moon_checkbox, earth_checkbox, mars_checkbox, jupiter_checkbox], checkbox)))
mars_checkbox = vp.checkbox(text='Mars', bind=lambda checkbox: (update_gravity(mars_g), uncheck_other_checkboxes([moon_checkbox, earth_checkbox, mars_checkbox, jupiter_checkbox], checkbox)))
jupiter_checkbox = vp.checkbox(text='Jupiter', bind=lambda checkbox: (update_gravity(jupiter_g), uncheck_other_checkboxes([moon_checkbox, earth_checkbox, mars_checkbox, jupiter_checkbox], checkbox)))

#Setting Grpahing objects
graph = vp.graph(width=1000, height=500, title="Energy vs Time")
ke1_plot = vp.gcurve(color=vp.color.red, label='KE1')
ke2_plot = vp.gcurve(color=vp.color.green, label='KE2')
ke3_plot = vp.gcurve(color=vp.color.blue, label='KE3')
pe1_plot = vp.gcurve(color=vp.color.orange, label='PE1')
pe2_plot = vp.gcurve(color=vp.color.yellow, label='PE2')
pe3_plot = vp.gcurve(color=vp.color.purple, label='PE3')


graph.legend = True

#Defining checkboxes for gravitational
def uncheck_other_checkboxes(checkbox_list, checked_checkbox):
    for checkbox in checkbox_list:
        if checkbox != checked_checkbox:
            checkbox.checked = False
            
#Defining toggle trace function           
def toggle_trace():
    global trail_curve
    if trace_checkbox.checked:
        if trail_curve is None:
            trail_curve = vp.curve(color=vp.color.red, radius=0.05)
    else:
        if trail_curve is not None:
            trail_curve.clear()
            trail_curve = None

# Creating checkbox for trace
trace_checkbox = vp.checkbox(bind=toggle_trace, text='Trace On/Off', checked=False)


def update_gravity(new_g):
    global g
    g = new_g

def update_labels():
    m1_label.text = f"Mass 1 (kg): {slider_m1.value:.2f}"
    m2_label.text = f"Mass 2 (kg): {slider_m2.value:.2f}"
    m3_label.text = f"Mass 3 (kg): {slider_m3.value:.2f}"
    l1_label.text = f"Length 1 (m): {slider_l1.value:.2f}"
    l2_label.text = f"Length 2 (m): {slider_l2.value:.2f}"
    l3_label.text = f"Length 3 (m): {slider_l3.value:.2f}"

# Simulation loop
dt = 0.01
time_list = []
bob3_trail = None
trail_enabled = False

while True:
    update_gravity(g)
    vp.rate(100)   
     
    # Defining damping constants
    damping1 = 0
    damping2 = 0
    damping3 = 0
    
    #updating the labels
    update_labels()

    
    #Updating mass values accordin to the slider
    m1 = slider_m1.value
    m2 = slider_m2.value
    m3 = slider_m3.value
    l1 = slider_l1.value
    l2 = slider_l2.value
    l3 = slider_l3.value
    
    # Calculating angular accelerations with damping
    alpha1 = (-g*(2*m1 + m2)*np.sin(theta1) - m2*g*np.sin(theta1-2*theta2) - 2*np.sin(theta1-theta2)*m2*(omega2**2*l2 + omega1**2*l1*np.cos(theta1-theta2)) - damping1*omega1) / (l1*(2*m1 + m2 - m2*np.cos(2*theta1-2*theta2)))
    alpha2 = (2*np.sin(theta1-theta2)*((omega1**2)*l1*(m1+m2) + g*(m1+m2)*np.cos(theta1) + (omega2**2)*l2*m2*np.cos(theta1-theta2)) - damping2*omega2) / (l2*(2*m1 + m2 - m2*np.cos(2*theta1-2*theta2)))
    alpha3 = (2*np.sin(theta2-theta3)*((omega2**2)*l2*(m2+m3) + g*(m2+m3)*np.cos(theta2) + (omega3**2)*l3*m3*np.cos(theta2-theta3)) - damping3*omega3) / (l3*(2*m2 + m3 - m3*np.cos(2*theta2-2*theta3)))
    
    # Updating angular velocities with damping
    omega1 += alpha1*dt - damping1*omega1*dt
    omega2 += alpha2*dt - damping2*omega2*dt
    omega3 += alpha3*dt - damping3*omega3*dt
    
    # Updating angular positions
    theta1 += omega1*dt
    theta2 += omega2*dt
    theta3 += omega3*dt
    
    # Updating bob and rod positions
    bob1.pos = pivot + vp.vector(l1*np.sin(theta1), -l1*np.cos(theta1), 0)
    bob2.pos = bob1.pos + vp.vector(l2*np.sin(theta2), -l2*np.cos(theta2), 0)
    bob3.pos = bob2.pos + vp.vector(l3*np.sin(theta3), -l3*np.cos(theta3), 0)
    rod1.axis = bob1.pos - pivot
    rod2.pos = bob1.pos
    rod2.axis = bob2.pos - bob1.pos
    rod3.pos = bob2.pos
    rod3.axis = bob3.pos - bob2.pos
    
    #Update trail if user wants to see it
    if trail_curve is not None:
        trail_curve.append(bob3.pos)
    
    # Calculating kinetic energy with damping
    KE1 = 0.5 * m1 * l1**2 * omega1**2 * (1 - damping1*dt)
    KE2 = 0.5 * m2 * (l1**2 * omega1**2 + l2**2 * omega2**2 + 2*l1*l2*omega1*omega2*np.cos(theta1-theta2)) * (1 - damping2*dt)
    KE3 = 0.5 * m3 * (l2**2 * omega2**2 + l3**2 * omega3**2 + 2*l2*l3*omega2*omega3*np.cos(theta2-theta3)) * (1 - damping3*dt)
    
    # Calculating potential energy with damping
    PE1 = m1 * g * l1 * (1 - np.cos(theta1))
    PE2 = m2 * g * (l1*(1-np.cos(theta1)) + l2*(1-np.cos(theta2)))
    PE3 = m3 * g * (l2*(1-np.cos(theta2)) + l3*(1-np.cos(theta3)))
    
    t = dt * len(time_list)
    time_list.append(t)

  
    vp.rate(100)  # Set the frame rate
    ke1_plot.plot(t, KE1)
    ke2_plot.plot(t, KE2)
    ke3_plot.plot(t, KE3)
    pe1_plot.plot(t, PE1)
    pe2_plot.plot(t, PE2)
    pe3_plot.plot(t, PE3)