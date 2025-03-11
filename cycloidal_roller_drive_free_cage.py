"""
Cycloidal Drive with Intermediate Rollers and Free Cage animation script + export to dxf.

YouTube: https://youtube.com/MishinMachine
Instagram: https://instagram.com/mishinmachine

Shoutout to @Tomato1107 for his original cycloidal-drive animation scripts.
Source: https://github.com/Tomato1107/Cycloidal-Drive-Animation
"""

import ezdxf as ezdxf
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button
import numpy as np

OUT_FILE = __file__.rsplit("/", 1)[0] + "/wdirfc_output.dxf"

# ------------------------------------------------------------------------
# Initial values
# ------------------------------------------------------------------------
INITIAL_D = 45.8               # separator diameter
INITIAL_d = 3                # pin diameter
INITIAL_e = INITIAL_d * 0.2  # eccentricity of the pins
INITIAL_n = 29               # number of pins
INITIAL_isd = 5.0            # input shaft diameter
INITIAL_ced = 15.0           # cam eccentric diameter

RESOLUTION = 50             # Resolution per pin

# visual and  animation
DOT_SIZE = 4
INTERVAL = 50 # ms

# initial rotation angle
rotation_angle = 0

# ------------------------------------------------------------------------
# Plot setup
# ------------------------------------------------------------------------
fig, ax = plt.subplots()
ax.grid()

# Set initial limits
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_aspect('equal', 'box')

# Reserve space at the bottom for constols
plt.subplots_adjust(bottom=0.48)

base_dot, = ax.plot([], [], 'o', ms=DOT_SIZE, color='gray')
base_circle = plt.Circle((INITIAL_e, 0),
                          radius=INITIAL_D/2,
                          fill=False,
                          linestyle='--',
                          edgecolor='gray',
                          linewidth=0.6)
ax.add_patch(base_circle)

cam_dot,       = ax.plot([], [], 'ro', ms=DOT_SIZE)
cam_gear, = ax.plot([], [], 'r-', 
                    label='Cam Gear',
                    )

crown_gear, = ax.plot([], [], 'b-',
                      label='Crown Gear',
                      )

ishaft_dot, = ax.plot([], [], 'ko', ms=DOT_SIZE)
ishaft_circle = plt.Circle((0, 0),
                                radius=INITIAL_isd/2,
                                fill=False,
                                edgecolor='black',
                                linewidth=1.5)
ax.add_patch(ishaft_circle)


ecc_circle = plt.Circle((2*INITIAL_e, 0),
                        radius=INITIAL_ced/2,
                        fill=False,
                        edgecolor='red',
                        linewidth=1.5)
ax.add_patch(ecc_circle)

pins = []
pin_dot, = ax.plot([], [], 'go', ms=5)

# ------------------------------------------------------------------------
# Draw functions
# ------------------------------------------------------------------------
def draw_pins(n, e, D, d, rotate = 0, offset = [0, 0]):
    pins = []
    R = D / 2 # base radius
    r = d / 2 # pin radius

    pinCos = np.cos(rotate / n)
    pinSin = np.sin(rotate / n)

    angles = np.linspace(0, 2*np.pi, n, endpoint=False)
    for angle in angles:
        xT = R * np.cos(angle) # + e
        yT = R * np.sin(angle)

        # Rotate the pin
        x = xT * pinCos - yT * pinSin + offset[0]
        y = xT * pinSin + yT * pinCos + offset[1]

        c = plt.Circle((x, y),
                       radius=r,
                       fill=False,
                       edgecolor='green',
                       linewidth=1.5)
        pins.append(c)
    return pins

# Draw Hypocycloid or Epicycloid
def draw_cycloid(D, d, e, n, epicycloid = False, offset = [0,0], angle_rad = 0):
    sign = epicycloid and -1 or 1 # sign for outer or inner cycloid

    R = D / 2 # base radius
    r = d / 2 # pin radius
    
    Rn = (R / n)
    Rz = (n - sign * 1) * Rn
    Rzn = Rz + sign * Rn

    cosA = np.cos(angle_rad)
    sinA = np.sin(angle_rad)

    t = np.linspace(0, 2 * np.pi, RESOLUTION * n)
    
    sin_t = np.sin(t)
    cos_t = np.cos(t)
    
    factor = Rzn / Rn 
    sin_factor_t = np.sin(factor * t)
    cos_factor_t = np.cos(factor * t)
    
    # Compute the base cycloid coordinates
    xa = Rzn * cos_t - sign * e * cos_factor_t
    ya = Rzn * sin_t - e * sin_factor_t

    dxa = Rzn * (-sin_t + sign * (e / Rn) * sin_factor_t)
    dya = Rzn * (cos_t - (e / Rn) * cos_factor_t)

    denom = np.sqrt(dxa**2 + dya**2)
    # Offset the cycloid coordinates along the normal direction
    xT = xa + sign * (r / denom) * (-dya) 
    yT = ya + sign * (r / denom) * dxa
    
    # Rotate by angle_rad and then apply the offset
    x = xT * cosA - yT * sinA + offset[0]
    y = xT * sinA + yT * cosA + offset[1]

    return [x, y]

def draw_cam(D, d, e, n, angle = 0):
    # angle = angle + np.pi / (n - 1) # if cam is not eccentric
    off = 2 * e
    offX = off * np.cos(-angle) 
    offY = off * np.sin(-angle)
    
    return draw_cycloid(D, d, e, n, False, [offX, offY], angle / ((n-1)/2))

def draw_crown(D, d, e, n, angle = 0):
    return draw_cycloid(D, d, e, n, True, [0, 0], angle)

# ------------------------------------------------------------------------
# Updates
# ------------------------------------------------------------------------
def update_slider_limits ():
    D_val = D_slider.val
    d_val = d_slider.val
    e_val = e_slider.val


    min_D = d_val * 2
    D_slider.valmin = min_D
    D_slider.ax.set_xlim(D_slider.valmin, D_slider.valmax)
    if D_val < min_D:
        D_slider.set_val(min_D)
        D_val = min_D

    max_d = D_val / 2
    d_slider.valmax = max_d
    d_slider.ax.set_xlim(d_slider.valmin, d_slider.valmax)
    if d_slider.val > max_d:
        d_slider.set_val(max_d)
        d_val = max_d

    max_e = round((min(d_val, max_d) / 4), 1)
    e_slider.valmax = max_e
    e_slider.ax.set_xlim(e_slider.valmin, e_slider.valmax)
    if e_slider.val > max_e:
        e_slider.set_val(max_e)
        e_val = max_e

    n_max = np.ceil((np.pi * D_val) / (d_val))
    n_slider.valmax = n_max
    n_slider.ax.set_xlim(n_slider.valmin, n_slider.valmax)
    if n_slider.val > n_max:
        n_slider.set_val(n_max)
      
    isd_max = D_val - d_val - e_val*5 - 2
    isd_slider.valmax = isd_max
    isd_slider.ax.set_xlim(isd_slider.valmin, isd_slider.valmax)

    if isd_slider.val > isd_max:
        isd_slider.set_val(isd_max)

    ced_min = isd_slider.val + e_val * 5 
    ced_max = D_val - d_val - e_val * 2 - 1

    if ced_min > ced_max:
        ced_min = ced_max

    ced_slider.valmin = ced_min
    ced_slider.valmax = ced_max
    ced_slider.ax.set_xlim(ced_min, ced_slider.valmax)

    if ced_slider.val < ced_min:
        ced_slider.set_val(ced_min)

    if ced_slider.val > ced_max:
        ced_slider.set_val(ced_max)

def update(_):
    update_slider_limits()

    D_val = D_slider.val
    d_val = d_slider.val
    e_val = e_slider.val
    n_val = int(n_slider.val)
    isd_val = isd_slider.val
    ced_val = ced_slider.val

    angle = rotation_angle / 180 * np.pi 

    cosA = np.cos(angle)
    sinA = np.sin(angle)

    eccOffX = e_val * cosA
    eccOffY = e_val * -sinA

    # Separator circle
    base_r = D_val / 2
    base_circle.center = (eccOffX, eccOffY)
    base_circle.set_radius(base_r)
    base_circle.set_linestyle('--')  # dashed

    # Cam gear
    cam_points = draw_cam(D_val, d_val, e_val, n_val, angle )
    cam_gear.set_data(cam_points[0], cam_points[1])
    cam_dot.set_data([cam_points[0][0]], [cam_points[1][0]])

    # Crown gear
    crown_points = draw_crown(D_val, d_val, e_val, n_val)
    crown_gear.set_data(crown_points[0], crown_points[1])

    # Input shaft
    isr = isd_val / 2
    ishaft_circle.set_radius(isr)
    ishaft_dot.set_data([isr * cosA], [isr * -sinA])

    # Cam eccentric circle
    ecc_circle.center = (eccOffX * 2, eccOffY * 2)
    ecc_circle.set_radius(ced_val / 2)

    # Remove old pins from the plot
    for p in pins:
        p.remove()
    pins.clear()
    
    # Create new pins
    new_pins = draw_pins(n_val, e_val, D_val, d_val, angle, [eccOffX, eccOffY])
    pins.extend(new_pins)
    for p in pins:
        ax.add_patch(p)

    pinX = pins[0].center[0]
    pinY = pins[0].center[1]

    pin_r = d_val / 2
    base_dot.set_data([pinX], [pinY])
    pin_dot.set_data([pinX + pin_r * cosA ], [pinY + pin_r * -sinA])

    # Adjust plot limits so everything is visible
    body_D = D_val + d_val + e_val * 2

    lim = body_D/2 + 1
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    
    ax.set_aspect('equal', 'box')
    fig.canvas.draw_idle()

# ------------------------------------------------------------------------
# Sliders
# ------------------------------------------------------------------------
slider_width = 0.65
slider_height = 0.03

D_ax = plt.axes([0.2, 0.30, slider_width, slider_height])
d_ax = plt.axes([0.2, 0.25, slider_width, slider_height])
e_ax = plt.axes([0.2, 0.20, slider_width, slider_height])
n_ax = plt.axes([0.2, 0.15, slider_width, slider_height])
isd_ax = plt.axes([0.2, 0.10, slider_width, slider_height])
ced_ax = plt.axes([0.2, 0.05, slider_width, slider_height])

# D_ax = plt.axes([0.2, 0.0, slider_width, slider_height])
# d_ax = plt.axes([0.2, 0.0, slider_width, slider_height])
# e_ax = plt.axes([0.2, 0.0, slider_width, slider_height])
# n_ax = plt.axes([0.2, 0.0, slider_width, slider_height])
# isd_ax = plt.axes([0.2, 0.0, slider_width, slider_height])
# ced_ax = plt.axes([0.2, 0.0, slider_width, slider_height])

D_slider = Slider(D_ax, 'Pins base: D', 10, 200.0, valinit=INITIAL_D, valstep=0.5)
d_slider = Slider(d_ax, 'Pin dia: d',       1, 20,  valinit=INITIAL_d, valstep=0.5)
e_slider = Slider(e_ax, 'Eccentricity: e',      0.1, 5.0,  valinit=INITIAL_e, valstep=0.1)
n_slider = Slider(n_ax, 'Pins num: n',     3,   50,   valinit=INITIAL_n, valstep=1)
isd_slider = Slider(isd_ax, 'Input shaft dia: isd',   0.1, 20.0, valinit=INITIAL_isd, valstep=0.1)
ced_slider = Slider(ced_ax, 'Cam ecc dia: ced',   0.1, 20.0, valinit=INITIAL_ced, valstep=0.1)

D_slider.on_changed(update)
d_slider.on_changed(update)
e_slider.on_changed(update)
n_slider.on_changed(update)
isd_slider.on_changed(update)
ced_slider.on_changed(update)


# ------------------------------------------------------------------------
# Controls & Export
# ------------------------------------------------------------------------
axcolor = 'lightgoldenrodyellow'
reset_btn = Button(plt.axes([0.2, 0.37, 0.1, 0.05]), 'Reset', color="#EEEEEE", hovercolor='0.975')
export_btn = Button(plt.axes([0.71, 0.37, 0.14, 0.05]), 'Export dxf', color='#00FF00', hovercolor='0.975')
status_text = plt.text(0.5, 0.39, '', ha='center', va='center', fontsize=10, transform=plt.gcf().transFigure)

def reset(_): 
    global rotation_angle
    rotation_angle = 0

    D_slider.reset()    
    d_slider.reset()
    e_slider.reset()
    n_slider.reset()
    isd_slider.reset()
    ced_slider.reset()    

    status_text.set_text("")
    update(_)

reset_btn.on_clicked(reset)

def export(_):
    D_val = D_slider.val
    d_val = d_slider.val
    e_val = e_slider.val
    n_val = int(n_slider.val)
    isd_val = isd_slider.val
    ced_val = ced_slider.val

    status_text.set_text("Saving... ")

    doc = ezdxf.new("R2000")
    msp = doc.modelspace()

    # Separator circle
    msp.add_circle((e_val, 0), radius=D_val / 2)

    # Cam gear
    cam_points = draw_cam(D_val, d_val, e_val, n_val)
    msp.add_lwpolyline(np.stack(cam_points, axis=1))

    # Crown gear
    crown_points = draw_crown(D_val, d_val, e_val, n_val)
    msp.add_lwpolyline(np.stack(crown_points, axis=1))

    # Pins
    R = D_val / 2
    r = d_val / 2
    angles = np.linspace(0, 2*np.pi, n_val, endpoint=False)
    for angle in angles:
        x = R * np.cos(angle) + e_val
        y = R * np.sin(angle)
        c = plt.Circle((x, y),
                       radius=r,
                       fill=False,
                       edgecolor='green',
                       linewidth=1.5)
        msp.add_circle((x, y), radius=r)

    # Input shaft
    msp.add_circle((0, 0), radius=isd_val / 2)

    # Eccentric circle
    msp.add_circle((2 * e_val, 0), radius=ced_val / 2)

    doc.saveas(OUT_FILE)

    status_text.set_text("Saved to " + OUT_FILE)
    fig.canvas.draw_idle()

export_btn.on_clicked(export)


# ------------------------------------------------------------------------
# Initialize
# ------------------------------------------------------------------------
step = 20
def animate(_):
    global rotation_angle
    n = int(n_slider.val)
    rotation_angle = (rotation_angle + step) % (360 * n)
    update(_)

ani = animation.FuncAnimation(fig, animate, frames=int(360/step), interval=INTERVAL, repeat=True)

plt.show()
