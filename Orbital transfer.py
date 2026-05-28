import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

# CONSTANTS
mu = 398600       # Earth's gravitational parameter (km^3/s^2)
R_earth = 6371    # Earth's radius (km)

# FUNCTION: ORBIT FROM e, p
def generate_orbit(e, p, theta_start=0, theta_end=360):
    theta_deg = np.arange(theta_start, theta_end, 0.1)
    theta = np.radians(theta_deg)

    r = p / (1 + e * np.cos(theta))

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    return theta, r, x, y


# FUNCTION: VELOCITY COMPONENTS
def velocity(e, p, theta, r):
    h = np.sqrt(mu * p)

    v_r = (mu / h) * e * np.sin(theta)
    v_theta = h / r

    return v_r, v_theta


# FUNCTION: TRANSFORM TO 3D (ECI)
def transform_to_3d(x_pf, y_pf, inc_deg, omega_deg, raan_deg=0.0):
    """Transform perifocal (orbital plane) coordinates to 3D ECI-like frame.

    Rotation chain:  Rz(-RAAN) . Rx(-i) . Rz(-omega)

    Parameters:
        x_pf, y_pf : perifocal coordinates (periapsis along x)
        inc_deg    : orbital inclination (degrees)
        omega_deg  : argument of perigee (degrees)
        raan_deg   : RAAN (degrees), default=0
    Returns:
        x_3d, y_3d, z_3d  in the reference frame
    """
    i = np.radians(inc_deg)
    w = np.radians(omega_deg)
    O = np.radians(raan_deg)

    # Rotate by omega in orbital plane
    x_rot = x_pf * np.cos(w) - y_pf * np.sin(w)
    y_rot = x_pf * np.sin(w) + y_pf * np.cos(w)

    # Tilt by inclination (around x-axis)
    x_tilt = x_rot
    y_tilt = y_rot * np.cos(i)
    z_tilt = y_rot * np.sin(i)

    # Rotate by RAAN (around z-axis)
    x_3d = x_tilt * np.cos(O) - y_tilt * np.sin(O)
    y_3d = x_tilt * np.sin(O) + y_tilt * np.cos(O)
    z_3d = z_tilt

    return x_3d, y_3d, z_3d


# FUNCTION: TRUE -> MEAN ANOMALY
def true_to_mean_anomaly(theta, e):
    """Convert true anomaly to mean anomaly (both in radians).
    Uses eccentric anomaly as intermediate.  Valid for e < 1 only.
    """
    if e >= 1:
        return theta
    E = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(theta / 2),
                        np.sqrt(1 + e) * np.cos(theta / 2))
    M = E - e * np.sin(E)
    return M


# FUNCTION: CLASSIFY ORBIT
def classify_orbit(e, a, r_peri):
    """Classify orbit type based on eccentricity and geometry."""
    if e == 0:
        return "Circular"
    elif 0 < e < 1:
        if r_peri < R_earth:
            return "Elliptical (sub-surface periapsis -- not physical!)"
        elif a < 35786 + R_earth:
            return "Elliptical (LEO/MEO region)"
        else:
            return "Elliptical (HEO / near-GEO region)"
    elif e == 1:
        return "Parabolic (escape trajectory)"
    else:
        return "Hyperbolic (escape trajectory)"


def orbit_type_simple(e):
    if e == 0:
        return "Circular"
    elif 0 < e < 1:
        return "Elliptical"
    elif e == 1:
        return "Parabolic"
    else:
        return "Hyperbolic"


# FUNCTION: FULL ORBITAL ANALYSIS
def orbital_analysis(name, e, p, inc=0.0, omega=0.0, M0=0.0):
    """Compute and return all Keplerian elements and derived quantities."""
    a = p / (1 - e**2)                       # semi-major axis
    b = a * np.sqrt(1 - e**2)                # semi-minor axis
    r_peri = a * (1 - e)                     # periapsis radius
    r_apo = a * (1 + e)                      # apoapsis radius
    h = np.sqrt(mu * p)                      # specific angular momentum
    energy = -mu / (2 * a)                   # specific orbital energy
    T = 2 * np.pi * np.sqrt(a**3 / mu)       # orbital period (seconds)

    # Velocities at key points
    v_peri = np.sqrt(mu * (2/r_peri - 1/a))  # velocity at periapsis
    v_apo = np.sqrt(mu * (2/r_apo - 1/a))    # velocity at apoapsis

    # Altitudes above Earth surface
    alt_peri = r_peri - R_earth
    alt_apo = r_apo - R_earth

    orbit_class = classify_orbit(e, a, r_peri)

    return {
        'name': name,
        'e': e,
        'p': p,
        'a': a,
        'b': b,
        'r_peri': r_peri,
        'r_apo': r_apo,
        'alt_peri': alt_peri,
        'alt_apo': alt_apo,
        'h': h,
        'energy': energy,
        'T': T,
        'v_peri': v_peri,
        'v_apo': v_apo,
        'orbit_class': orbit_class,
        'orbit_type': orbit_type_simple(e),
        'inc': inc,
        'omega': omega,
        'M0': M0,
    }


def format_orbit_info(info):
    """Pretty-format orbital parameters."""
    return (f"  Orbit Type         : {info['orbit_class']}\n"
            f"  Eccentricity (e)   : {info['e']:.6f}\n"
            f"  Semi-latus rectum  : {info['p']:.2f} km\n"
            f"  Semi-major axis (a): {info['a']:.2f} km\n"
            f"  Semi-minor axis (b): {info['b']:.2f} km\n"
            f"  Periapsis radius   : {info['r_peri']:.2f} km  (altitude: {info['alt_peri']:.2f} km)\n"
            f"  Apoapsis radius    : {info['r_apo']:.2f} km  (altitude: {info['alt_apo']:.2f} km)\n"
            f"  Inclination (i)    : {info['inc']:.2f} deg\n"
            f"  Arg of perigee (w) : {info['omega']:.2f} deg\n"
            f"  Mean anomaly (M0)  : {info['M0']:.2f} deg\n"
            f"  Ang. momentum (h)  : {info['h']:.2f} km^2/s\n"
            f"  Specific energy    : {info['energy']:.4f} km^2/s^2\n"
            f"  Orbital period     : {info['T']:.2f} s  ({info['T']/60:.2f} min)  ({info['T']/3600:.4f} hr)\n"
            f"  V at periapsis     : {info['v_peri']:.4f} km/s\n"
            f"  V at apoapsis      : {info['v_apo']:.4f} km/s")


# ORBIT 1 (INITIAL)
e1 = 0.2
p1 = 8000
i1 = 28.5         # inclination (deg)  -- typical Cape Canaveral launch
omega1 = 45.0     # argument of perigee (deg)
M0_1 = 0.0        # mean anomaly at epoch (deg)

theta1, r1, x1, y1 = generate_orbit(e1, p1)
v_r1, v_t1 = velocity(e1, p1, theta1, r1)
orb1 = orbital_analysis("Initial Orbit", e1, p1, i1, omega1, M0_1)

# 3D coordinates
x1_3d, y1_3d, z1_3d = transform_to_3d(x1, y1, i1, omega1)


# ORBIT 2 (FINAL / TARGET)
e2 = 0.6
p2 = 15000
i2 = 51.6         # inclination (deg)  like ISS
omega2 = 120.0    # argument of perigee (deg)
M0_2 = 0.0        # mean anomaly at epoch (deg)

theta2, r2, x2, y2 = generate_orbit(e2, p2)
v_r2, v_t2 = velocity(e2, p2, theta2, r2)
orb2 = orbital_analysis("Final Orbit (Target)", e2, p2, i2, omega2, M0_2)

# 3D coordinates
x2_3d, y2_3d, z2_3d = transform_to_3d(x2, y2, i2, omega2)


# HOHMANN TRANSFER ORBIT
# Transfer from periapsis of orbit 1 to apoapsis of orbit 2
r_peri = orb1['r_peri']     # periapsis of orbit 1 (departure at theta=0)
r_apo = orb2['r_apo']       # apoapsis of orbit 2 (arrival at theta=180)

a_transfer = (r_peri + r_apo) / 2
e_transfer = (r_apo - r_peri) / (r_apo + r_peri)
p_transfer = a_transfer * (1 - e_transfer**2)

# Transfer orbit orientation: departs in orbit 1's plane, arrives in orbit 2's plane
# (plane change is combined with the second burn at apoapsis)
i_t = i1         # transfer orbit starts in orbit 1 plane
omega_t = omega1
orb_t = orbital_analysis("Hohmann Transfer Orbit", e_transfer, p_transfer,
                          i_t, omega_t, 0.0)

# Transfer orbit -- only the active half (periapsis -> apoapsis) in perifocal
theta_t, r_t, x_t, y_t = generate_orbit(e_transfer, p_transfer,
                                          theta_start=0, theta_end=180)

# Full transfer ellipse (for reference)
_, _, x_t_full, y_t_full = generate_orbit(e_transfer, p_transfer)

# 3D coords -- transfer arc with interpolated plane change
frac_t = np.linspace(0, 1, len(x_t))
inc_interp  = i1 + frac_t * (i2 - i1)
omega_interp = omega1 + frac_t * (omega2 - omega1)
_inc_rad  = np.radians(inc_interp)
_omega_rad = np.radians(omega_interp)

_xr = x_t * np.cos(_omega_rad) - y_t * np.sin(_omega_rad)
_yr = x_t * np.sin(_omega_rad) + y_t * np.cos(_omega_rad)
x_t_3d = _xr
y_t_3d = _yr * np.cos(_inc_rad)
z_t_3d = _yr * np.sin(_inc_rad)

# Full transfer ellipse 3D (shown in orbit 1's plane for reference)
x_tf_3d, y_tf_3d, z_tf_3d = transform_to_3d(x_t_full, y_t_full, i1, omega1)


# DELTA-V CALCULATIONS
# At departure (periapsis, theta = 0)
v1_dep  = orb1['v_peri']          # satellite speed on orbit 1
v_t_dep = orb_t['v_peri']         # required on transfer orbit
delta_v1 = v_t_dep - v1_dep       # prograde burn (in-plane)

# At arrival (apoapsis, theta = 180)
v_t_arr = orb_t['v_apo']          # on transfer orbit at apoapsis
v2_arr  = orb2['v_apo']           # required on orbit 2 at apoapsis

# Pure speed change (no plane change)
delta_v2_speed = v2_arr - v_t_arr

# Inclination change
delta_i_deg = abs(i2 - i1)
delta_i_rad = np.radians(delta_i_deg)

# Combined delta-v at arrival (speed change + plane change in one burn)
dv2_combined = np.sqrt(v_t_arr**2 + v2_arr**2 - 2 * v_t_arr * v2_arr * np.cos(delta_i_rad))

# Plane-change-only delta-v (for comparison)
dv_plane_only = 2 * v_t_arr * np.sin(delta_i_rad / 2)

# Total
total_dv = abs(delta_v1) + dv2_combined

# Transfer time = half the period of the transfer ellipse
t_transfer = orb_t['T'] / 2

# Escape velocity at departure
v_escape_dep = np.sqrt(2 * mu / r_peri)


# PRINT FULL RESULTS
dc = "This is a perfectly circular orbit." if orb2['e'] == 0 else \
     "Nearly circular orbit (low eccentricity)." if orb2['e'] < 0.1 else \
     "Moderately elliptical orbit." if orb2['e'] < 0.5 else \
     f"Highly elliptical orbit -- large diff between periapsis and apoapsis.\n         Periapsis: {orb2['r_peri']:.2f} km,  Apoapsis: {orb2['r_apo']:.2f} km\n         Ratio: {orb2['r_apo']/orb2['r_peri']:.2f}"

warn_pa = f"!! WARNING: Periapsis alt ({orb2['alt_peri']:.0f} km) is very low! Drag will be significant." if orb2['alt_peri'] < 200 else \
          "[i] Periapsis is in LEO (Low Earth Orbit) region." if orb2['alt_peri'] < 2000 else \
          "[i] Periapsis is in MEO (Medium Earth Orbit) region." if orb2['alt_peri'] < 35786 else "[i] Periapsis is above GEO altitude."

warn_ap = f"[i] Apoapsis extends beyond GEO ({orb2['alt_apo']:.0f} km)." + ("\n     Resembles Molniya-type highly elliptical orbit." if orb2['e'] > 0.5 else "") if orb2['alt_apo'] > 35786 else ""

output_str = f"""
       ORBITAL TRANSFER SIMULATION


  1. INITIAL ORBIT
{format_orbit_info(orb1)}


  2. FINAL ORBIT (TARGET)
{format_orbit_info(orb2)}


  3. HOHMANN TRANSFER ORBIT
{format_orbit_info(orb_t)}


  4. VELOCITY ANALYSIS
  > At Departure (periapsis of initial orbit, theta = 0 deg):
      V on orbit 1    : {v1_dep:.4f} km/s
      V on transfer   : {v_t_dep:.4f} km/s
      dV1 (burn 1)    : {delta_v1:+.4f} km/s  ({'prograde' if delta_v1 > 0 else 'retrograde'})
      V_escape here   : {v_escape_dep:.4f} km/s

  > At Arrival (apoapsis of final orbit, theta = 180 deg):
      V on transfer   : {v_t_arr:.4f} km/s
      V on orbit 2    : {v2_arr:.4f} km/s
      dV2 speed only  : {delta_v2_speed:+.4f} km/s  ({'prograde' if delta_v2_speed > 0 else 'retrograde'})

  > Plane Change at Arrival:
      Inclination change  : {delta_i_deg:.2f} deg
      dV plane-only       : {dv_plane_only:.4f} km/s
      dV combined (speed+plane): {dv2_combined:.4f} km/s
      Savings vs sep burns: {(abs(delta_v2_speed)+dv_plane_only) - dv2_combined:.4f} km/s


  5. TRANSFER SUMMARY
  Total dV             : {total_dv:.4f} km/s
  Transfer time (TOF)  : {t_transfer:.2f} s  ({t_transfer/60:.2f} min)  ({t_transfer/3600:.4f} hr)

  > Orbit change        : {orb1['orbit_type']} (e={orb1['e']:.4f})  -->  {orb_t['orbit_type']}  -->  {orb2['orbit_type']} (e={orb2['e']:.4f})
  > Periapsis alt change: {orb1['alt_peri']:.2f} km  -->  {orb2['alt_peri']:.2f} km  (d = {orb2['alt_peri']-orb1['alt_peri']:+.2f} km)
  > Apoapsis alt change : {orb1['alt_apo']:.2f} km  -->  {orb2['alt_apo']:.2f} km  (d = {orb2['alt_apo']-orb1['alt_apo']:+.2f} km)
  > Energy change       : {orb1['energy']:.4f}  -->  {orb2['energy']:.4f} km^2/s^2  (d = {orb2['energy'] - orb1['energy']:+.4f})
  > Period change       : {orb1['T']/60:.2f} min  -->  {orb2['T']/60:.2f} min  (d = {(orb2['T']-orb1['T'])/60:+.2f} min)
  > Inclination change  : {i1:.2f} deg  -->  {i2:.2f} deg  (di = {i2-i1:+.2f} deg)


  6. DESTINATION ORBIT CLASSIFICATION
  After the transfer, the satellite is now in:
      >>  {orb2['orbit_class']}
      >>  Eccentricity: {orb2['e']:.4f}
      >>  Inclination : {orb2['inc']:.2f} deg
      >>  Arg perigee : {orb2['omega']:.2f} deg
      >>  {dc}

  {warn_pa}
  {warn_ap}

  Saving plots and outputs to file...
"""

save_dir = r"E:\Coding work\Antigravity\S1\Orbital_Mechanics"
os.makedirs(save_dir, exist_ok=True)
with open(os.path.join(save_dir, "orbital_transfer_results.txt"), "w", encoding="utf-8") as f:
    f.write(output_str)

# PLOT 1: STATIC ORBITS (2D perifocal)
fig1, ax1 = plt.subplots(figsize=(8, 8))

ax1.plot(x1, y1, 'dodgerblue', linewidth=1.5, label="Initial Orbit")
ax1.plot(x2, y2, 'limegreen', linewidth=1.5, label="Final Orbit")

ax1.plot(x_t, y_t, '-', color='orange', linewidth=2, label="Hohmann Transfer (active)")

# Mark departure and arrival points
ax1.plot(x_t[0], y_t[0], 'v', color='red', markersize=10,
         label=f"Departure (dV1={delta_v1:+.3f} km/s)")
ax1.plot(x_t[-1], y_t[-1], '^', color='magenta', markersize=10,
         label=f"Arrival (dV2_comb={dv2_combined:.3f} km/s)")

# Draw Earth
earth_circle = plt.Circle((0, 0), R_earth, color='deepskyblue', alpha=0.3,
                           label="Earth (to scale)")
ax1.add_patch(earth_circle)
ax1.scatter(0, 0, color='blue', s=30, zorder=5)

# Annotations
ax1.annotate(f'Periapsis\n{orb1["r_peri"]:.0f} km',
             xy=(x_t[0], y_t[0]), xytext=(x_t[0]+2000, y_t[0]+3000),
             fontsize=8, arrowprops=dict(arrowstyle='->', color='red'),
             color='red', fontweight='bold')
ax1.annotate(f'Apoapsis\n{orb2["r_apo"]:.0f} km',
             xy=(x_t[-1], y_t[-1]), xytext=(x_t[-1]-8000, y_t[-1]+3000),
             fontsize=8, arrowprops=dict(arrowstyle='->', color='magenta'),
             color='magenta', fontweight='bold')

ax1.set_xlabel("X (km)")
ax1.set_ylabel("Y (km)")
ax1.set_title("Hohmann Transfer Between Two Elliptical Orbits (Perifocal View)")
ax1.legend(loc='upper right', fontsize=7)
ax1.set_aspect("equal")
ax1.grid(alpha=0.3)

plt.tight_layout()
fig1.savefig(os.path.join(save_dir, "plot1_perifocal.png"), dpi=300)
plt.show(block=True)
print("Plot 1 saved.")


# PLOT 2: VELOCITY VECTORS
fig2, ax2 = plt.subplots(figsize=(8, 8))

ax2.plot(x1, y1, 'dodgerblue', linewidth=1.5, label="Initial Orbit")

# Plot velocity vectors (sample points)
step = 200
for idx in range(0, len(x1), step):
    vx = v_r1[idx]*np.cos(theta1[idx]) - v_t1[idx]*np.sin(theta1[idx])
    vy = v_r1[idx]*np.sin(theta1[idx]) + v_t1[idx]*np.cos(theta1[idx])
    v_mag = np.sqrt(vx**2 + vy**2)

    ax2.quiver(x1[idx], y1[idx], vx, vy, scale=5000, color='red', width=0.004)
    ax2.annotate(f'{v_mag:.2f} km/s', xy=(x1[idx], y1[idx]),
                 xytext=(5, 10), textcoords='offset points',
                 fontsize=6, color='darkred')

ax2.scatter(0, 0, color='deepskyblue', s=200, zorder=5, edgecolors='blue')

ax2.set_title("Velocity Vectors Along Initial Orbit")
ax2.set_aspect("equal")
ax2.grid(alpha=0.3)
ax2.legend()

plt.tight_layout()
fig2.savefig(os.path.join(save_dir, "plot2_velocity_vectors.png"), dpi=300)
plt.show(block=True)
print("Plot 2 saved.")


# PLOT 3: VELOCITY MAGNITUDE PROFILE
fig3, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, t, vr, vt, orb, title in zip(
    axes, 
    [theta1, theta2], [v_r1, v_r2], [v_t1, v_t2], [orb1, orb2], 
    ["Initial Orbit", "Final Orbit (Target)"]):
    
    ax.plot(np.degrees(t), np.sqrt(vr**2 + vt**2), 'dodgerblue' if 'Init' in title else 'limegreen', lw=1.5, label='|V|')
    ax.plot(np.degrees(t), vr, '--', color='orange', lw=1, label='V_r (radial)')
    ax.plot(np.degrees(t), vt, '--', color='green', lw=1, label='V_theta (tangential)')
    ax.axhline(y=orb['v_peri'], color='red', ls=':', alpha=0.5, label=f"V_peri = {orb['v_peri']:.3f}")
    ax.axhline(y=orb['v_apo'], color='purple', ls=':', alpha=0.5, label=f"V_apo = {orb['v_apo']:.3f}")
    ax.set(xlabel="True Anomaly theta (deg)", ylabel="Velocity (km/s)", title=f"Velocity Profile -- {title}")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

plt.tight_layout()
fig3.savefig(os.path.join(save_dir, "plot3_velocity_profile.png"), dpi=300)
plt.show(block=True)
print("Plot 3 saved.")

# ANIMATION: DUAL-PANEL
#   Left  = 2D perifocal
#   Right = 3D orbital-elements diagram
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Phase 1 -- coast on initial orbit (one full revolution, perifocal)
phase1_x = x1
phase1_y = y1

# Phase 2 -- transfer arc (perifocal)
phase2_x = x_t
phase2_y = y_t

# Phase 3 -- coast on final orbit starting from theta=180 (perifocal)
theta3, r3, x3_phase, y3_phase = generate_orbit(e2, p2, theta_start=180, theta_end=540)

# 3D transforms per phase
p1_x3d, p1_y3d, p1_z3d = transform_to_3d(phase1_x, phase1_y, i1, omega1)

_frac2 = np.linspace(0, 1, len(phase2_x))
_i2 = np.radians(i1 + _frac2 * (i2 - i1))
_w2 = np.radians(omega1 + _frac2 * (omega2 - omega1))
_xr2 = phase2_x * np.cos(_w2) - phase2_y * np.sin(_w2)
_yr2 = phase2_x * np.sin(_w2) + phase2_y * np.cos(_w2)
p2_x3d = _xr2
p2_y3d = _yr2 * np.cos(_i2)
p2_z3d = _yr2 * np.sin(_i2)

p3_x3d, p3_y3d, p3_z3d = transform_to_3d(x3_phase, y3_phase, i2, omega2)

# Concatenate -- perifocal
anim_x_pf = np.concatenate([phase1_x, phase2_x, x3_phase])
anim_y_pf = np.concatenate([phase1_y, phase2_y, y3_phase])

# Concatenate -- 3D
anim_x3d = np.concatenate([p1_x3d, p2_x3d, p3_x3d])
anim_y3d = np.concatenate([p1_y3d, p2_y3d, p3_y3d])
anim_z3d = np.concatenate([p1_z3d, p2_z3d, p3_z3d])

n1 = len(phase1_x)
n2 = len(phase2_x)
n3 = len(x3_phase)

skip = 5
anim_x_pf = anim_x_pf[::skip]
anim_y_pf = anim_y_pf[::skip]
anim_x3d  = anim_x3d[::skip]
anim_y3d  = anim_y3d[::skip]
anim_z3d  = anim_z3d[::skip]
n1_ds = n1 // skip
n2_ds = n2 // skip
total_frames = len(anim_x3d)

# FIGURE: left = 2-D perifocal | right = 2-D diagram
fig4 = plt.figure(figsize=(19, 9))
ax4a = fig4.add_subplot(1, 2, 1)
ax4b = fig4.add_subplot(1, 2, 2, projection='3d')

# LEFT PANEL (2D perifocal animation)
ax4a.set_title("Orbit Transfer Animation", fontsize=11, fontweight='bold')
ax4a.plot(x1, y1, 'dodgerblue', linewidth=1, alpha=0.4, label="Initial Orbit")
ax4a.plot(x2, y2, 'limegreen', linewidth=1, alpha=0.4, label="Final Orbit")
ax4a.plot(x_t, y_t, '--', color='orange', linewidth=1, alpha=0.4, label="Transfer Arc")

earth_2d = plt.Circle((0, 0), R_earth, color='deepskyblue', alpha=0.3)
ax4a.add_patch(earth_2d)
ax4a.scatter(0, 0, color='blue', s=30, zorder=5, label="Earth")

ax4a.plot(x_t[0], y_t[0], 'v', color='red', markersize=8, alpha=0.6)
ax4a.plot(x_t[-1], y_t[-1], '^', color='magenta', markersize=8, alpha=0.6)

all_x = np.concatenate([x1, x2, x_t])
all_y = np.concatenate([y1, y2, y_t])
margin = 1.15
ax4a.set_xlim(min(all_x)*margin, max(all_x)*margin)
ax4a.set_ylim(min(all_y)*margin, max(all_y)*margin)
ax4a.set_xlabel("X (km)")
ax4a.set_ylabel("Y (km)")
ax4a.set_aspect("equal")
ax4a.legend(loc='upper right', fontsize=7)
ax4a.grid(alpha=0.3)

# Animated 2-D elements
trail_2d, = ax4a.plot([], [], '-', linewidth=2, alpha=0.7)
sat_2d,   = ax4a.plot([], [], 'ro', markersize=8, zorder=10)
phase_text = ax4a.text(0.02, 0.95, '', transform=ax4a.transAxes,
                       fontsize=10, fontweight='bold',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
info_text  = ax4a.text(0.02, 0.02, '', transform=ax4a.transAxes,
                       fontsize=7.5, fontfamily='monospace',
                       verticalalignment='bottom',
                       bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

# RIGHT PANEL (isometric orbital-elements diagram)
ax4b.set_title("Orbital Elements Diagram", fontsize=11, fontweight='bold')
try:
    ax4b.set_aspect('equal')
except NotImplementedError:
    pass
ax4b.axis('off')
ax4b.set_facecolor('none')

# scaling for the diagram
r_plane = orb2['r_apo'] * 0.45       # radius of the plane disks

# Equatorial reference plane (gray filled ellipse)
_tp = np.linspace(0, 2*np.pi, 120)
eq_x3 = r_plane * np.cos(_tp)
eq_y3 = r_plane * np.sin(_tp)
eq_z3 = np.zeros_like(_tp)
ax4b.add_collection3d(Poly3DCollection([list(zip(eq_x3, eq_y3, eq_z3))], color='silver', alpha=0.25))
ax4b.plot(eq_x3, eq_y3, eq_z3, color='gray', lw=0.7, alpha=0.5, zorder=1)

# Label the equatorial plane
ax4b.text(r_plane*0.55, -r_plane*0.7, 0, 'Plane of\nreference', fontsize=8, color='gray',
          ha='center', style='italic', alpha=0.8, zorder=2)

# Orbital plane for orbit 1 (yellow/gold, tilted by i1)
i1r = np.radians(i1)
op1_x3 = r_plane * np.cos(_tp)
op1_y3 = r_plane * np.sin(_tp) * np.cos(i1r)
op1_z3 = r_plane * np.sin(_tp) * np.sin(i1r)
ax4b.add_collection3d(Poly3DCollection([list(zip(op1_x3, op1_y3, op1_z3))], color='gold', alpha=0.15))
ax4b.plot(op1_x3, op1_y3, op1_z3, color='goldenrod', lw=0.7, alpha=0.5, zorder=2)

# Label the orbital plane
ax4b.text(-r_plane*0.4, r_plane*0.6, r_plane*0.3*np.sin(i1r), 'Orbital\nplane', fontsize=8, color='goldenrod',
          ha='center', style='italic', fontweight='bold', zorder=3)

# Orbital plane for orbit 2 (green tint, tilted by i2)
i2r = np.radians(i2)
op2_x3 = r_plane * np.cos(_tp)
op2_y3 = r_plane * np.sin(_tp) * np.cos(i2r)
op2_z3 = r_plane * np.sin(_tp) * np.sin(i2r)
ax4b.add_collection3d(Poly3DCollection([list(zip(op2_x3, op2_y3, op2_z3))], color='limegreen', alpha=0.07))
ax4b.plot(op2_x3, op2_y3, op2_z3, color='green', lw=0.5, alpha=0.3, zorder=2)

# Orbit curves
ax4b.plot(x1_3d, y1_3d, z1_3d, 'dodgerblue', lw=1.5, alpha=0.4, zorder=4)
ax4b.plot(x2_3d, y2_3d, z2_3d, 'limegreen', lw=1.5, alpha=0.4, zorder=4)
ax4b.plot(x_t_3d, y_t_3d, z_t_3d, '--', color='orange', lw=1.2, alpha=0.4, zorder=4)

# Earth at center
earth_r2d = R_earth * 0.8  # slightly smaller for visual clarity
_ec = np.linspace(0, 2*np.pi, 60)
_ex3 = earth_r2d * np.cos(_ec)
_ey3 = earth_r2d * np.sin(_ec)
_ez3 = np.zeros_like(_ec)
ax4b.add_collection3d(Poly3DCollection([list(zip(_ex3, _ey3, _ez3))], color='deepskyblue', alpha=0.4))
ax4b.plot(_ex3, _ey3, _ez3, color='blue', lw=0.5, zorder=5)
cx0, cy0, cz0 = 0, 0, 0
ax4b.plot([0], [0], [0], 'o', color='blue', markersize=4, zorder=6)

# Reference direction arrow (red, +X)
ref_len = r_plane * 1.15
ax4b.quiver(0, 0, 0, ref_len, 0, 0, length=1, arrow_length_ratio=0.1, color='red', lw=1.5, zorder=6)
ax4b.text(ref_len*1.08, 0, 0, 'Reference\ndirection', fontsize=7, color='red',
          ha='left', va='center', fontweight='bold', zorder=7)

# Ascending node line + marker
node_r = r_plane * 0.9
ax4b.plot([0, node_r], [0, 0], [0, 0], color='darkgreen',
          ls='-.', lw=1.2, alpha=0.6, zorder=3)
ax4b.text(node_r*0.85, -r_plane*0.06, 0, 'Ascending\nnode',
          fontsize=7, color='darkgreen', ha='center', va='top',
          fontweight='bold', zorder=7)


# Satellite point (black dot)
sat_diag, = ax4b.plot([0], [0], [0], 'ko', markersize=8, zorder=10)
sat_label = ax4b.text(0, 0, r_plane*0.05, 'Satellite', fontsize=7,
                      ha='center', va='bottom', fontweight='bold', zorder=10)

# Trail 3D
trail_3d, = ax4b.plot([], [], [], '-', linewidth=2.5, alpha=0.8, zorder=9)

# Angle readout box
ang_text = ax4b.text2D(0.97, 0.03, '', transform=ax4b.transAxes,
                     fontsize=8.5, fontfamily='monospace',
                     ha='right', va='bottom',
                     bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow',
                               alpha=0.85, edgecolor='gray'),
                     zorder=10)

# Set axis limits for right panel
pad = r_plane * 0.35
ax4b.set_xlim([-r_plane-pad, r_plane+pad])
ax4b.set_ylim([-r_plane-pad, r_plane+pad])
ax4b.set_zlim([-r_plane-pad, r_plane+pad])

# Trail length
trail_len = 80

def update(frame):
    # perifocal coords
    cx_pf = anim_x_pf[frame]
    cy_pf = anim_y_pf[frame]
    r_now = np.sqrt(cx_pf**2 + cy_pf**2)
    theta_now = np.arctan2(cy_pf, cx_pf)

    # LEFT PANEL (2D perifocal)
    sat_2d.set_data([cx_pf], [cy_pf])
    start = max(0, frame - trail_len)
    trail_2d.set_data(anim_x_pf[start:frame+1], anim_y_pf[start:frame+1])
    
    # 3D Trail
    trail_3d.set_data(anim_x3d[start:frame+1], anim_y3d[start:frame+1])
    trail_3d.set_3d_properties(anim_z3d[start:frame+1])

    # phase logic
    if frame < n1_ds:
        phase_text.set_text("Phase 1: Coasting on Initial Orbit")
        trail_2d.set_color('dodgerblue')
        v_now  = np.sqrt(mu * (2/r_now - 1/orb1['a']))
        e_now  = e1
        i_now  = i1
        w_now  = omega1
        orb_label = f"Orbit 1 ({orb1['orbit_type']}, e={e1:.2f})"
    elif frame < n1_ds + n2_ds:
        phase_text.set_text("Phase 2: Hohmann Transfer + Plane Change")
        trail_2d.set_color('orange')
        v_now  = np.sqrt(mu * (2/r_now - 1/orb_t['a']))
        e_now  = e_transfer
        frac   = (frame - n1_ds) / max(n2_ds, 1)
        i_now  = i1 + frac * (i2 - i1)
        w_now  = omega1 + frac * (omega2 - omega1)
        orb_label = f"Transfer (e={e_transfer:.4f})"
    else:
        phase_text.set_text("Phase 3: Coasting on Final Orbit")
        trail_2d.set_color('limegreen')
        v_now  = np.sqrt(mu * (2/r_now - 1/orb2['a']))
        e_now  = e2
        i_now  = i2
        w_now  = omega2
        orb_label = f"Orbit 2 ({orb2['orbit_type']}, e={e2:.2f})"

    M_now = np.degrees(true_to_mean_anomaly(theta_now, e_now)) % 360
    theta_deg_now = np.degrees(theta_now) % 360
    
    trail_3d.set_color(trail_2d.get_color())

    info_text.set_text(
        f"r = {r_now:.0f} km | alt = {r_now - R_earth:.0f} km\n"
        f"|V| = {v_now:.3f} km/s\n"
        f"{orb_label}"
    )

    # RIGHT PANEL (3D diagram updates)
    sx, sy, sz = anim_x3d[frame], anim_y3d[frame], anim_z3d[frame]
    sat_diag.set_data([sx], [sy])
    sat_diag.set_3d_properties([sz])
    sat_label.set_position((sx, sy))
    sat_label.set_3d_properties(sz + r_plane*0.04, zdir='z')

    # Angle readout text box
    ang_text.set_text(
        f"True anomaly  v = {theta_deg_now:6.1f} deg\n"
        f"Mean anomaly  M = {M_now:6.1f} deg\n"
        f"Arg perigee   w = {w_now:6.1f} deg\n"
        f"Inclination   i = {i_now:6.1f} deg\n"
        f"Eccentricity  e = {e_now:.4f}"
    )

    return (sat_2d, trail_2d, phase_text, info_text,
            sat_diag, sat_label,
            trail_3d, ang_text)

ani = FuncAnimation(fig4, update, frames=total_frames,
                    interval=25, blit=False)

fig4.suptitle("Orbit Transfer -- Animation + Orbital Elements", fontsize=14, fontweight='bold', y=0.99)
plt.tight_layout(rect=[0, 0, 1, 0.97])
print("\n[Plot 4] Opening animation window...")
plt.show(block=False)
plt.pause(0.5)

ans = input("The animation is playing. Do you want to save it as MP4? (It will freeze momentarily while saving) (y/n): ").strip().lower()
if ans == 'y':
    print(f"Saving animation to orbit_animation.mp4 (this may take a few seconds)...")
    try:
        ani.save(os.path.join(save_dir, "orbit_animation.mp4"), writer="ffmpeg", fps=24)
        print("Saved animation.")
    except Exception as e:
        print(f"Error saving animation (is ffmpeg installed?): {e}")

print("(Animation will continue to play. Close the window to finish running the script)")
plt.show(block=True)
print(f"\nAll operations completed. Saved files (if any) are located in: {save_dir}")
