# Vieri's part: RK4 for a Hohmann transfer
# Uses functions, loops, NumPy arrays, and basic plotting.
# Units: distance in km, time in seconds, velocity in km/s.

import numpy as np
import matplotlib.pyplot as plt


# Differential equations for motion under gravity.
# u stores [x position, y position, x velocity, y velocity].
def orbit(u, mu):
    x = u[0]
    y = u[1]
    vx = u[2]
    vy = u[3]
    r = np.sqrt(x**2 + y**2)

    ax = -mu * x / r**3
    ay = -mu * y / r**3
    return np.array([vx, vy, ax, ay])


# Classical fourth-order Runge-Kutta method.
def rk4(mu, r1, r2, N):
    a = (r1 + r2) / 2
    T = np.pi * np.sqrt(a**3 / mu)
    h = T / N

    t = np.linspace(0, T, N + 1)
    u = np.zeros((N + 1, 4))

    # First burn: start at (r1, 0) with the transfer-orbit velocity.
    v_start = np.sqrt(mu * (2 / r1 - 1 / a))
    u[0] = [r1, 0, 0, v_start]

    for i in range(N):
        k1 = orbit(u[i], mu)
        k2 = orbit(u[i] + h * k1 / 2, mu)
        k3 = orbit(u[i] + h * k2 / 2, mu)
        k4 = orbit(u[i] + h * k3, mu)

        u[i + 1] = u[i] + h * (k1 + 2*k2 + 2*k3 + k4) / 6

    return t, u


def main():
    # Example inputs. Radii are measured from Earth's center.
    # Switch r1 and r2 for an inward transfer.
    mu = 398600.4418
    r1 = 7000.0
    r2 = 14000.0
    N = 200

    a = (r1 + r2) / 2
    T = np.pi * np.sqrt(a**3 / mu)

    # Ideal burn calculations from circular velocity and vis-viva.
    # Positive delta-v speeds up; negative delta-v slows down.
    v_circular1 = np.sqrt(mu / r1)
    v_transfer1 = np.sqrt(mu * (2 / r1 - 1 / a))
    v_transfer2 = np.sqrt(mu * (2 / r2 - 1 / a))
    v_circular2 = np.sqrt(mu / r2)
    dv1 = v_transfer1 - v_circular1
    dv2 = v_circular2 - v_transfer2

    t, u = rk4(mu, r1, r2, N)

    # At half the ellipse period, the exact arrival point is (-r2, 0).
    # This reference comes from orbital geometry, not another numerical solver.
    error = np.sqrt((u[N, 0] + r2)**2 + u[N, 1]**2)

    # Apply the second burn along the local prograde tangent.
    # Keep numerical arrival error instead of replacing it with the exact state.
    x_end = u[N, 0]
    y_end = u[N, 1]
    r_end = np.sqrt(x_end**2 + y_end**2)
    vx_after = u[N, 2] + dv2 * (-y_end / r_end)
    vy_after = u[N, 3] + dv2 * (x_end / r_end)
    radial_velocity = (x_end * vx_after + y_end * vy_after) / r_end

    print('Time of flight (s):', T)
    print('First burn (km/s):', dv1)
    print('Second burn (km/s):', dv2)
    print('Total delta-v (km/s):', abs(dv1) + abs(dv2))
    print('Arrival position error (km):', error)
    print('Radial velocity after second burn (km/s):', radial_velocity)

    # Convergence study: double N, which halves h.
    steps = [25, 50, 100, 200, 400, 800]
    errors = np.zeros(len(steps))
    h_values = np.zeros(len(steps))
    orders = np.zeros(len(steps))
    orders[0] = np.nan  # No previous result for the first row.

    print('\nSteps      h (s)       Error (km)      Order')
    for j in range(len(steps)):
        n = steps[j]
        times, states = rk4(mu, r1, r2, n)
        h_values[j] = T / n
        errors[j] = np.sqrt((states[n, 0] + r2)**2 + states[n, 1]**2)

        if j > 0:
            orders[j] = np.log(errors[j - 1] / errors[j]) / np.log(2)

        print('%5d   %10.4f   %12.5e   %6.3f' %
              (n, h_values[j], errors[j], orders[j]))

    # Save tables in the folder from which the program is run.
    table = np.zeros((len(steps), 4))
    table[:, 0] = steps
    table[:, 1] = h_values
    table[:, 2] = errors
    table[:, 3] = orders
    np.savetxt('basic_convergence.csv', table, delimiter=',',
               header='steps,h_seconds,endpoint_error_km,order', comments='')

    trajectory = np.zeros((N + 1, 5))
    trajectory[:, 0] = t
    trajectory[:, 1:] = u
    np.savetxt('basic_trajectory.csv', trajectory, delimiter=',',
               header='time_s,x_km,y_km,vx_km_s,vy_km_s', comments='')

    # Plot the two circular orbits and the numerical transfer path.
    theta = np.linspace(0, 2*np.pi, 300)
    plt.figure()
    plt.plot(r1*np.cos(theta), r1*np.sin(theta), '--', label='Initial orbit')
    plt.plot(r2*np.cos(theta), r2*np.sin(theta), '--', label='Final orbit')
    plt.plot(u[:, 0], u[:, 1], label='RK4 transfer')
    plt.plot(-r2, 0, 'ko', label='Exact arrival')
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')
    plt.title('Hohmann transfer using RK4')
    plt.axis('equal')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('basic_orbit.png', dpi=150)

    plt.figure()
    plt.loglog(h_values, errors, 'o-', label='RK4 error')
    reference = errors[-1] * (h_values / h_values[-1])**4
    plt.loglog(h_values, reference, '--', label='Fourth-order slope')
    plt.xlabel('Time step h (s)')
    plt.ylabel('Arrival position error (km)')
    plt.title('RK4 convergence')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('basic_convergence.png', dpi=150)
    plt.show()



if __name__ == '__main__':
    main()
