import numpy as np
import matplotlib.pyplot as plt

# Stellar and planetary parameters
R_star = 1.1   # in solar radii
R_planet = 0.41 # in solar radii 
a = 2.77      # semi-major axis in solar radii

# Range of eccentricities to test
eccentricities = np.linspace(0, 0.55, 100)

# Different arguments of periastron to compare (degrees)
omega_values_deg = np.linspace(0, 180, 100)
omega_values_rad = np.radians(omega_values_deg)

# Inclination range (degrees)
inclinations_deg = np.linspace(0, 180, 1000)  # typical transit inclinations
inclinations_rad = np.radians(inclinations_deg)

#Calculate the probability of seeing a secondary eclipse
loop_counter = 0
no_eclipse_counter = 0
for i in range(len(inclinations_rad)):
    for w in range(len(omega_values_rad)):
        for e in range(len(eccentricities)):
            loop_counter += 1
            eclipse_value = (a/R_star)*np.cos(inclinations_rad[i])*((1-((eccentricities[e])**2))/(1-((eccentricities[e])*np.sin(omega_values_rad[w]))))
            transit_value = (a/R_star)*np.cos(inclinations_rad[i])*((1-((eccentricities[e])**2))/(1+((eccentricities[e])*np.sin(omega_values_rad[w]))))
            if (transit_value > (1 + (R_planet/R_star))):
                if (eclipse_value > (1 + (R_planet/R_star))):
                    no_eclipse_counter +=1
probability_no_secondary_eclipse = no_eclipse_counter/loop_counter
print(probability_no_secondary_eclipse)
exit()

# Example: compute probability for each eccentricity, omega, averaging over inclinations
plt.figure(figsize=(8,5))
for omega, omega_deg in zip(omega_values_rad, omega_values_deg):
    probs = []
    for e in eccentricities:
        # Treat inclinations as uniform in cos(i)
        misses = [P_no_secondary_incl(R_star, R_planet, a, e, omega, i) for i in inclinations_rad]
        prob = np.mean(misses)  # fraction of inclinations where secondary is missed
        probs.append(prob)
    plt.plot(eccentricities, probs, label=f'ω = {omega_deg}°')

plt.xlabel('Eccentricity')
plt.ylabel('Probability of missing secondary eclipse')
plt.title('Probability of Missing Secondary Eclipse vs Eccentricity\n(inclination included)')
plt.grid(True)
plt.legend()
plt.show()
