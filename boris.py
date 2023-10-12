import numpy as np
import matplotlib.pyplot as plt
import time
# Initilization
start_time = time.time()
time_interval = 1.0
num_iter = 10 
dt = time_interval/num_iter
pos = np.zeros((3, num_iter))
vel = np.zeros((3, num_iter))
vel[:,0] = [0.1, 0.0, 0.02]
# Given Electric and Magnetic Field
def magnetic_field(pos):
     return np.array([0.0, 0.0, 1.0])
def electric_field(pos):
     return np.array([0.0, 0.0, 0.0])
# Relativistic Boris Solver
for i in range(0, num_iter-1, 1):
     pos[:,i+1] = pos[:,i] + dt*vel[:,i]
     U_vec = vel[:,i]/np.sqrt(1-np.linalg.norm(vel[:,i])**2)
     t_vec =  (dt/2.0)*magnetic_field(pos[:,i+1])/np.sqrt(1+np.linalg.norm(U_vec)**2)
     s_vec =  (dt/2.0)*electric_field(pos[:,i+1])
     u_vec = U_vec + s_vec
     w_vec = u_vec + np.cross(u_vec, t_vec)
     u_vec = u_vec + (2.0/(1.0+np.linalg.norm(t_vec)**2))*np.cross(w_vec, t_vec)
     vel[:,i+1]= (u_vec + s_vec)/np.sqrt(1+np.linalg.norm(u_vec)**2)

#Exact Solution
time_set = np.linspace(0,time_interval,num_iter-1)
omega = -1.0*np.sqrt(1-0.1**2-0.02**2)
x_exact = (0.1/omega)*np.sin(omega*time_set)
y_exact = (0.1/omega)*(1-np.cos(omega*time_set))
z_exact = 0.02*time_set

end_time = time.time()
mse_error = np.sqrt(np.linalg.norm(x_exact-pos[0,0:num_iter-1])**2+np.linalg.norm(y_exact-pos[1,0:num_iter-1])**2+np.linalg.norm(z_exact-pos[2,0:num_iter-1])**2)
print(end_time-start_time)
print(mse_error)

# Plot
plt.plot(pos[0,:],pos[1,:])
plt.plot(x_exact[:],y_exact[:])
plt.show()

