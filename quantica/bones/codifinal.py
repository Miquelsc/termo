
# coding: utf-8

# In[1]:


from __future__ import division, print_function

import numpy as np
import math
import matplotlib.pyplot as plt
import time
import os




# In[2]:


#
# Equacions trascendentals pels casos parell i senar
# i els diferents valors de l'energia:
#		 0 < E < V0  --> *_l
#		 V0 < E, E>0 --> *_g
#        V0 < E < 0  --> *_n
#
def _even_l(E):
	return (math.sqrt(V0-E))*math.tanh(k*(math.sqrt(V0-E))*l)*math.sin(k*(math.sqrt(E))*(L-l)) + 		(math.sqrt(E))*math.cos(k*(math.sqrt(E))*(L-l))

def _even_g(E):
	return (math.sqrt(E-V0))*math.sin(k*(math.sqrt(E-V0))*l)*math.sin(k*(math.sqrt(E))*(L-l)) - 		(math.sqrt(E))*math.cos(k*(math.sqrt(E-V0))*l)*math.cos(k*(math.sqrt(E))*(L-l))

def _even_n(E):
    return math.sqrt(E-V0)*math.tanh(k*math.sqrt(-E)*(l-L))*math.sin(k*math.sqrt(E-V0)*l) +         math.sqrt(-E)*math.cos(k*math.sqrt(E-V0)*l)

def _odd_l(E):
	return (math.sqrt(V0-E))*math.sin(k*(math.sqrt(E))*(L-l)) + 		(math.sqrt(E))*math.tanh(k*(math.sqrt(V0-E))*l)*math.cos(k*(math.sqrt(E))*(L-l))

def _odd_g(E):
	return (math.sqrt(E-V0))*math.cos(k*(math.sqrt(E-V0))*l)*math.sin(k*(math.sqrt(E))*(L-l)) + 		(math.sqrt(E))*math.sin(k*(math.sqrt(E-V0))*l)*math.cos(k*(math.sqrt(E))*(L-l))
    
def _odd_n(E):
    return math.sqrt(E-V0)*math.tanh(k*math.sqrt(-E)*(l-L))*math.cos(k*math.sqrt(E-V0)*l) -         math.sqrt(-E)*math.sin(k*math.sqrt(E-V0)*l)

# Funcions per guardar les energies en un arxiu de text i per llegir-les d'aquest

def save_energies(E):
	with open(FILE_ENERGIES, 'w') as outf:
		for j in range(len(E)):
			outf.write('%d\t%.6g\n' % (j, E[j])) # Compte xifres significatives.

def read_energies(file_name):
	Ep = []
	with open(file_name) as f:
		for line in f:
			Ep.append(float(line.split('\t')[1].strip()))
	return np.array(Ep)

# Funcio per trobar els valors propis de les energies
Eparell=[]
Esenar=[]

def find_roots():
	print(V0)
	E0 = min(0,V0)
	E = E0 + dE
    
	Ep = [] # energia dels estats
	j = 0 # numero d'estats

    # Per diferenciar el cas del pou al de la barrera
	if V0>0 or V0==0:
		last_even, last_odd = _even_l(0), _odd_l(0)
	elif V0<0:    
		last_even, last_odd = _even_n(E0), _odd_n(E0)

	print('Start root finding...', end=' ')
	start = current_milli_time()

	while E < V0 and j < N:
		e, o = _even_l(E), _odd_l(E)
		
    

		if e * last_even < 0: # canvi de signe, arrel trobada
			Ep.append(E)
			Eparell.append(E)
           
			j+=1
        
		if o * last_odd < 0: 
			Ep.append(E)
			Esenar.append(E)
            
			j+=1

		last_even, last_odd = e, o
		E += dE

	while E<0 and j < N:
		e, o = _even_n(E), _odd_n(E)
		

		if e * last_even < 0: # canvi de signe, arrel trobada
			Ep.append(E)
            
			Eparell.append(E)
			j+=1

		if o * last_odd < 0: # canvi de signe, arrel trobada
			Ep.append(E)
			Esenar.append(E)
			j+=1

		last_even, last_odd = e, o
		E += dE

	last_even, last_odd = _even_g(max(0,V0)), _odd_g(max(0,V0))
	while j < N:
		
		e, o = _even_g(E), _odd_g(E)

		if e * last_even < 0: # canvi de signe, arrel trobada
			Ep.append(E)
			Eparell.append(E)
			j+=1

		if o * last_odd < 0: # canvi de signe, arrel trobada
			Ep.append(E)
			Esenar.append(E)
			j+=1

		last_even, last_odd = e, o
		E += dE

	print('OK (%.2f s)' % ((current_milli_time() - start) / 1000))

	return sorted(Ep)



# In[42]:





# In[3]:


###################################################################################
###################################################################################
####################									   ########################
####################				 PART II			   ########################
####################   Definicio de les funcions propies   ########################
####################									   ########################
###################################################################################
###################################################################################

# Definició de les funcions d'ona pel cas parell i pel senar

def _phi_even_l(reg, E, x):
	if reg == 1:
		return np.sin(k*(np.sqrt(E))*(x+L))
	elif reg == 2:
		return np.sin(k*(np.sqrt(E))*(L-l))*np.cosh(k*(np.sqrt(V0-E))*x)/(np.cosh(k*(np.sqrt(V0-E))*l))
	elif reg == 3:
		return -np.sin(k*(np.sqrt(E))*(x-L))

def _phi_even_g(reg, E, x):
	if reg == 1:
		return np.sin(k*(np.sqrt(E))*(x+L))
	elif reg == 2:
		return np.sin(k*(np.sqrt(E))*(L-l))*np.cos(k*(np.sqrt(E-V0))*x)/(np.cos(k*(np.sqrt(E-V0))*l))
	elif reg == 3:
		return -np.sin(k*(np.sqrt(E))*(x-L))
    
def _phi_even_n(reg, E, x):
	if reg == 1:
		return np.sinh(k*np.sqrt(-E)*(x+L))
	elif reg == 2:
		return np.sinh(k*(np.sqrt(-E))*(L-l))*np.cos(k*(np.sqrt(E-V0))*x)/(np.cos(k*(np.sqrt(E-V0))*l))
	elif reg == 3:
		return -np.sinh(k*(np.sqrt(-E))*(x-L))

def _phi_odd_l(reg, E, x):
	if reg == 1:
		return np.sin(k*(np.sqrt(E))*(x+L))
	elif reg == 2:
		return -np.sin(k*(np.sqrt(E))*(L-l))*np.sinh(k*(np.sqrt(V0-E))*x)/(np.sinh(k*(np.sqrt(V0-E))*l))
	elif reg == 3:
		return np.sin(k*(np.sqrt(E))*(x-L))

def _phi_odd_g(reg, E, x):
	if reg == 1:
		return np.sin(k*(np.sqrt(E))*(x+L))
	elif reg == 2:
		return -np.sin(k*(np.sqrt(E))*(L-l))*np.sin(k*(np.sqrt(E-V0))*x)/(np.sin(k*(np.sqrt(E-V0))*l))
	elif reg == 3:
		return np.sin(k*(np.sqrt(E))*(x-L))
    
def _phi_odd_n(reg, E, x):
	if reg == 1:
		return np.sinh(k*(np.sqrt(-E))*(x+L))
	elif reg == 2:
		return -np.sinh(k*(np.sqrt(-E))*(L-l))*np.sin(k*(np.sqrt(E-V0))*x)/(np.sin(k*(np.sqrt(E-V0))*l))
	elif reg == 3:
		return np.sinh(k*(np.sqrt(-E))*(x-L))

def phi_odd(reg, E, x):
	if E<V0:
		return _phi_odd_l(reg, E, x)
	elif E>0:
		return _phi_odd_g(reg, E, x)
	else:
		return _phi_odd_n(reg, E, x)

def phi_even(reg, E, x):
	if E<V0:
		return _phi_even_l(reg, E, x)
	elif E>0:
		return _phi_even_g(reg, E, x)
	else:
		return _phi_even_n(reg, E, x)


def evaluate_wave_function(Ep):
	# matriu que contindrà totes les funcions propies, cadascuna en una fila:
	PHI = np.zeros((N, Nx))

	# defineix les tres regions diferents de x
	x1, x2, x3 = np.linspace(-L, -l, N1), np.linspace(-l, l, N2), np.linspace(l, L, N3) 

	for j in range(N): # bucle en tots els estats
		E = Ep[j] 
		if j % 2 == 0:
			PHI[j,:N1] = phi_even(1, E, x1)
			PHI[j,N1:N2 + N1] = phi_even(2, E, x2)
			PHI[j,N1 + N2:N3 + N2 + N1] = phi_even(3, E, x3)
		else:
			PHI[j,:N1] = phi_odd(1, E, x1)
			PHI[j,N1:N2 + N1] = phi_odd(2, E, x2)
			PHI[j,N1 + N2:N3 + N2 + N1] = phi_odd(3, E, x3)

		# normalitzacio de la funció d'ona
		PHI[j] /= np.sqrt(np.sum(PHI[j] * PHI[j]))

	return PHI



# In[4]:


##################################################################################
##################################################################################
##################											  	##################
##################				   PART III						##################
##################		definicio de la funcio gaussiana		##################
##################			  i aplicacio del kick				##################
##################											  	##################
##################################################################################
##################################################################################
import scipy.integrate as integrate

# Definim ara una funcio gaussiana

def gaussiana(x):
	return np.exp( - (x - xi)**2 / (4 * sigmax**2) )



# In[15]:


var('V a q')
V=0.5
a=2

f(x)=(1+V*sin(k*sqrt(x-V)*a)**2/(4*x*(x-V)))**(-1)
            


# In[18]:


f(V+1/k**2 *1/(4*sigmax**2))


# In[13]:


(1+0.5*k**2)^(-1)


# In[74]:


from __future__ import division, print_function

import numpy as np
import math
import matplotlib.pyplot as plt
import time
import os

current_milli_time = lambda: int(round(time.time() * 1000))

k = 5.12665			  	# Factor d'unitats (  eV^(-1/2)*nm^(-1)  )
hbar = 6.582e-4			# h barra (eV · ps)

L = 60.0				   	# Mitja longitud de la caixa (nm)#amb aquestes dades peta entre 26 i 27 #0.9e-6eV
xi =-L/2			   	# Posicio inicial del paquet (nm)
l = 1				   	# Mitja amplada de la barrera (nm) ho he posat a 0 per un truquis
sigmax = 6			# Incertesa inicial del paquet (nm)

T = 0.15				   	# Energia cinetica (eV)
V0 = 0.5				# Barrera de potencial (eV)

Nx = 1024				# Numero de particions en x
dx = 2 * L / Nx 	   	# Pas en x
N1, N2, N3 = int((L - l) / dx), int(2 * l/ dx), int((L - l) / dx) # Numero de particions abans, dins i despres de la barrera o del pou

dt = 0.0005			   	# Pas de temps (ps)
Nt = 400
# Numero de passos de temps

N = 2*128				# Numero d'estats propis que utilitzarem
dE = (math.pi/(2*k*L))**2/20		  		# Precisio en energies (eV)
FILE_ENERGIES = 'energies_{}.txt'.format([L,l,V0,N,dE])
FILE_PHI = 'phi_{}.txt'.format([L,l,V0,N])

x = np.linspace(-L, L, Nx)

# Evalua les energies si no estan desades. Si ja estan calculades, les llegeix del arxiu
#Ep=np.array([(n*math.pi/(2*k*L))**2 for n in range(1,N+1)])
# Evalua les energies si no estan desades. Si ja estan calculades, les llegeix del arxiu
if os.path.exists(FILE_ENERGIES):
	print('Reading energies')
	Ep = read_energies(FILE_ENERGIES)
else:
    print('Evaluating energies')
    Ep= find_roots()
    save_energies(Ep)


N=len(Ep)
start=current_milli_time()
print('Evaluating wave functions')
phi = evaluate_wave_function(Ep)/np.sqrt(dx)

print('OK (%.2f s)' % ((current_milli_time() - start) / 1000))

# La normalitzem de -L a +L
integral = integrate.quad(lambda x: (gaussiana(x))**2, -L, L)
Norm = np.sqrt((integral[0]))

# I ens definim un vector que li direm gauss, on hi posarem els valors de la gaussiana
# a cada posicio de la caixa

gauss = gaussiana(x) / Norm


# Aplicacio del kick, i separacio de la gaussiana en part real i part imaginaria

gaussr = np.cos(k*np.sqrt(T)*x)*gauss
gaussi = np.sin(k*np.sqrt(T)*x)*gauss

# Plots de la part real, la part imaginaria i el modul al quadrat de la gaussiana

coefs=np.zeros((N,Nt),dtype=np.complex_)

print('Calculant coeficients')
start=current_milli_time()
for j in range(N):

    coefs[j][0]=np.sum([phi[j][i]*gaussr[i] for i in range(Nx)])*dx+1j*np.sum([phi[j][i]*gaussi[i] for i in range(Nx)])*dx
print('OK (%.2f s)' % ((current_milli_time() - start) / 1000))    







Ep=np.array(Ep)/hbar
print('Evolucionant coeficients')
start=current_milli_time()
for i in range(N):
    for j in range(Nt):
        coefs[i][j]=coefs[i][0]*(np.cos(Ep[i]*j*dt)-1j*np.sin(Ep[i]*j*dt))
print('OK (%.2f s)' % ((current_milli_time() - start) / 1000)) 
        
print('Desfent el trencament')
start=current_milli_time()        
gauss_2=[reduce(lambda x,y:np.add(x,y),[coefs[j][ll]*phi[j] for j in range(N)]) for ll in range(Nt)]
print('OK (%.2f s)' % ((current_milli_time() - start) / 1000)) 


Paquet=np.array([abs(i) for i in gauss_2])

PROB=Paquet**2




            
    
            
            


# In[83]:


for i in range(0,10):
    y = phi[i]
    fig = plt.figure()
    ax = plt.axes(xlim=(-L, L), ylim=(-1, 1))
    line, = ax.plot([], [], lw=2)
    line.set_data(x, y)
    ax.set_title(u'Funció {}'.format(i+1))
    plt.plot()
    plt.savefig('funcio{}.eps'.format(i))


# In[142]:



import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(-L,L), ylim=(0, 1))
line, = ax.plot([], [], lw=2)
l=1
x = np.linspace(-L, L, Nx)

# initialization function: plot the background of each frame

line.set_data([], [])

    
k=2978
y = Paquet[k]
line.set_data(x, y)
ax.set_xlabel('temps {} ps'.format(float(k*dt)))
    
    


# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so 
plt.savefig('gran.eps')
        


# In[79]:



import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(-L,L), ylim=(0, 1))
line, = ax.plot([], [], lw=2)

x = np.linspace(-L, L, Nx)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    ax.set_xlabel('temps 0')
    
    
    plt.fill([-l,l,l,-l], [0,0,V0,V0], 'green', alpha=0.2)
    ax.plot(np.linspace(0,0,100),np.linspace(0,1,100),linestyle='--')

   
    return line,

# animation function.  This is called sequentially
def animate(i):
    
    y = Paquet[i]
    line.set_data(x, y)
    ax.set_xlabel('temps {} ps'.format(float(i*dt)))
    
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=(400), interval=200, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html

anim.save('patadamenuda2.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
        
        


# In[11]:


from matplotlib.animation import FFMpegWriter

class FasterFFMpegWriter(FFMpegWriter):
    '''FFMpeg-pipe writer bypassing figure.savefig.'''
    def __init__(self, **kwargs):
        '''Initialize the Writer object and sets the default frame_format.'''
        super().__init__(**kwargs)
        self.frame_format = 'argb'

    def grab_frame(self, **savefig_kwargs):
        '''Grab the image information from the figure and save as a movie frame.

        Doesn't use savefig to be faster: savefig_kwargs will be ignored.
        '''
        try:
            # re-adjust the figure size and dpi in case it has been changed by the
            # user.  We must ensure that every frame is the same size or
            # the movie will not save correctly.
            self.fig.set_size_inches(self._w, self._h)
            self.fig.set_dpi(self.dpi)
            # Draw and save the frame as an argb string to the pipe sink
            self.fig.canvas.draw()
            self._frame_sink().write(self.fig.canvas.tostring_argb()) 
        except (RuntimeError, IOError) as e:
            out, err = self._proc.communicate()
            raise IOError('Error saving animation to file (cause: {0}) '
                      'Stdout: {1} StdError: {2}. It may help to re-run '
                      'with --verbose-debug.'.format(e, out, err)) 


# In[69]:





# In[69]:





# In[197]:


l


# In[143]:


np.sinh(3)


# In[176]:


def T(V0,E,k,a):
    x=(1+V0**2*math.sin(k*a)**2/(4*E*(E-V0)))
    return 1/x


# In[181]:


E=0.5
V0=0.5

a=2
m=k**2*hbar**2/2
kk=np.sqrt(2*m*(E-V0)/hbar**2)
#T(V0,E,kk,a)
sss=(1+m*V0**2*a**2/(2*E*hbar**2))


# In[184]:


k**2*0.5*4/4


# In[80]:


mig=int(Nx/2)
R=np.array([np.sum(PROB[i][:mig]) for i in range(Nt)])
T=np.array([np.sum(PROB[i][mig:]) for i in range(Nt)])


# In[81]:


R*=dx
T*=dx


# In[82]:


tl=139

t=np.linspace(0,Nt,Nt)
fig = plt.figure()
ax = plt.axes(xlim=(0, Nt), ylim=(0, 1))
ax.set_xlabel('temps/dt')
ax.set_ylabel('probabilitat')
ax.plot(t, T, color='green')
ax.plot(t, R, color='red')
ax.plot(np.linspace(tl,tl,100),np.linspace(0,1,100),linestyle='--')
ax.text(0.9, 0.80, 'R='+'{:4.2f}'.format(R[tl]), horizontalalignment='center',
verticalalignment='center', transform=ax.transAxes,bbox=dict(edgecolor='red',facecolor='red' ,alpha=0.5))
ax.text(0.9, 0.90, 'T='+'{:4.2f}'.format(T[tl]), horizontalalignment='center',
verticalalignment='center', transform=ax.transAxes,bbox=dict(edgecolor='green',facecolor='green' ,alpha=0.5))
ax.legend()
plt.savefig('patadamenuda2.eps',bbox_inches='tight')


# In[123]:


T


# In[94]:


abs(Paquet)

