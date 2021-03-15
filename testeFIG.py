import numpy as np
import matplotlib.pyplot as plt
np.seterr(divide='ignore', invalid='ignore')


# Image formation - SIS with circular source

phi = np.arange(0.0, 2.0*np.pi, 0.0001)
xc1=np.cos(phi) #curva critica para eixo x1
xc2=np.sin(phi) #curva critica para eixo x2

#def xmenos(r0, s, e, phie):
#   return np.where(1+(1/(1-e*np.cos(2*(phi-phie))))*(s*np.cos(phi)-e*s*np.cos(phi-2*phie)-np.sqrt((r0**2)*(1-e*np.cos(2*(phi-phie)))-(s**2)*(1-e**2)*np.sin(phi)**2)))

#def xmais(r0, s, e, phie):
#    return np.where(1+(1/(1-e*np.cos(2*(phi-phie))))*(s*np.cos(phi)-e*s*np.cos(phi-2*phie)+np.sqrt((r0**2)*(1-e*np.cos(2*(phi-phie)))-(s**2)*(1-e**2)*np.sin(phi)**2)))

r0   = float(input("R0 ="))
s    = float(input("S ="))
theta   = float(input("theta ="))*(np.pi/180)
e    = float(input("e ="))
phie = float(input("phie ="))*(np.pi/180) 

xmenos = 1+(1/(1-e*np.cos(2*(phi-phie))))*(s*np.cos(theta-phi)-e*s*np.cos(theta+phi-2*phie)-np.sqrt((r0**2)*(1-e*np.cos(2*(phi-phie)))-(s**2)*(1-e**2)*np.sin(theta-phi)**2))

xmais = 1+(1/(1-e*np.cos(2*(phi-phie))))*(s*np.cos(theta-phi)-e*s*np.cos(theta+phi-2*phie)+np.sqrt((r0**2)*(1-e*np.cos(2*(phi-phie)))-(s**2)*(1-e**2)*np.sin(theta-phi)**2))

xmed = 1+(1/(1-e*np.cos(2*(phi-phie))))*(s*np.cos(theta-phi)-e*s*np.cos(theta+phi-2*phie))
#subplot(nrows,nlines,position)


#====================================================================================================
# Image Formation - Lens Plane
#====================================================================================================


ax2 = plt.subplot(122)
plt.plot(np.cos(phi), np.sin(phi),'k--',0,0,'ko')
plt.plot(xmed*np.cos(phi), xmed*np.sin(phi), 'grey', linestyle='--', lw = 1, label=r'$\bar{x}$')
#plt.plot(s*np.cos(theta)+(r0/(np.sqrt(1-e)))*np.cos(phi)*np.cos(phie)-(r0/(np.sqrt(1+e)))*np.sin(phi)*np.sin(phie),s*np.sin(theta)+(r0/(np.sqrt(1-e)))*np.cos(phi)*np.sin(phie)+(r0/(np.sqrt(1+e)))*np.sin(phi)*np.cos(phie),'purple')
titlefont = {
        'family' : 'serif',
        'color'  : 'black',
  #      'weight' : 'bold',
        'size'   : 16,
        }
plt.plot(xc1, xc2, 'k--', lw=1, label='$x_c$')
plt.plot(xmenos*np.cos(phi), xmenos*np.sin(phi), 'b-', lw=1, label='$x^{-}$')
plt.plot(xmais*np.cos(phi), xmais*np.sin(phi), 'r-', lw = 1, label='$x^{+}$') 
#plt.plot((1+(1/(1-e*np.cos(2*(phi-phie))))*(s*np.cos(phi)-e*s*np.cos(phi-2*phie)))*np.cos(phi), (1+(1/(1-e*np.cos(2*(phi-phie))))*(s*np.cos(phi)-e*s*np.cos(phi-2*phie)))*np.sin(phi) ,  'k--', lw=0.5, label=r'$\bar{x}$')
titlefont = {
        'family' : 'serif',
        'color'  : 'black',
  #      'weight' : 'bold',
        'size'   : 16,
        }

ax2.set_title("Plano da Lente", # title
             va='bottom',                   # some space below the title
             fontdict = titlefont           # set the font properties
             )        
ax2.grid(False, which='both')
ax2.axis([-2,2, -2, 2])
ax2.axhline(y=0, color='0.9',ls='--',lw=1)
ax2.axvline(x=0, color='0.9',ls='--',lw=1)
plt.xlabel('$x_1$',fontsize=20)
plt.ylabel('$x_2$',fontsize=20)
ax2.legend(loc="lower left")           # legend location

plt.show() 

