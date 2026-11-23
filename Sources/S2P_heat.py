import numpy as np
import matplotlib.pyplot as plt

class Strata2P_sim:
    def __init__(self, simulant = 'JSC1', r_heater=.3175, r_chamber=2,dr=.01,t=3600,dt=1,debug=False,
                inner_bc_type="Type 1",
                 inner_bc_value= 250.,
                 outer_bc_type="Type 2",
                 outer_bc_value= 0.):
        #Bring in values:
        self.debug=debug
        self.simulant=simulant
        self.r_heater=r_heater
        self.r_chamber=r_chamber
        self.t=t
        self.dt=dt
        self.dr=dr
        self.inner_bc_type=inner_bc_type
        self.inner_bc_value=inner_bc_value
        self.outer_bc_type=outer_bc_type
        self.outer_bc_value=outer_bc_value
        # new for the control_function also housekeeping just in case we use control_function later
        # self.target_T=target_T
        # self.heater_power=heater_power
        # target_T=250,heater_power=40

        #Bring in simulant parameters:
        if (self.simulant=='JSC1'):
            self.D = 1e-7
            self.K = .15
        
        if (self.simulant=='sand'):
            self.D = 1e-5
            self.K = 1.5
        
        self.factor = self.D * self.dt / (self.dr ** 2)
        print(f'Stability factor: {self.factor} (must be <= 0.5)')

        
        #Construct required objects:
        self.Time = np.arange(0,self.t,self.dt)
        self.Rad  = np.arange(self.r_heater,self.r_chamber,self.dr)
        self.Temp = np.zeros((len(self.Rad), len(self.Time)))
        
    def details(self):
        print('debug   : {}'.format(self.debug))
        print('simulant: {}'.format(self.simulant))
        print('heater d: {} cm'.format(2*self.r_heater))
        print('radius  : {} cm'.format(self.r_chamber))
        print('dr      : {} cm'.format(self.dr))
        print('time lim: {} s'.format(self.t))
        print('dt      : {} cm'.format(self.dt))

    def run(self):
        #Step through TIME j:
        if(self.debug):
            print('Begin Running Simulation')
        for j in range(0,len(self.Time)-1):
            self.experiment_control()
            self.handle_boundaries(j)
            self.time_step(j) #Use timestep Temp[:,j] to compute Temp[:,j+1] for Temp[r,t]
        return

    #Step through SPACE i at this TIME j:
    def time_step(self,j):
        if(self.debug):
            print(f'timestep {j}')
        for i in range(1, len(self.Rad) - 1):
            self.Temp[i, j+1] = self.Temp[i, j] + self.factor * (self.Temp[i-1, j] - 2*self.Temp[i, j] + self.Temp[i+1, j]) + 1/(2*self.Rad[i])*self.D*self.dt/self.dr * (self.Temp[i+1,j] - self.Temp[i-1,j])
        return

    #Control experiment by setting appropriate BC RULES. For example,
    #  once heater temp gets to a set_temp, switch from dT/dr BC to
    #  a constant T BC (neumann->dirachlet)
    def experiment_control(self):
        if(self.debug):
            print('Experiment Control')
        return #TODO
        
    #Each timestep, implement the BC RULES (as updated by expriment_control())
    def handle_boundaries(self,j):
        if(self.debug):
            print('Handle Boundaries')

        if(self.inner_bc_type == "Type 1"):
            self.Temp[0,j] = self.inner_bc_value
        elif(self.inner_bc_type == "Type 2"):
            self.Temp[0,j] = self.Temp[1,j]+self.inner_bc_value
        else:
            print("Invalid inner BC type")
            return 1
            
        if(self.outer_bc_type == "Type 1"):
            self.Temp[-1,j] = self.outer_bc_value
        elif(self.outer_bc_type == "Type 2"):
            self.Temp[-1,j] = self.Temp[-2,j]+self.outer_bc_value
        else:
            print("Invalid outer BC type")
            return 1
        
        return 0

    def plot(self):
        #Plot the heat:
        fig,ax=plt.subplots()
        Rax, Tax = np.meshgrid(self.Rad, self.Time) #<- sets up axis grids for imshow
        p=ax.imshow(np.transpose(self.Temp),aspect='auto',origin='lower',extent=[self.Rad[0]*100, self.Rad[-1]*100, 0, self.Time[-1]],vmin=0,vmax=None)
        fig.colorbar(p,ax=ax,label=r'Temp ($\degree C$)')
        ax.set_ylabel('Time (s)')
        ax.set_xlabel('r (cm)')
        